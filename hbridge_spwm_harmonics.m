function R = hbridge_spwm_harmonics(fc, f0, Vdc, M, opts)
%HBRIDGE_SPWM_HARMONICS  H桥单极性SPWM输出电压谐波分析(精确法 + 采样FFT法)
%
%  调制规律与原 Python 脚本 h_bridge_svpwm.py 完全一致:
%      phase = frac(t*fc)
%      c(t)  = 2*|2*phase - 1| - 1          三角载波, 峰值 +1 / 谷值 -1
%      v(t)  = M*sin(2*pi*f0*t)             调制波
%      左桥臂 Q1 = ( v > c)     右桥臂 Q3 = (-v > c)
%      Q2 = 1-Q1, Q4 = 1-Q3
%      Vout  = Vdc*(Q1 - Q3) = Vdc * sign(v) .* (|v| > |c|)
%  即三电平输出 {+Vdc, 0, -Vdc}, 零电平段为续流。
%
%  谐波计算提供两种互相独立、可交叉验证的方法:
%    (1) 精确法 : 解析求出每个开关跳变时刻(二分法收敛到机器精度),
%                 再把分段常数波形逐段解析积分得到傅里叶系数。
%                 不含任何时间网格量化误差, 是谐波幅值的"真值"。
%    (2) 采样FFT法: 与原 Python 脚本相同的等步长采样(DT, 默认1us)+FFT,
%                 用于评估原仿真步长带来的谐波误差。
%
%  用法:
%      R = hbridge_spwm_harmonics();                          % 默认 fc=2400, f0=50, Vdc=100, M=0.8
%      R = hbridge_spwm_harmonics(2400, 50, 100, 0.8);
%      R = hbridge_spwm_harmonics(2400, 50, 100, 0.8, Ncyc=10, Kmax=4000, Rload=0.5, Lload=2e-3);
%
%  返回结构体 R 主要字段:
%      R.harm, R.freq       谐波次数 / 频率(Hz), 下标 h 对应 h 次谐波
%      R.Amp                各次谐波幅值(峰值,V), 精确法
%      R.Ph_deg             各次谐波相位(度)
%      R.DC                 直流分量(V)
%      R.THD_50/R.THD_2fc/R.THD_all   不同求和上限的电压THD(相对基波有效值)
%      R.WTHD               1/h 加权THD(更接近电流畸变程度)
%      R.RMS_time/R.RMS_spec 波形有效值 / 频谱有效值(帕塞瓦尔校验)
%      R.Amp_fft, R.err_rel_fft  1us采样FFT的幅值及其相对精确法的误差
%      R.Iamp, R.THD_i_*    负载电流各次谐波幅值(A) 及其THD (需 Rload>0)
%      R.exact.transitions  精确开关跳变时刻(s)
%
%  说明: 当 fc/f0 为偶数整数时, 输出波形具有半波对称性, 只存在奇次谐波。
%
%  作者: 谐波分析工具 (MATLAB R2025b)
%  日期: 2026
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    fc  (1,1) double {mustBePositive} = 2400      % 载波(开关)频率 Hz
    f0  (1,1) double {mustBePositive} = 50        % 基波(调制波)频率 Hz
    Vdc (1,1) double {mustBePositive} = 100       % 直流母线电压 V
    M   (1,1) double {mustBeNonnegative} = 0.8    % 调制比
    opts.Ncyc  (1,1) double {mustBePositive, mustBeInteger} = 10  % 分析周期数
    opts.DT    (1,1) double {mustBePositive} = 1e-6               % 采样步长 s
    opts.Kmax  (1,1) double {mustBePositive, mustBeInteger} = 4000 % 最高谐波次数
    opts.Rload (1,1) double {mustBeNonnegative} = 0               % 负载电阻 ohm
    opts.Lload (1,1) double {mustBeNonnegative} = 0               % 负载电感 H
    opts.UseSampledFFT (1,1) logical = true
    opts.Verbose       (1,1) logical = true
end

T     = opts.Ncyc / f0;
ratio = fc / f0;

R = struct();
R.fc = fc; R.f0 = f0; R.Vdc = Vdc; R.M = M;
R.T = T; R.Ncyc = opts.Ncyc; R.Kmax = opts.Kmax; R.DT = opts.DT;
R.carrierRatio  = ratio;
R.isIntegerRatio = abs(ratio - round(ratio)) < 1e-9;
R.isEvenRatio    = R.isIntegerRatio && mod(round(ratio), 2) == 0;
R.isPeriodicRecord = abs(opts.Ncyc*ratio - round(opts.Ncyc*ratio)) < 1e-9;
if ~R.isPeriodicRecord
    warning('hbridge_spwm:harmonics', ...
        ['记录长度 %d 个基波周期不是载波周期的整数倍 (Ncyc*fc/f0 = %g),\n' ...
         '        频谱会出现泄漏, 请令 fc 为 f0 的整数倍, 或令 Ncyc*fc/f0 为整数。'], ...
        opts.Ncyc, opts.Ncyc*ratio);
end
R.Rload = opts.Rload; R.Lload = opts.Lload;

% ================= 1) 精确法 =================
[tr, segA, segB, segLev] = exact_waveform(fc, f0, M, Vdc, T);
[ck, ck0] = exact_fourier(segA, segB, segLev, T, f0, opts.Kmax);

R.exact.transitions = tr;
R.exact.nTransitions = numel(tr);
R.exact.nSegments    = numel(segLev);
R.exact.segA = segA; R.exact.segB = segB; R.exact.segLev = segLev;
R.exact.DC   = real(ck0);
R.exact.RMS  = sqrt(sum(segLev.^2 .* (segB - segA)) / T);
R.RMS_time   = R.exact.RMS;   % 波形有效值(精确)

R.harm   = (1:opts.Kmax).';
R.freq   = R.harm * f0;
R.Amp    = 2 * abs(ck);            % 峰值幅值
R.Ph_deg = angle(ck) * 180 / pi;   % 相位(度), 余弦形式
R.DC     = real(ck0);
R.RMS_spec = sqrt(R.DC^2 + sum(R.Amp.^2) / 2);

% ---- 畸变指标 ----
[R.THD_50,  R.RMS_h_50]  = local_thd(R.Amp, [2, min(50, opts.Kmax)]);
[R.THD_2fc, R.RMS_h_2fc] = local_thd(R.Amp, [2, min(floor(2*fc/f0), opts.Kmax)]);
[R.THD_4fc, R.RMS_h_4fc] = local_thd(R.Amp, [2, min(floor(4*fc/f0), opts.Kmax)]);
[R.THD_all, R.RMS_h_all] = local_thd(R.Amp, [2, opts.Kmax]);
if R.Amp(1) > 0
    R.WTHD = sqrt(sum((R.Amp(2:end) ./ R.harm(2:end)).^2)/2) / (R.Amp(1)/sqrt(2));
else
    R.WTHD = NaN;
end

% ---- 基波理论值 ----
R.A1_theory = M * Vdc;                       % 线性调制区: A1 = M*Vdc
R.A1_error  = abs(R.Amp(1) - R.A1_theory);

% ---- 负载电流谐波 (把输出看成经 R-L 负载的电压源) ----
if opts.Rload > 0
    R.Zmag = abs(opts.Rload + 1j*2*pi*f0*R.harm*opts.Lload);
    R.Iamp = R.Amp ./ R.Zmag;
    R.THD_i_50  = local_thd(R.Iamp, [2, min(50, opts.Kmax)]);
    R.THD_i_2fc = local_thd(R.Iamp, [2, min(floor(2*fc/f0), opts.Kmax)]);
    R.THD_i_all = local_thd(R.Iamp, [2, opts.Kmax]);
else
    R.Zmag = []; R.Iamp = [];
    R.THD_i_50 = NaN; R.THD_i_2fc = NaN; R.THD_i_all = NaN;
end

% ================= 2) 采样FFT法 (复现原Python仿真) =================
if opts.UseSampledFFT
    N  = round(T / opts.DT);
    dt = T / N;                        % 保证恰好 Ncyc 个整基波周期
    t  = (0:N-1).' * dt;
    Vout = Vdc * spwm_out_level(t, fc, f0, M);

    Y    = fft(Vout);
    half = floor(N/2);
    mag  = 2*abs(Y(1:half+1)) / N;
    mag(1) = abs(Y(1)) / N;            % 直流

    R.sampled.t = t; R.sampled.Vout = Vout; R.sampled.dt = dt; R.sampled.N = N;
    R.sampled.f = (0:half).' / T;
    R.sampled.mag = mag;
    R.sampled.RMS = sqrt(mean(Vout.^2));

    Af  = nan(opts.Kmax, 1);
    idx = R.harm * opts.Ncyc + 1;      % 第h次谐波落在 bin = h*Ncyc
    ok  = idx <= half + 1;
    Af(ok) = mag(idx(ok));
    R.Amp_fft   = Af;
    R.err_rel_fft = abs(Af - R.Amp) ./ max(R.Amp, 1e-12);
    R.Amp_fft_err_abs = abs(Af - R.Amp);

    % 只把幅值 > 1% Vdc 的谐波作为"显著谐波"统计误差 (微小谐波的相对误差无意义)
    sig = R.Amp > 0.01 * Vdc;
    R.sig_mask = sig;
    if any(sig)
        [R.Amp_fft_maxerr, kk] = max(R.err_rel_fft(sig));
        hs = find(sig);
        R.Amp_fft_maxerr_h = hs(kk);
    else
        R.Amp_fft_maxerr = 0; R.Amp_fft_maxerr_h = NaN;
    end
    R.Amp_fft_maxerr_abs = max(R.Amp_fft_err_abs);
    R.Amp_fft_fund_err   = abs(Af(1) - R.Amp(1)) / R.Amp(1);
    R.Amp_fft_cluster_err = NaN;
    if R.isIntegerRatio
        hc = 2 * round(ratio);                       % 2fc 附近的谐波群
        hwin = max(1, hc-7) : min(numel(Af), hc+7);
        m2 = sig(hwin);                              % 只看其中显著的谐波
        if any(m2)
            R.Amp_fft_cluster_err = max(R.err_rel_fft(hwin(m2)));
        end
    end
else
    R.Amp_fft = []; R.err_rel_fft = [];
    R.Amp_fft_maxerr = NaN; R.sampled = struct();
end

if opts.Verbose
    local_verbose(R);
end
end % ====================== 主函数结束 ======================

% ------------------------------------------------------------------
%  精确波形: 求跳变时刻 + 分段常数表示
% ------------------------------------------------------------------
function [tr, segA, segB, segLev] = exact_waveform(fc, f0, M, Vdc, T)
tdummy = [0; T];
if M > 0
    dtScan = max(1/(1000*fc), T/4e6);        % 扫描步长 <= 载波周期/1000
    ts = (0:dtScan:T).';
    if ts(end) < T - eps(T)
        ts(end+1) = T;
    end
    lv  = spwm_out_level(ts, fc, f0, M);
    chg = find(lv(1:end-1) ~= lv(2:end));
    if ~isempty(chg)
        a0 = ts(chg); b0 = ts(chg+1);
        tr2 = refine_roots(@(t) spwm_level_residual(t, fc, f0, M), a0, b0);
        tdummy = [0; tr2(:); T];
    end
end
tr = unique(tdummy);
tr = tr([true; diff(tr) > 1e-13*T]);          % 去掉重合点

segA = tr(1:end-1); segB = tr(2:end);
keep = (segB - segA) > 0;
segA = segA(keep); segB = segB(keep);
segLev = Vdc * spwm_out_level((segA + segB)/2, fc, f0, M);
end

% ------------------------------------------------------------------
%  跳变时刻细化: 向量化二分法 (要求区间端点异号)
% ------------------------------------------------------------------
function tr = refine_roots(F, a, b)
fa = F(a); fb = F(b);
ok = (fa .* fb) < 0;
tr = 0.5*(a + b);
if any(ok)
    lo = a(ok); hi = b(ok); flo = fa(ok);
    for it = 1:80
        mid  = 0.5*(lo + hi);
        fm   = F(mid);
        same = (flo .* fm) > 0;
        lo(same) = mid(same); flo(same) = fm(same);
        hi(~same) = mid(~same);
    end
    tr(ok) = 0.5*(lo + hi);
end
% 端点同号的罕见情况(跳变恰落在扫描格点上或相切): 用局部极小定位
for k = find(~ok).'
    tr(k) = fminbnd(@(t) abs(F(t)), a(k), b(k));
end
end

% ------------------------------------------------------------------
%  逐段解析积分求傅里叶系数
%      c_h = (1/T) * sum_j lev_j * dt_j * exp(-1j*pi*h*f0*(a_j+b_j)) * sinc(h*f0*dt_j)
% ------------------------------------------------------------------
function [ck, ck0] = exact_fourier(a, b, lev, T, f0, Kmax)
dt  = b - a;
mid = a + b;
ck  = zeros(Kmax, 1);
blk = 500;                                   % 分块, 控制内存
for i0 = 1:blk:Kmax
    idx = i0:min(i0+blk-1, Kmax);
    hh  = idx(:);
    E = exp(-1j*pi*(hh*f0)*mid.');           % numel(idx) x Nseg
    S = sinc((hh*f0)*dt.');
    ck(idx) = (E.*S) * (lev.*dt) / T;
end
ck0 = sum(lev.*dt) / T;
end

% ------------------------------------------------------------------
%  畸变指标
% ------------------------------------------------------------------
function [thd, rms_h] = local_thd(Amp, ordRange)
if numel(Amp) < 2 || Amp(1) <= 0
    thd = NaN; rms_h = NaN; return
end
h1 = max(2, ordRange(1));
h2 = min(ordRange(2), numel(Amp));
if h2 < h1
    thd = 0; rms_h = 0; return
end
rms_h = sqrt(sum(Amp(h1:h2).^2) / 2);
thd   = rms_h / (Amp(1)/sqrt(2));
end

% ------------------------------------------------------------------
function local_verbose(R)
fprintf('\n');
fprintf('=========== H桥单极性SPWM 输出电压谐波分析 (精确解析法) ===========\n');
fprintf('  载波(开关)频率 fc : %g Hz\n', R.fc);
fprintf('  基波频率     f0  : %g Hz\n', R.f0);
fprintf('  载波比     fc/f0 : %g %s\n', R.carrierRatio, ...
    ternary(R.isIntegerRatio, '(整数: 仅奇次谐波, 纯正弦级数)', '(非整数!)'));
if ~R.isPeriodicRecord
    fprintf('  [警告] 记录长度不是载波周期的整数倍, 频谱存在泄漏!\n');
end
fprintf('  调制比       M   : %g\n', R.M);
fprintf('  直流母线   Vdc   : %g V\n', R.Vdc);
fprintf('  分析时长         : %g ms (%d 个基波周期, %g 个开关周期)\n', ...
    R.T*1000, R.Ncyc, R.T*R.fc);
fprintf('  精确跳变时刻数   : %d\n', R.exact.nTransitions);
fprintf('  基波幅值 A1      : %.6f V  (理论 M*Vdc = %.6f V, 误差 %.2e V)\n', ...
    R.Amp(1), R.A1_theory, R.A1_error);
fprintf('  直流分量         : %.3e V\n', R.DC);
fprintf('  输出有效值       : %.4f V (波形) / %.4f V (频谱, <=%g kHz)\n', ...
    R.RMS_time, R.RMS_spec, R.Kmax*R.f0/1000);
fprintf('  ---- 电压谐波畸变(相对基波有效值) ----\n');
fprintf('  THD (2~%d次)             : %8.3f %%\n', min(50,R.Kmax), 100*R.THD_50);
fprintf('  THD (2~%d次, <=2fc)       : %8.3f %%\n', min(floor(2*R.fc/R.f0),R.Kmax), 100*R.THD_2fc);
fprintf('  THD (2~%d次, <=4fc)       : %8.3f %%\n', min(floor(4*R.fc/R.f0),R.Kmax), 100*R.THD_4fc);
fprintf('  THD (2~%d次, <=%g kHz)   : %8.3f %%\n', R.Kmax, R.Kmax*R.f0/1000, 100*R.THD_all);
fprintf('  WTHD (1/h加权)           : %8.3f %%\n', 100*R.WTHD);
if R.Rload > 0
    fprintf('  ---- 负载电流谐波 (R=%g ohm, L=%g H) ----\n', R.Rload, R.Lload);
    fprintf('  电流THD (2~50次)         : %8.3f %%\n', 100*R.THD_i_50);
    fprintf('  电流THD (2~%d次, <=2fc)   : %8.3f %%\n', ...
        min(floor(2*R.fc/R.f0),R.Kmax), 100*R.THD_i_2fc);
    fprintf('  电流THD (2~%d次)         : %8.3f %%\n', R.Kmax, 100*R.THD_i_all);
end
if isfield(R,'Amp_fft') && ~isempty(R.Amp_fft)
    fprintf('  ---- 与原1us等步长采样FFT对比 (评估原Python仿真精度) ----\n');
    fprintf('  基波相对误差                : %.4f %%\n', 100*R.Amp_fft_fund_err);
    fprintf('  显著谐波(>1%%Vdc)最大相对误差  : %.3f %% (h=%d)\n', ...
        100*R.Amp_fft_maxerr, R.Amp_fft_maxerr_h);
    if isfinite(R.Amp_fft_cluster_err)
        fprintf('  2fc谐波群最大相对误差        : %.3f %%\n', 100*R.Amp_fft_cluster_err);
    end
    fprintf('  最大绝对误差                : %.4f V (%.3f %% Vdc)\n', ...
        R.Amp_fft_maxerr_abs, 100*R.Amp_fft_maxerr_abs/R.Vdc);
end
fprintf('==================================================================\n');
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
