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
%  死区 (opts.Tdead, 例如 4e-6 = 4us):
%     开通延迟 td / 关断立即; 死区期间桥臂电位由续流二极管按电流方向钳位。
%     每个载波周期的平均误差电压为 dV = -2*sign(i)*td*fc*Vdc (两个桥臂都斩波)。
%     后果: 原本严格为 0 的 3、5、7... 次谐波出现(幅值 ~ 8*td*fc*Vdc/(pi*h)),
%           基波幅值下降(约为 dV 的基波分量), 窄脉冲(< td)丢失。
%     电流方向由"电压谐波 -> R-L 电流"自洽迭代确定; 无负载信息时按基波
%     sin(2*pi*f0*t - phi0) 的符号近似(phi0 = atan(2*pi*f0*L/R))。
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
%  用法(死区):
%      R = hbridge_spwm_harmonics(2400, 50, 100, 0.8, Tdead=4e-6, Rload=0.5, Lload=2e-3);
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
    opts.Lload (1,1) double {mustBeNonnegative} = 0                % 负载电感 H
    opts.Tdead (1,1) double {mustBeNonnegative} = 0                % 固定开关死区 s (如 4e-6)
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
R.Tdead = opts.Tdead;
R.td_over_Tc = opts.Tdead * fc;        % 死区占载波周期的比例
R.td_periods = opts.Tdead * fc;        % 每开关周期损失的有效导通比例

% ================= 1) 精确法 =================
if opts.Tdead > 0
    % 含死区: 断点=各桥臂跳变时刻 ∪ (跳变+td) ∪ 电流过零点; 电流符号自洽迭代
    phi0 = atan2(2*pi*f0*opts.Lload, max(opts.Rload, eps));
    [tr, segA, segB, segLev, dtinfo] = ...
        exact_waveform_deadtime(fc, f0, M, Vdc, T, opts.Tdead, phi0, ...
                                opts.Rload, opts.Lload, opts.Kmax);
    R.deadtime = dtinfo;
else
    [tr, segA, segB, segLev] = exact_waveform(fc, f0, M, Vdc, T);
    R.deadtime = struct('converged', true, 'nIter', 0, 'iterated', false, ...
                        'phi_deg', NaN);
end
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
R.A1_theory = M * Vdc;                       % 线性调制区(无死区): A1 = M*Vdc
R.A1_error  = abs(R.Amp(1) - R.A1_theory);
R.A1_drop   = R.A1_theory - R.Amp(1);        % 死区造成的基波损失
R.A1_drop_pct = 100 * R.A1_drop / R.A1_theory;

% ---- 死区引入的低次谐波 (理想时严格为 0) ----
if opts.Tdead > 0
    R.deadtime.dV_per_Tc = 2 * opts.Tdead * fc * Vdc;    % 每载波周期平均误差幅值
    R.deadtime.low_order_est = zeros(0,2);
    lo = 3:2:min(25, opts.Kmax);
    est = 8 * opts.Tdead * fc * Vdc ./ (pi * lo).^1;      % 8*td*fc*Vdc/(pi*h)
    R.deadtime.low_order_est = [lo(:), est(:), R.Amp(lo)];
end

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

% ---- 死区"电流方向由基波决定"假设的量化核对 ----
if opts.Tdead > 0
    if opts.Rload > 0
        tg = linspace(0, T, 20001).';            % 固定网格, 便于跨工况比较
        iFull = local_current_eval(tg, ck, opts.Rload, opts.Lload, 2*pi*f0);
        upFull = iFull > 0;
        upFund = sin(2*pi*f0*tg - R.deadtime.phi_rad) > 0;
        R.deadtime.sign_mismatch_pct = 100 * mean(upFull ~= upFund);
    else
        R.deadtime.sign_mismatch_pct = NaN;
    end
    R.deadtime.phi_rad_final = R.deadtime.phi_rad;
end

% ================= 2) 采样FFT法 (复现原Python仿真) =================
if opts.UseSampledFFT
    N  = round(T / opts.DT);
    dt = T / N;                        % 保证恰好 Ncyc 个整基波周期
    t  = (0:N-1).' * dt;
    if opts.Tdead > 0
        Vout = Vdc * spwm_deadtime_out_level(t, fc, f0, M, opts.Tdead, ...
                                             R.deadtime.phi_rad);
    else
        Vout = Vdc * spwm_out_level(t, fc, f0, M);
    end

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
%  含死区的精确波形 (分段常数)
%    死区模型: 开通延迟 td / 关断立即, 死区期间由续流二极管按电流方向钳位
%    断点集合 = 各桥臂理想跳变时刻 ∪ (跳变时刻+td) ∪ 电流过零点
%
%    电流方向的自洽条件: 续流钳位取决于电流符号, 而电流又由(含死区的)电压
%    决定。由于 R-L 负载把开关纹波滤得很干净, 电流过零点几乎完全由**基波**
%    决定, 故以基波滞后相位 phi 为自洽变量做标量固定点迭代(收敛快且良态):
%        phi = angle(Z1) - angle(ck_1) - pi/2
%    其中 angle(ck_1) 是含死区后输出电压基波的相位。(开关纹波确实会在基波
%    过零附近使电流在载波周期内来回变号, 这属于"过零畸变"区域; 本模型用基波
%    方向表示该区域, 并用 sign_mismatch 指标量化该近似的误差。)
% ------------------------------------------------------------------
function [tr, segA, segB, segLev, info] = exact_waveform_deadtime( ...
        fc, f0, M, Vdc, T, td, phi0, Rload, Lload, Kmax) %#ok<INUSD>
t1 = -2*td;  t2 = T;
[eAll, eL, eR] = spwm_ideal_edges(fc, f0, M, t1, t2);

inT = @(x) x(x >= 0 & x <= T);
brkBase = unique([0; inT(eAll); inT(eL) + td; inT(eR) + td; T]);
brkBase = brkBase(brkBase >= 0 & brkBase <= T);

Z1ang    = angle(Rload + 1j*2*pi*f0*Lload);      % 基波阻抗角
w        = 2*pi*f0;
phi      = phi0;
iterated = Rload > 0;
converged = ~iterated;
nIter = 0;
maxIter = 15;
tol = 1e-9;                                      % rad

for it = 1:maxIter
    nIter = it;
    brk = unique([brkBase; local_crossings(phi, w, T)]);
    brk = sort(brk(:));
    a = brk(1:end-1);  b = brk(2:end);
    keep = (b - a) > 0;
    a = a(keep);  b = b(keep);
    if isempty(a)
        a = 0; b = T;
    end
    lev = Vdc * spwm_deadtime_out_level((a + b)/2, fc, f0, M, td, phi);

    if ~iterated
        break                                    % 无负载信息: 用基波符号, 不迭代
    end
    % 只需基波相位即可更新 phi (不必算全谱, 因此迭代代价极小)
    ck1 = exact_fourier(a, b, lev, T, f0, 1);
    phiNew = Z1ang - angle(ck1(1)) - pi/2;
    d = mod(phiNew - phi + pi, 2*pi) - pi;       % 归一化到 (-pi, pi]
    phi = phi + 0.8*d;                           % 松弛迭代
    if abs(d) < tol
        converged = true; break
    end
end

segA = a; segB = b; segLev = lev;
tr = unique([0; segA(:); segB(:); T]);

info = struct();
info.converged   = converged;
info.iterated    = iterated;
info.nIter       = nIter;
info.td          = td;
info.phi_rad     = phi;
info.phi_deg     = phi * 180/pi;
info.phi0_deg    = angle(Rload + 1j*2*pi*f0*Lload) * 180/pi;
info.nEdges      = numel(eAll);
info.nSegments   = numel(segA);
end

% ------------------------------------------------------------------
%  基波电流的过零点 (解析):  sin(2*pi*f0*t - phi) = 0
% ------------------------------------------------------------------
function tz = local_crossings(phi, w, T)
k1 = ceil((-phi)/pi);
k2 = floor((w*T - phi)/pi);
tz = ((k1:k2).' * pi + phi) / w;
tz = tz(tz > 0 & tz < T);
end

% ------------------------------------------------------------------
%  由电压谐波系数重建 R-L 负载电流 (相量法), 用于核对电流方向假设
%     v_h(t) = 2*Re(ck_h e^{j h w t});  i_h = v_h / (R + j h w L)
% ------------------------------------------------------------------
function i = local_current_eval(t, ck, R, L, w)
Ih = ck ./ (R + 1j*(1:numel(ck)).'*w*L);
t  = t(:).';
i  = zeros(size(t));
hh = (1:numel(Ih)).';
blk = 400;                                       % 分块, 控制内存
for i0 = 1:blk:numel(Ih)
    idx = i0:min(i0+blk-1, numel(Ih));
    i = i + (2*Ih(idx)).' * exp(1j*(w*hh(idx))*t);   % (1xn)*(nxNt) = 1xNt
end
i = i(:);
end

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
if R.Tdead > 0
    fprintf('  开关死区     td  : %g us  (占载波周期 %.2f %%)\n', ...
        R.Tdead*1e6, 100*R.td_over_Tc);
end
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
if R.Tdead > 0
    fprintf('  ---- 死区影响 (td = %g us) ----\n', R.Tdead*1e6);
    fprintf('  每载波周期平均误差电压        : %.4f V  (= 2*td*fc*Vdc)\n', ...
        R.deadtime.dV_per_Tc);
    fprintf('  基波损失                      : %.4f V (%.4f %%)\n', ...
        R.A1_drop, R.A1_drop_pct);
    if R.deadtime.iterated
        fprintf('  电流方向自洽迭代(基波相位)    : %d 次, %s (phi = %.3f deg)\n', ...
            R.deadtime.nIter, ternary(R.deadtime.converged, '已收敛', '未收敛(!)'), ...
            R.deadtime.phi_deg);
        if isfinite(R.deadtime.sign_mismatch_pct)
            fprintf('  电流符号假设偏差(全谐波核对)  : %.4f %% (纹波在过零附近改变方向)\n', ...
                R.deadtime.sign_mismatch_pct);
        end
    else
        fprintf('  电流符号                      : 未迭代(按基波近似, 建议给出 Rload/Lload)\n');
    end
    fprintf('  低次谐波 (理论应为 0) 实测/估计:\n');
    fprintf('     h      幅值/V     8*td*fc*Vdc/(pi*h)/V\n');
    for k = 1:size(R.deadtime.low_order_est,1)
        fprintf('   %4d  %11.6f  %11.6f\n', R.deadtime.low_order_est(k,1), ...
            R.deadtime.low_order_est(k,3), R.deadtime.low_order_est(k,2));
    end
end
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
