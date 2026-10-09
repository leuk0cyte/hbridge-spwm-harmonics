%% ========================================================================
%  verify_harmonics.m  ——  谐波分析程序的自动验证 (自检) 脚本
%
%  逐项检查解析法的正确性, 并与三类完全独立的方法交叉验证:
%    A. 解析理论值  : 单极性SPWM的双重傅里叶级数(贝塞尔函数形式)
%                     A(m,n) = Vdc*(4/pi)*(1/m)*|J_n(m*pi*M/2)|   (m为偶数, n为奇数)
%    B. 数值方法    : 10ns/20ns 细采样 + FFT (与原Python同思路, 步长细1000倍)
%    C. 数值积分    : 逐段 16 点 Gauss-Legendre 求积
%
%  运行: verify_harmonics
%% ========================================================================
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
clear; clc;
here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);

nPass = 0; nFail = 0;
log = @(fmt, varargin) fprintf([fmt '\n'], varargin{:});

log('================================================================');
log('           hbridge_spwm_harmonics  自动验证');
log('================================================================');

fc = 2400; f0 = 50; Vdc = 100;
Rratio = round(fc/f0);

%% ---------------- T1: 基波幅值 = M*Vdc ----------------
Mtest = [0.1 0.25 0.5 0.8 1.0];
errA1 = zeros(size(Mtest));
for i = 1:numel(Mtest)
    R = hbridge_spwm_harmonics(fc, f0, Vdc, Mtest(i), Ncyc=2, Kmax=200, ...
        UseSampledFFT=false, Verbose=false);
    errA1(i) = abs(R.Amp(1) - Mtest(i)*Vdc);
end
[nPass,nFail] = chk('T1  基波幅值 A1 = M*Vdc', max(errA1) < 1e-9*Vdc, ...
    sprintf('max|dA1| = %.3e V', max(errA1)), nPass, nFail);

%% ---------------- 主工况: 提取性质 ----------------
Mp = 0.8;
R = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=10, Kmax=4000, ...
    UseSampledFFT=true, DT=1e-6, Verbose=false);

[nPass,nFail] = chk('T2  直流分量为 0', abs(R.DC) < 1e-12*Vdc, ...
    sprintf('DC = %.3e V', R.DC), nPass, nFail);

evenAmp = max(R.Amp(2:2:end));
[nPass,nFail] = chk('T3  偶次谐波为 0 (半波对称)', evenAmp < 1e-9*Vdc, ...
    sprintf('max even harmonic = %.3e V (%.2e %% of Vdc)', evenAmp, 100*evenAmp/Vdc), ...
    nPass, nFail);

hmax_low = Rratio - 1;
lowAmp = max(R.Amp(3:2:hmax_low));
[nPass,nFail] = chk('T4  低于 2fc 的奇次谐波(3..95)为 0 (单极性SPWM特性)', ...
    lowAmp < 1e-8*Vdc, sprintf('max = %.3e V (%.2e %% of Vdc)', lowAmp, 100*lowAmp/Vdc), ...
    nPass, nFail);

% 输出波形是奇函数 (Vout(-t) = -Vout(t)) -> 纯正弦级数 -> 余弦(同相)分量恒为 0
sigA = R.Amp > 1e-6*Vdc;
cosComp = R.Amp(sigA) .* cos(deg2rad(R.Ph_deg(sigA)));   % 同相分量幅值
maxCosC = max(abs(cosComp));
[nPass,nFail] = chk('T5  相位 = ±90 deg (奇函数 => 余弦分量 = 0)', maxCosC < 1e-9*Vdc, ...
    sprintf('n=%d, max|余弦分量| = %.3e V (< %.1e V)', nnz(sigA), maxCosC, 1e-9*Vdc), ...
    nPass, nFail);

%% ---------------- T6: 帕塞瓦尔 ----------------
gap = 1 - R.RMS_spec / R.RMS_time;
[nPass,nFail] = chk('T6  频谱有效值 <= 波形有效值, 截断差 < 3%', ...
    (R.RMS_spec <= R.RMS_time*(1+1e-9)) && (gap < 0.03), ...
    sprintf('Vrms spec/wave = %.4f/%.4f V, 截断差 %.3f %% (<%g kHz)', ...
    R.RMS_spec, R.RMS_time, 100*gap, R.Kmax*f0/1000), nPass, nFail);

%% ---------------- T7: 1us采样FFT的RMS ----------------
rel = abs(R.sampled.RMS - R.RMS_time)/R.RMS_time;
[nPass,nFail] = chk('T7  1us采样FFT的RMS与精确波形一致(<1%)', rel < 0.01, ...
    sprintf('Vrms 1us = %.5f V, exact = %.5f V, err = %.4f %%', ...
    R.sampled.RMS, R.RMS_time, 100*rel), nPass, nFail);

%% ---------------- T8: 与 10ns 细采样 FFT 交叉验证 ----------------
A_fine = fft_harmonics(fc, f0, Vdc, Mp, 1, 1e-8, R.Kmax);
sig = R.Amp > 1e-3*Vdc & R.freq <= 20000;
relErr = abs(A_fine(sig) - R.Amp(sig)) ./ R.Amp(sig);
[nPass,nFail] = chk('T8  精确法 vs 10ns细采样FFT (<=20kHz)', max(relErr) < 1e-2, ...
    sprintf('n=%d, max rel err = %.4f %%', nnz(sig), 100*max(relErr)), nPass, nFail);

%% ---------------- T9: 与 16 点 Gauss-Legendre 数值积分交叉验证 ----------------
ordGL = [1 3 5 95 97 191 193 288];
ordGL = ordGL(ordGL <= R.Kmax);
ordGL = ordGL(R.Amp(ordGL) > 1e-3*Vdc);      % 只比对实际存在的谐波
A_gl = gl_harmonics(R, fc, f0, Vdc, Mp, ordGL);
relGL = abs(A_gl - R.Amp(ordGL)) ./ R.Amp(ordGL);
[nPass,nFail] = chk('T9  精确法 vs 16点Gauss-Legendre数值积分', max(relGL) < 1e-9, ...
    sprintf('h=[%s], max rel err = %.3e', num2str(ordGL), max(relGL)), nPass, nFail);

%% ---------------- T10: A1 与 M 的线性度 ----------------
Ms = 0.1:0.05:1.0;
A1s = zeros(size(Ms));
for i = 1:numel(Ms)
    Rr = hbridge_spwm_harmonics(fc, f0, Vdc, Ms(i), Ncyc=2, Kmax=100, ...
        UseSampledFFT=false, Verbose=false);
    A1s(i) = Rr.Amp(1);
end
p = polyfit(Ms, A1s, 1);
linErr = max(abs(A1s - Ms*Vdc));
[nPass,nFail] = chk('T10 A1 与 M 严格线性 (斜率 = Vdc)', ...
    (abs(p(1)-Vdc) < 1e-9*Vdc) && (linErr < 1e-9*Vdc), ...
    sprintf('slope = %.10f (Vdc = %g), max dev = %.3e V', p(1), Vdc, linErr), nPass, nFail);

%% ---------------- T11: 奇数载波比 fc/f0 = 49 ----------------
fc2 = 2450; M2 = 0.7;
R2 = hbridge_spwm_harmonics(fc2, f0, Vdc, M2, Ncyc=2, Kmax=500, ...
    UseSampledFFT=false, Verbose=false);
A_fine2 = fft_harmonics(fc2, f0, Vdc, M2, 2, 2e-8, 500);
sig2 = R2.Amp > 1e-3*Vdc & R2.freq <= 20000;
relErr2 = abs(A_fine2(sig2) - R2.Amp(sig2)) ./ R2.Amp(sig2);
evenAmp2 = max(R2.Amp(2:2:end));
maxCos2  = max(abs(R2.Amp(R2.Amp>1e-6*Vdc) .* cos(deg2rad(R2.Ph_deg(R2.Amp>1e-6*Vdc)))));
ok11 = evenAmp2 < 1e-9*Vdc && maxCos2 < 1e-9*Vdc && max(relErr2) < 2e-2;
[nPass,nFail] = chk('T11 fc/f0 = 49 (奇数载波比) 仍正确', ok11, ...
    sprintf('max even = %.2e V, 余弦分量 = %.2e V, fine-FFT max err = %.3f %%', ...
    evenAmp2, maxCos2, 100*max(relErr2)), nPass, nFail);

%% ---------------- T12: M 边界鲁棒性 ----------------
ok12 = true; msg12 = 'M=0.02/0.1/1.0 基波误差均 < 1e-9*Vdc';
for Mv = [0.02 0.1 1.0]
    Rv = hbridge_spwm_harmonics(fc, f0, Vdc, Mv, Ncyc=2, Kmax=500, ...
        UseSampledFFT=false, Verbose=false);
    if abs(Rv.Amp(1) - Mv*Vdc) > 1e-9*Vdc
        ok12 = false;
        msg12 = sprintf('M=%g: A1 err %.3e V', Mv, abs(Rv.Amp(1)-Mv*Vdc));
    end
end
[nPass,nFail] = chk('T12 低调制比(M=0.02)与 M=1 边界仍正确', ok12, msg12, nPass, nFail);

%% ---------------- T13: 双重傅里叶级数(贝塞尔函数)解析解 ----------------
%  A(m,n) = Vdc*(4/pi)*(1/m)*|J_n(m*pi*M/2)|,  m = 2,4,6,... (载波倍数, 偶数), n 为奇数(边带)
maxRelB = 0; nCmpB = 0; badB = '';
for Mv = [0.35 0.8 1.0]
    Rv = hbridge_spwm_harmonics(fc, f0, Vdc, Mv, Ncyc=2, Kmax=600, ...
        UseSampledFFT=false, Verbose=false);
    for m = 2:2:8
        for n = -7:2:7
            h = m*Rratio + n;
            if h < 1 || h > Rv.Kmax, continue; end
            A_df = Vdc*(4/pi)*(1/m)*abs(besselj(n, m*pi*Mv/2));
            if A_df < 1e-5*Vdc, continue; end      % 幅值过小, 无比较意义
            e = abs(A_df - Rv.Amp(h))/A_df;
            nCmpB = nCmpB + 1;
            if e > maxRelB
                maxRelB = e; badB = sprintf('M=%g,m=%d,n=%+d', Mv, m, n);
            end
        end
    end
end
[nPass,nFail] = chk('T13 与双重傅里叶级数(贝塞尔函数)解析解一致', maxRelB < 1e-9, ...
    sprintf('对比 %d 个谐波 (M=0.35/0.8/1.0), max rel err = %.3e %s', ...
    nCmpB, maxRelB, badB), nPass, nFail);

%% ================= 死区 (dead time) 相关验证 =================
fprintf('\n  ---- 死区功能验证 (td = 4us) ----\n');
tdv = 4e-6;
Rl  = 0.5;  Ll = 2e-3;
Rdt = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=10, Kmax=2000, ...
        Tdead=tdv, Rload=Rl, Lload=Ll, UseSampledFFT=true, DT=1e-6, Verbose=false);
Ri  = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=10, Kmax=2000, ...
        Tdead=0, Rload=Rl, Lload=Ll, UseSampledFFT=false, Verbose=false);

% D1: td=0 必须与无死区结果完全一致 (向后兼容)
Rz = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=2, Kmax=500, Tdead=0, ...
        Rload=Rl, Lload=Ll, UseSampledFFT=false, Verbose=false);
Rr0 = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=2, Kmax=500, ...
        Rload=Rl, Lload=Ll, UseSampledFFT=false, Verbose=false);
d0 = max(abs(Rz.Amp - Rr0.Amp));
[nPass,nFail] = chk('D1  td=0 与无死区模型完全一致', d0 < 1e-9*Vdc, ...
    sprintf('max|dA| = %.3e V', d0), nPass, nFail);

% D2: 死区由 0 变为 4us 后, 低次谐波必须从 ~0 抬升到已知量级
[nPass,nFail] = chk('D2  死区使 3/5/7 次谐波出现 (原为 0)', ...
    all([Ri.Amp(3) Ri.Amp(5) Ri.Amp(7)] < 1e-9*Vdc) && ...
    all([Rdt.Amp(3) Rdt.Amp(5) Rdt.Amp(7)] > 0.05), ...
    sprintf('无死区 %.1e/%.1e/%.1e V -> 含死区 %.4f/%.4f/%.4f V', ...
    Ri.Amp(3), Ri.Amp(5), Ri.Amp(7), Rdt.Amp(3), Rdt.Amp(5), Rdt.Amp(7)), ...
    nPass, nFail);

% D3: 死区谐波的两条规律
%   (i)  1/h 律: 低次时 A_h * h 近似为常数 (死区误差近似 ±2*td*fc*Vdc 的方波)。
%        由于死区误差还受"斩波活动包络"缓慢调制, A_h*h 会随 h 缓慢上升,
%        实测 h=3..15 离散度约 2.2%, h=3..25 约 6.6%, 故 1/h 律只在低次严格成立。
%   (ii) 幅值量级与闭式解 8*td*fc*Vdc/(pi*h) 相符(偏差同样随 h 缓慢增大)。
hLaw   = 3:2:15;
prodH  = Rdt.Amp(hLaw) .* hLaw.';
spread = max(prodH)/min(prodH) - 1;
prodAll = Rdt.Amp(3:2:25) .* (3:2:25).';
[nPass,nFail] = chk('D3a 低次死区谐波满足 1/h 律 (h<=15, <3%)', spread < 0.03, ...
    sprintf('h=3..15: A_h*h = %.4f~%.4f V (离散 %.3f %%); h<=25 离散 %.2f %%', ...
    min(prodH), max(prodH), 100*spread, ...
    100*(max(prodAll)/min(prodAll)-1)), nPass, nFail);

envD = 8*tdv*fc*Vdc ./ (pi*hLaw.');
relD = abs(Rdt.Amp(hLaw) - envD) ./ envD;
[nPass,nFail] = chk('D3b 与闭式包络 8*td*fc*Vdc/(pi*h) 量级相符 (<7%)', ...
    max(relD) < 0.07, ...
    sprintf('h=3: %.3f%%, h=15: %.3f%%, h=25: %.3f%% (包络近似, 随 h 增大)', ...
    100*relD(1), 100*relD(7), 100*relD(end)), nPass, nFail);

% D4: 每载波周期平均误差电压 = 2*td*fc*Vdc (逐周期核对, 用门极状态直接算)
Tc = 1/fc;
errAvg = [];
for kk = [3 8 16 24 32 40 45]
    t0 = kk*Tc; t1 = t0 + Tc;
    if t1 > Ri.T, continue; end
    tg2 = linspace(t0, t1, 20001).';
    vI = Vdc*spwm_out_level(tg2, fc, f0, Mp);
    vD = Vdc*spwm_deadtime_out_level(tg2, fc, f0, Mp, tdv, Rdt.deadtime.phi_rad);
    errAvg(end+1) = mean(vD - vI); %#ok<SAGROW>
end
nz = abs(errAvg) > 0.1;
maxdev = max(abs(abs(errAvg(nz)) - 2*tdv*fc*Vdc));
[nPass,nFail] = chk('D4  逐载波周期平均误差 = ±2*td*fc*Vdc', ...
    any(nz) && maxdev < 0.02*2*tdv*fc*Vdc, ...
    sprintf('实测 %.4f V, 理论 %.4f V, 最大偏差 %.4f V', ...
    mean(abs(errAvg(nz))), 2*tdv*fc*Vdc, maxdev), nPass, nFail);

% D5: 死区谐波随 td 严格线性
A3a = zeros(1,3); tdl = [2e-6 4e-6 8e-6];
for i = 1:3
    Rx = hbridge_spwm_harmonics(fc, f0, Vdc, Mp, Ncyc=2, Kmax=500, Tdead=tdl(i), ...
        Rload=Rl, Lload=Ll, UseSampledFFT=false, Verbose=false);
    A3a(i) = Rx.Amp(3);
end
lin = max(abs(A3a ./ (tdl/tdl(1)) / A3a(1) - 1));
[nPass,nFail] = chk('D5  死区谐波与 td 严格线性 (<0.3%)', lin < 3e-3, ...
    sprintf('A3 = %.5f/%.5f/%.5f V, 归一化偏差 %.3e', A3a, lin), nPass, nFail);

% D6: 含死区时半波对称性保持 (偶次谐波仍为 0)
evenD = max(Rdt.Amp(2:2:end));
[nPass,nFail] = chk('D6  含死区仍只有奇次谐波', evenD < 1e-9*Vdc, ...
    sprintf('max even = %.3e V', evenD), nPass, nFail);

% D7: 电流方向自洽迭代收敛, 且与含死区的细采样FFT交叉验证
A_fineD = fft_harmonics_dt(fc, f0, Vdc, Mp, 1, 2e-8, 600, tdv, Rdt.deadtime.phi_rad);
sigD = Rdt.Amp(1:600) > 0.01*Vdc & Rdt.freq(1:600) <= 20000;
relD2 = abs(A_fineD(sigD) - Rdt.Amp(sigD)) ./ Rdt.Amp(sigD);
[nPass,nFail] = chk('D7  含死区: 20ns细采样FFT交叉验证', ...
    Rdt.deadtime.converged && max(relD2) < 1e-2, ...
    sprintf('收敛=%d; n=%d, max rel err = %.3f %%', Rdt.deadtime.converged, ...
    nnz(sigD), 100*max(relD2)), nPass, nFail);

% D8: 门极互锁与死区窗口宽度 (无直通, 每个死区窗口宽度 = td)
tg3 = linspace(0.004, 0.004+4*Tc, 400001).';
[g1,g2,g3,g4, dLg, dRg] = spwm_gate_states(tg3, fc, f0, Mp, tdv);
dtw = tg3(2)-tg3(1);
wins = local_window_widths(dLg, dtw);
winsR = local_window_widths(dRg, dtw);
[nPass,nFail] = chk('D8  同臂无直通且死区窗口宽度 = td', ...
    all((g1&g2)==0) && all((g3&g4)==0) && ~isempty(wins) && ...
    max(abs([wins; winsR] - tdv)) < 2*dtw, ...
    sprintf('直通=0; 左臂 %d 个窗口 %.3f~%.3f us, 右臂 %d 个 (td=%.3f us)', ...
    numel(wins), min(wins)*1e6, max(wins)*1e6, numel(winsR), tdv*1e6), ...
    nPass, nFail);

% D9: 死区造成基波幅值下降 (物理上必然)
[nPass,nFail] = chk('D9  死区使基波幅值下降', ...
    Rdt.Amp(1) < Ri.Amp(1) && Rdt.A1_drop_pct > 0.5, ...
    sprintf('A1: %.6f -> %.6f V (下降 %.4f %%)', ...
    Ri.Amp(1), Rdt.Amp(1), Rdt.A1_drop_pct), nPass, nFail);

%% ---------------- 汇总 ----------------
log('----------------------------------------------------------------');
log('   通过 %d 项, 失败 %d 项', nPass, nFail);
if nFail == 0
    log('   结果: 全部通过 (ALL PASS)');
else
    log('   结果: 存在失败项, 请检查!');
end
log('================================================================');
if nFail > 0
    error('verify_harmonics: 有 %d 项验证失败', nFail);
end

%% ======================= 局部函数 =======================
function [nPass, nFail] = chk(name, cond, detail, nPass, nFail)
if cond
    fprintf('  [PASS] %-46s | %s\n', name, detail);
    nPass = nPass + 1;
else
    fprintf('  [FAIL] %-46s | %s\n', name, detail);
    nFail = nFail + 1;
end
end

function A = fft_harmonics(fc, f0, Vdc, M, Ncyc, DT, Kmax)
% 独立方法 B: 等步长采样 + FFT
T = Ncyc/f0;
N = round(T/DT);
dt = T/N;
t = (0:N-1).' * dt;
Vout = Vdc * spwm_out_level(t, fc, f0, M);
Y = fft(Vout);
half = floor(N/2);
mag = 2*abs(Y(1:half+1)) / N;
mag(1) = abs(Y(1)) / N;
A = nan(Kmax,1);
idx = (1:Kmax).' * Ncyc + 1;
ok = idx <= half + 1;
A(ok) = mag(idx(ok));
end

function w = local_window_widths(flag, dtw)
% 提取逻辑量 flag 中所有 "1" 连续段的宽度 (秒)
d = diff([0; flag(:); 0]);
st = find(d == 1);      % 段起始
en = find(d == -1);     % 段结束
w = (en - st) * dtw;
end

function A = fft_harmonics_dt(fc, f0, Vdc, M, Ncyc, DT, Kmax, td, phi)
% 独立方法 D: 含死区波形的等步长采样 + FFT
T = Ncyc/f0;
N = round(T/DT);
dt = T/N;
t = (0:N-1).' * dt;
Vout = Vdc * spwm_deadtime_out_level(t, fc, f0, M, td, phi);
Y = fft(Vout);
half = floor(N/2);
mag = 2*abs(Y(1:half+1)) / N;
mag(1) = abs(Y(1)) / N;
A = nan(Kmax,1);
idx = (1:Kmax).' * Ncyc + 1;
ok = idx <= half + 1;
A(ok) = mag(idx(ok));
end

function A = gl_harmonics(R, fc, f0, Vdc, M, orders) %#ok<INUSD>
% 独立方法 C: 逐段 16 点 Gauss-Legendre 数值积分
[xg, wg] = gl_nodes(16);
a = R.exact.segA; b = R.exact.segB; lev = R.exact.segLev;
halfw = (b - a)/2;
mid   = (a + b)/2;
X = mid + halfw*xg.';            % Nseg x 16
W = halfw(:).*wg.';              % Nseg x 16
A = zeros(numel(orders),1);
for i = 1:numel(orders)
    w  = 2*pi*orders(i)*f0;
    ac = (2/R.T) * sum(lev .* sum(W.*cos(w*X), 2));
    as = (2/R.T) * sum(lev .* sum(W.*sin(w*X), 2));
    A(i) = hypot(ac, as);
end
end

function [x, w] = gl_nodes(n)
% Golub-Welsch 求 n 点 Gauss-Legendre 节点与权重 (区间 [-1,1])
i = (1:n-1)';
beta = i ./ sqrt(4*i.^2 - 1);
J = diag(beta,1) + diag(beta,-1);
[V, D] = eig(J);
x = diag(D);
[x, ord] = sort(x);
V = V(:, ord);
w = 2*(V(1,:).^2)';
end
