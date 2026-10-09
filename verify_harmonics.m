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
