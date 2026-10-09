function S = harmonic_report(R, topN)
%HARMONIC_REPORT  打印H桥单极性SPWM谐波分析报告, 并返回关键指标结构体
%
%   S = harmonic_report(R)        打印全部内容
%   S = harmonic_report(R, 20)    指定"最大谐波"列表长度
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

if nargin < 2 || isempty(topN)
    topN = 15;
end

H   = R.harm;  A = R.Amp;  F = R.freq;  PH = R.Ph_deg;
pct = 100 * A / max(A(1), eps);

fprintf('\n');
fprintf('=================== 谐波分析报告 (精确解析法) ===================\n');
fprintf('  载波频率 fc = %g Hz   基波 f0 = %g Hz   调制比 M = %g   Vdc = %g V\n', ...
    R.fc, R.f0, R.M, R.Vdc);
fprintf('  载波比 fc/f0 = %g %s\n', R.carrierRatio, ...
    ternary_str(R.isIntegerRatio, '(整数: 仅奇次谐波, 纯正弦级数)', '(非整数, 存在频谱泄漏!)'));
fprintf('  基波幅值 A1 = %.6f V  (理论 M*Vdc = %.6f V)\n', A(1), R.A1_theory);
fprintf('  输出有效值 = %.4f V   直流分量 = %.3e V\n', R.RMS_time, R.DC);
fprintf('  -----------------------------------------------------------\n');
fprintf('  THD (<=50次)      = %8.3f %%\n', 100*R.THD_50);
fprintf('  THD (<=%d次, 2fc)  = %8.3f %%\n', min(floor(2*R.fc/R.f0),R.Kmax), 100*R.THD_2fc);
fprintf('  THD (<=%d次, 4fc)  = %8.3f %%\n', min(floor(4*R.fc/R.f0),R.Kmax), 100*R.THD_4fc);
fprintf('  THD (<=%d次, %gkHz)= %8.3f %%\n', R.Kmax, R.Kmax*R.f0/1000, 100*R.THD_all);
fprintf('  WTHD (1/h 加权)   = %8.3f %%\n', 100*R.WTHD);
if R.Rload > 0
    fprintf('  电流THD (R=%g ohm, L=%g H, <=50次) = %8.3f %%\n', ...
        R.Rload, R.Lload, 100*R.THD_i_50);
    fprintf('  电流THD (<=%d次, %gkHz)             = %8.3f %%\n', ...
        R.Kmax, R.Kmax*R.f0/1000, 100*R.THD_i_all);
end
fprintf('================================================================\n');

% ---------------- 低次谐波表 ----------------
fprintf('\n--- 低次谐波 (h = 1 ~ %d) ---\n', min(41, R.Kmax));
fprintf('    h      f/Hz       幅值/V     占基波%%     相位/deg\n');
fprintf('  -----------------------------------------------------\n');
for k = 1:min(41, R.Kmax)
    if A(k) > 1e-6*R.Vdc
        fprintf('  %3d  %9.1f  %11.6f  %9.4f  %10.2f\n', k, F(k), A(k), pct(k), PH(k));
    else
        fprintf('  %3d  %9.1f  %11.6f  %9.4f  %10s\n', k, F(k), A(k), pct(k), '--');
    end
end
fprintf('  (幅值小于 1e-6*Vdc 的谐波视为 0, 其相位无意义, 以 -- 表示)\n');

% ---------------- 幅值最大的谐波 ----------------
[As, is] = sort(A, 'descend');
ns = min(topN, numel(As));
fprintf('\n--- 幅值最大的 %d 个谐波 ---\n', ns);
fprintf('    h      f/Hz       幅值/V     占基波%%\n');
fprintf('  ---------------------------------------\n');
for k = 1:ns
    h = is(k);
    fprintf('  %3d  %9.1f  %11.6f  %9.4f\n', h, F(h), A(h), pct(h));
end

% ---------------- 载波谐波群 ----------------
Rr = round(R.carrierRatio);
S.clusters = struct('m', {}, 'center_h', {}, 'harm', {}, 'amp', {});
if R.isIntegerRatio
    fprintf('\n--- 载波谐波群: h = %d*m ± 奇数 (单极性SPWM的谐波只出现在偶数倍载波附近) ---\n', Rr);
    fprintf('  %6s %10s %10s %12s\n', 'h', 'f/Hz', '幅值/V', '占基波%');
    fprintf('  ---------------------------------------------\n');
    for m = 2:2:8
        hc = m * Rr;
        if hc - 8 > R.Kmax, break; end
        fprintf('  ---- m = %d (中心 %g Hz) ----\n', m, hc*R.f0);
        offs = -7:2:7;
        hh = hc + offs;
        hh = hh(hh >= 1 & hh <= R.Kmax);
        for h = hh
            fprintf('  %6d %10.1f %10.5f %12.4f\n', h, F(h), A(h), pct(h));
        end
        S.clusters(end+1) = struct('m', m, 'center_h', hc, 'harm', hh, 'amp', A(hh)); %#ok<AGROW>
    end
    fprintf('  ---- 中心 %d Hz (h=%d, m=2 中心点) ----\n', 2*R.fc, 2*Rr);
    if 2*Rr <= R.Kmax
        fprintf('  %6d %10.1f %10.5f %12.4f\n', 2*Rr, F(2*Rr), A(2*Rr), pct(2*Rr));
    end
end

% ---------------- 返回值 ----------------
S.A1 = A(1);
S.A1_theory = R.A1_theory;
S.DC = R.DC;
S.RMS_time = R.RMS_time;
S.THD_50 = R.THD_50; S.THD_2fc = R.THD_2fc; S.THD_4fc = R.THD_4fc;
S.THD_all = R.THD_all; S.WTHD = R.WTHD;
S.THD_i_50 = R.THD_i_50; S.THD_i_all = R.THD_i_all;
S.top_harm = is(1:ns);
S.top_amp  = As(1:ns);
S.low_order_pct = pct(1:min(41, R.Kmax));
if ~isempty(R.Amp_fft)
    S.fft_max_rel_err = R.Amp_fft_maxerr;
end
end

function s = ternary_str(c, a, b)
if c, s = a; else, s = b; end
end
