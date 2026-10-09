function S = harmonic_report(R, topN, Rref)
%HARMONIC_REPORT  打印H桥单极性SPWM谐波分析报告, 并返回关键指标结构体
%
%   S = harmonic_report(R)              打印全部内容
%   S = harmonic_report(R, 20)          指定"最大谐波"列表长度
%   S = harmonic_report(R, 20, Rref)    与无死区的参考工况 Rref 对比(死区影响)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

if nargin < 2 || isempty(topN)
    topN = 15;
end
doCmp = (nargin >= 3) && ~isempty(Rref);

H   = R.harm;  A = R.Amp;  F = R.freq;  PH = R.Ph_deg;
pct = 100 * A / max(A(1), eps);

fprintf('\n');
fprintf('=================== 谐波分析报告 (精确解析法) ===================\n');
fprintf('  载波频率 fc = %g Hz   基波 f0 = %g Hz   调制比 M = %g   Vdc = %g V\n', ...
    R.fc, R.f0, R.M, R.Vdc);
fprintf('  载波比 fc/f0 = %g %s\n', R.carrierRatio, ...
    ternary_str(R.isIntegerRatio, '(整数: 仅奇次谐波, 纯正弦级数)', '(非整数, 存在频谱泄漏!)'));
fprintf('  基波幅值 A1 = %.6f V  (理论 M*Vdc = %.6f V)\n', A(1), R.A1_theory);
if isfield(R,'Tdead') && R.Tdead > 0
    fprintf('  开关死区 td = %g us (占载波周期 %.2f %%)\n', R.Tdead*1e6, 100*R.td_over_Tc);
    fprintf('  每载波周期平均误差电压 = %.4f V (= 2*td*fc*Vdc)\n', R.deadtime.dV_per_Tc);
    fprintf('  死区造成基波损失 = %.4f V (%.4f %%)\n', R.A1_drop, R.A1_drop_pct);
end
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

% ---------------- 死区影响: 新出现的低次谐波 ----------------
if isfield(R,'Tdead') && R.Tdead > 0
    fprintf('\n--- 死区引入的低次谐波 (无死区时严格为 0) ---\n');
    fprintf('    h      f/Hz      含死区/V    占基波%%   8*td*fc*Vdc/(pi*h)   理论误差%%\n');
    fprintf('  ---------------------------------------------------------------------\n');
    lo = 3:2:min(25, R.Kmax);
    est = 8 * R.Tdead * R.fc * R.Vdc ./ (pi * lo);
    for k = 1:numel(lo)
        h = lo(k);
        fprintf('  %3d  %9.1f  %11.6f  %9.4f  %18.6f  %11.3f\n', ...
            h, F(h), A(h), pct(h), est(k), 100*(A(h)/est(k) - 1));
    end
    fprintf(['  说明: 死区误差电压每载波周期为 -2*sign(i)*td*fc*Vdc, 其傅里叶\n' ...
             '        奇次谐波幅值 ~ 8*td*fc*Vdc/(pi*h), 故低次谐波按 1/h 衰减。\n']);
end

% ---------------- 与无死区工况对比 ----------------
if doCmp
    fprintf('\n--- 死区前后对比 (参考工况 td = 0) ---\n');
    fprintf('  %-28s %14s %14s %12s\n', '指标', '无死区', '含死区', '变化');
    fprintf('  %-28s %14.4f %14.4f %11.3f%%\n', '基波幅值 / V', ...
        Rref.Amp(1), A(1), 100*(A(1)/Rref.Amp(1)-1));
    fprintf('  %-28s %14.4f %14.4f %12s\n', 'THD (<=50次) / %', ...
        100*Rref.THD_50, 100*R.THD_50, ...
        local_fmt_change(Rref.THD_50, R.THD_50));
    fprintf('  %-28s %14.4f %14.4f %12.3f%%\n', 'THD (<=2fc) / %', ...
        100*Rref.THD_2fc, 100*R.THD_2fc, 100*(R.THD_2fc/Rref.THD_2fc-1));
    fprintf('  %-28s %14.4f %14.4f %12.3f%%\n', 'THD (全部) / %', ...
        100*Rref.THD_all, 100*R.THD_all, 100*(R.THD_all/Rref.THD_all-1));
    fprintf('  %-28s %14.4f %14.4f %12.3f%%\n', 'WTHD / %', ...
        100*Rref.WTHD, 100*R.WTHD, 100*(R.WTHD/Rref.WTHD-1));
    if R.Rload > 0 && Rref.Rload > 0
        fprintf('  %-28s %14.4f %14.4f %12s\n', '电流THD (<=50次) / %', ...
            100*Rref.THD_i_50, 100*R.THD_i_50, ...
            local_fmt_change(Rref.THD_i_50, R.THD_i_50));
        fprintf('  %-28s %14.4f %14.4f %12.3f%%\n', '电流THD (全部) / %', ...
            100*Rref.THD_i_all, 100*R.THD_i_all, ...
            100*(R.THD_i_all/Rref.THD_i_all-1));
    end
    fprintf('  %-28s %14.2e %14.2e\n', '3次谐波幅值 / V', Rref.Amp(3), A(3));
    fprintf('  %-28s %14.2e %14.2e\n', '5次谐波幅值 / V', Rref.Amp(5), A(5));
    fprintf('  %-28s %14.2e %14.2e\n', '7次谐波幅值 / V', Rref.Amp(7), A(7));
    S.cmp.A1_ratio = A(1)/Rref.Amp(1);
    S.cmp.THD50_ratio = R.THD_50/max(Rref.THD_50,eps);
    S.cmp.A3 = A(3); S.cmp.A5 = A(5); S.cmp.A7 = A(7);
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

function v = local_ratio_change(base, new)
% 相对变化(%%); 若基准本身 ~0 (如无死区时的低次谐波THD), 返回 NaN 以便显示为 --
if abs(base) < 1e-9
    v = NaN;
else
    v = 100*(new/base - 1);
end
end

function s = local_fmt_change(base, new)
% 变化列文本: 基准~0 时显示 "从 ~0 抬升"
v = local_ratio_change(base, new);
if isnan(v)
    s = '从 ~0 抬升';
else
    s = sprintf('%.3f%%', v);
end
end
