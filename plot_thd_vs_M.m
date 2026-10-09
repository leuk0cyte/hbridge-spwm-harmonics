function f = plot_thd_vs_M(S)
%PLOT_THD_VS_M  图4: 调制比 M 对基波幅值/THD/WTHD/电流THD的影响
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
set_cjk_font();

f = figure('Color','w','Position',[80 60 1150 760]);
tl = tiledlayout(f, 2, 2, 'TileSpacing','compact');

% ---- (a) 基波幅值 ----
nexttile; hold on; grid on; box on
plot(S.M, S.A1_theory, 'k-', 'LineWidth', 1.5, 'DisplayName','理论 A_1 = M\cdotV_{dc}');
plot(S.M, S.A1, 'ro', 'MarkerSize', 6, 'LineWidth', 1, 'DisplayName','精确法计算值');
xlabel('调制比 M'); ylabel('基波幅值 / V');
title('(a) 基波幅值与调制比 (线性调制区应为直线)');
legend('Location','northwest');
err = max(abs(S.A1 - S.A1_theory));
text(0.35, 0.25, sprintf('最大偏差 = %.3e V', err), 'Units','normalized', 'FontSize', 9);

% ---- (b) 电压THD ----
nexttile; hold on; grid on; box on
plot(S.M, 100*S.THD_2fc,  's-', 'LineWidth', 1.2, 'DisplayName','THD (\leq2f_c)');
plot(S.M, 100*S.THD_4fc,  'd-', 'LineWidth', 1.2, 'DisplayName','THD (\leq4f_c)');
plot(S.M, 100*S.THD_all,  '^-', 'LineWidth', 1.2, ...
    'DisplayName', sprintf('THD (\\leq%g kHz)', S.Kmax*S.f0/1000));
xlabel('调制比 M'); ylabel('电压 THD / %');
set(gca,'YScale','log'); set(gca,'YLimMode','auto');
title('(b) 电压 THD 与调制比 (M 越小, 基波越小, THD 越大)');
legend('Location','northeast');
text(0.03, 0.06, sprintf('THD (\\leq50次) = 0 %%  (低次谐波理论为 0, 实测 < %.0e %%)', ...
    max(100*S.THD_50)), 'Units','normalized', 'FontSize', 9);

% ---- (c) WTHD ----
nexttile; hold on; grid on; box on
plot(S.M, 100*S.WTHD, 'd-', 'LineWidth', 1.5, 'Color', [0.85 0.35 0.1]);
xlabel('调制比 M'); ylabel('WTHD / %');
set(gca,'YScale','log');
title('(c) 1/h 加权 THD (更能反映电感负载电流畸变)');

% ---- (d) 电流THD ----
nexttile; hold on; grid on; box on
if any(isfinite(S.THD_i_50))
    plot(S.M, 100*S.THD_i_all, 's-', 'LineWidth', 1.2, ...
        'DisplayName', sprintf('电流THD (\\leq%g kHz)', S.Kmax*S.f0/1000));
    set(gca,'YScale','log'); set(gca,'YLimMode','auto');
    legend('Location','northeast');
    text(0.03, 0.06, sprintf('电流THD (\\leq50次) = 0 %%  (实测 < %.0e %%)', ...
        max(100*S.THD_i_50)), 'Units','normalized', 'FontSize', 9);
else
    text(0.5,0.5,'未设置负载 (Rload=0)', 'Units','normalized','HorizontalAlignment','center');
end
xlabel('调制比 M'); ylabel('电流 THD / %');
title('(d) 负载电流 THD (R-L 负载)');

title(tl, sprintf('调制比 M 扫描结果  (fc=%g Hz, f0=%g Hz, Vdc=%g V)', S.fc, S.f0, S.Vdc), ...
    'FontWeight','bold');
end
