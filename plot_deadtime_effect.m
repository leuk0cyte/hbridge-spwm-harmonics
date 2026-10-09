function f = plot_deadtime_effect(Rdt, R0, tds, tsweep)
%PLOT_DEADTIME_EFFECT  图6: 死区对谐波的影响
%   (a) 理想 vs 含死区 低次谐波对比 (死区把 3,5,7... 次谐波从 0 抬起来)
%   (b) 含死区频谱 + 理论 8*td*fc*Vdc/(pi*h) 包络
%   (c) 死区谐波随 td 的线性增长
%   (d) 死区造成的基波损失与 THD 增长随 td 的变化
%
%   f = plot_deadtime_effect(Rdt, R0)                      tsweep 由函数内部扫描
%   f = plot_deadtime_effect(Rdt, R0, [0 1 2 4 8]*1e-6)
%   f = plot_deadtime_effect(Rdt, R0, tds, tsweep)

arguments
    Rdt struct
    R0  struct
    tds (1,:) double = [1 2 4 8]*1e-6
    tsweep = []
end
set_cjk_font();

td = Rdt.Tdead;
fc = Rdt.fc; f0 = Rdt.f0; Vdc = Rdt.Vdc;

f = figure('Color','w','Position',[60 30 1250 820]);
tl = tiledlayout(f, 2, 2, 'TileSpacing','compact');

% ---- (a) 低次谐波对比 ----
nexttile; hold on; grid on; box on
lo = 1:2:31;
b = bar([lo(:), lo(:)], [R0.Amp(lo), Rdt.Amp(lo)], 1.0, 'grouped');
b(1).FaceColor = [0.55 0.70 0.85]; b(1).EdgeColor = 'none';
b(2).FaceColor = [0.85 0.25 0.20]; b(2).EdgeColor = 'none';
xticks(lo); xtickangle(0);
xlabel('谐波次数 h'); ylabel('幅值 / V  (线性坐标)');
title(sprintf('(a) 低次谐波: 无死区(蓝) vs 含 %g \\mus 死区(红)', td*1e6));
legend({'无死区 (理论为 0)', sprintf('含死区 td=%g \\mus', td*1e6)}, ...
    'Location','northeast');
ymax = max(Rdt.Amp(lo));
ylim([0, ymax*1.35]);
idealMax = max(R0.Amp(lo(2:end)));      % 只取低次谐波(不含基波)
text(0.30, 0.60, sprintf(['无死区的 3,5,7... 次谐波幅值 <= %.1e V\n' ...
    '(在 0~%.1f V 的线性坐标下不可见, 即"不存在")'], idealMax, ymax), ...
    'Units','normalized','FontSize',9);
text(0.30, 0.42, sprintf('含死区 3 次谐波 = %.4f V (占基波 %.3f%%)', ...
    Rdt.Amp(3), 100*Rdt.Amp(3)/Rdt.Amp(1)), 'Units','normalized','FontSize',9);

% ---- (b) 频谱与理论包络 ----
nexttile; hold on; grid on; box on
sel = Rdt.freq <= min(3*fc, Rdt.Kmax*f0);
stem(Rdt.freq(sel)/1000, Rdt.Amp(sel), 'filled', 'MarkerSize', 3, ...
    'Color', [0.85 0.25 0.20], 'DisplayName', sprintf('含死区 td=%g \\mus', td*1e6));
h = (1:max(Rdt.harm)).';
env = 8*td*fc*Vdc ./ (pi*h);
selE = h >= 1 & h <= 51;
stairs(Rdt.freq(selE)/1000, env(selE), 'k--', 'LineWidth', 1.2, ...
    'DisplayName', '理论包络 8t_df_cV_{dc}/(\pih)');
set(gca,'YScale','log');
ylim([1e-4, max(Rdt.Amp)*1.5]);
xlabel('频率 / kHz'); ylabel('幅值 / V (对数)');
title('(b) 含死区频谱与死区谐波理论包络');
legend('Location','northeast');

% ---- (c) 死区谐波随 td 线性增长 ----
nexttile; hold on; grid on; box on
A3 = zeros(size(tds)); A1loss = zeros(size(tds)); thd50 = zeros(size(tds));
for i = 1:numel(tds)
    Rt = hbridge_spwm_harmonics(fc, f0, Vdc, Rdt.M, Ncyc=2, Kmax=500, ...
        Tdead=tds(i), Rload=Rdt.Rload, Lload=Rdt.Lload, ...
        UseSampledFFT=false, Verbose=false);
    A3(i) = Rt.Amp(3); A1loss(i) = Rt.A1_drop_pct; thd50(i) = 100*Rt.THD_50;
end
plot(tds*1e6, A3, 'o-', 'LineWidth', 1.4, 'Color', [0.85 0.25 0.20]);
plot(tds*1e6, (8*fc*Vdc/(3*pi))*tds, 'k--', 'LineWidth', 1.0);
xlabel('死区 t_d / \mus'); ylabel('3 次谐波幅值 / V');
title('(c) 3次谐波随死区线性增长 (虚线=理论斜率 8f_cV_{dc}/(3\pi))');
legend({'精确计算', '理论'}, 'Location','northwest');

% ---- (d) 基波损失与THD ----
nexttile; hold on; grid on; box on
yyaxis left
plot(tds*1e6, A1loss, 's-', 'LineWidth', 1.4);
ylabel('基波幅值损失 / %');
yyaxis right
plot(tds*1e6, thd50, '^-', 'LineWidth', 1.4);
ylabel('THD (\leq50次) / %');
xlabel('死区 t_d / \mus');
title(sprintf('(d) 基波损失(左轴) 与 低次谐波THD(右轴) 随死区变化  (M=%g)', Rdt.M));

title(tl, sprintf(['死区对H桥单极性SPWM输出电压谐波的影响  ' ...
    '(fc=%g Hz, f0=%g Hz, Vdc=%g V, M=%g, td=%g \\mus)'], ...
    fc, f0, Vdc, Rdt.M, td*1e6), 'FontWeight','bold');
end