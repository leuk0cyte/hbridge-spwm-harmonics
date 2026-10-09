function f = plot_spectrum_detail(R)
%PLOT_SPECTRUM_DETAIL  图3: 频谱细节
%   (a) 0~2kHz  低次谐波(理论上仅剩基波)
%   (b) 2fc 附近谐波群
%   (c) 4fc 附近谐波群
%   (d) 幅值最大的20个谐波(柱状)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
set_cjk_font();

f = figure('Color','w','Position',[60 30 1250 780]);
tl = tiledlayout(f, 2, 2, 'TileSpacing','compact');

% ---- (a) 低频段 ----
nexttile; hold on; grid on; box on
sel = R.freq <= 2000;
stem(R.freq(sel)/1000, R.Amp(sel), 'filled', 'MarkerSize', 4, 'Color', [0.85 0.15 0.15]);
set(gca,'YScale','log'); ylim([1e-9 max(R.Amp)*1.5]);
xlim([0 2]);
xlabel('频率 / kHz'); ylabel('幅值 / V (对数)');
title('(a) 0~2 kHz: 低次谐波 (单极性SPWM 中 3,5,7... 理论上为 0)');
lowAmp = max(R.Amp(3:2:min(numel(R.Amp), floor(2000/R.f0))), [], 'omitnan');
text(0.35, 0.72, sprintf('低次奇次谐波最大幅值 = %.3e V (%.2e %%基波)', ...
    lowAmp, 100*lowAmp/R.Amp(1)), 'Units','normalized', 'FontSize', 9);

% ---- (b) 2fc 附近 ----
nexttile; hold on; grid on; box on
hb = local_win(R, 2*R.fc, 8*R.f0);
local_cluster_plot(R, hb);
title(sprintf('(b) 2fc = %g Hz 附近谐波群 (fc/f0=%g, 边带间隔 2f0 = %g Hz)', ...
    2*R.fc, R.carrierRatio, 2*R.f0));

% ---- (c) 4fc 附近 ----
nexttile; hold on; grid on; box on
hc = local_win(R, 4*R.fc, 8*R.f0);
local_cluster_plot(R, hc);
title(sprintf('(c) 4fc = %g Hz 附近谐波群', 4*R.fc));

% ---- (d) 最大谐波柱状图 ----
nexttile; hold on; grid on; box on
[As, is] = sort(R.Amp, 'descend');
n = min(20, numel(As));
b = bar(1:n, As(1:n), 0.7, 'FaceColor', [0.25 0.55 0.85]);
b.EdgeColor = 'none';
xticks(1:n);
xticklabels(arrayfun(@(h) sprintf('h=%d', h), is(1:n), 'UniformOutput', false));
xtickangle(60); set(gca,'FontSize',7);
ylabel('幅值 / V'); xlabel('谐波次数 h (按幅值降序,  频率 = h \times 50 Hz)');
title('(d) 幅值最大的 20 个谐波');

title(tl, sprintf('输出电压频谱细节  (M=%g, fc=%g Hz, f0=%g Hz)', R.M, R.fc, R.f0), ...
    'FontWeight','bold');
end

function hb = local_win(R, fcenter, halfwidth)
f1 = fcenter - halfwidth; f2 = fcenter + halfwidth;
hb = find(R.freq >= f1 & R.freq <= f2);
end

function local_cluster_plot(R, hb)
if isempty(hb)
    return
end
[freqs, order] = sort(R.freq(hb));
hb = hb(order);
stem(freqs/1000, R.Amp(hb), 'filled', 'MarkerSize', 4, 'Color', [0.85 0.15 0.15]);
for k = 1:numel(hb)
    if R.Amp(hb(k)) > 1e-4*R.Vdc
        text(freqs(k)/1000, R.Amp(hb(k)), sprintf('h=%d', R.harm(hb(k))), ...
            'FontSize', 7, 'VerticalAlignment','bottom', 'HorizontalAlignment','center');
    end
end
xlabel('频率 / kHz'); ylabel('幅值 / V');
xlim([min(freqs) max(freqs)]/1000);
ylim([0 max(R.Amp(hb))*1.25 + eps]);
end
