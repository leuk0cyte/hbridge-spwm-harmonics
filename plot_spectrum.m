function f = plot_spectrum(R, Fmax)
%PLOT_SPECTRUM  图2: 电压谐波频谱 (线性谱 / 对数谱 / 与原1us采样FFT的误差)
%
%   f = plot_spectrum(R)               默认显示到 2.5*fc
%   f = plot_spectrum(R, 12000)        指定最高显示频率 Hz
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    R struct
    Fmax (1,1) double = 0
end
if Fmax <= 0
    Fmax = 2.5 * R.fc;
end
set_cjk_font();

sel  = R.freq <= Fmax;
fHz  = R.freq(sel) / 1000;      % kHz
Aex  = R.Amp(sel);
hasFFT = ~isempty(R.Amp_fft);
if hasFFT
    Afft = R.Amp_fft(sel);
    Err  = R.err_rel_fft(sel) * 100;   % %
end

f = figure('Color','w','Position',[80 40 1150 880]);
tl = tiledlayout(f, 3, 1, 'TileSpacing','compact');

% ---- (a) 线性幅值谱 ----
nexttile; hold on; grid on; box on
stem(fHz, Aex, 'filled', 'MarkerSize', 3, 'Color', [0.85 0.15 0.15], 'DisplayName','精确解析法');
if hasFFT
    plot(fHz, Afft, 'o', 'MarkerSize', 4, 'LineWidth', 1, ...
        'Color', [0.10 0.35 0.90], 'DisplayName','1\mus 采样FFT法 (原Python仿真步长)');
end
for m = 2:2:8
    xline(m*R.fc/1000, ':', sprintf('%dfc', m), 'Color', [0.4 0.4 0.4], ...
        'LabelOrientation','horizontal', 'FontSize', 8, 'HandleVisibility','off');
end
xlim([0 Fmax/1000]); ylim([0 max(Aex)*1.12 + eps]);
xlabel('频率 / kHz'); ylabel('谐波幅值 / V');
title(sprintf('(a) 输出电压谐波幅值谱  (fc=%g Hz, f0=%g Hz, M=%g, Vdc=%g V)', ...
    R.fc, R.f0, R.M, R.Vdc));
legend('Location','northeast');

% ---- (b) 对数幅值谱 ----
nexttile; hold on; grid on; box on
stem(fHz(Aex>0), Aex(Aex>0), 'filled', 'MarkerSize', 3, 'Color', [0.85 0.15 0.15], 'DisplayName','精确解析法');
if hasFFT
    plot(fHz(Afft>0), Afft(Afft>0), 'o', 'MarkerSize', 4, 'Color', [0.10 0.35 0.90], ...
        'DisplayName','1\mus 采样FFT法');
end
for m = 2:2:8
    xline(m*R.fc/1000, ':', sprintf('%dfc', m), 'Color', [0.4 0.4 0.4], ...
        'LabelVerticalAlignment','bottom', 'FontSize', 8, 'HandleVisibility','off');
end
set(gca, 'YScale','log');
xlim([0 Fmax/1000]); ylim([1e-6 10^(ceil(log10(max(Aex)*1.2)))]);
xlabel('频率 / kHz'); ylabel('谐波幅值 / V (对数)');
title('(b) 对数幅值谱 (可见低次谐波几乎为零, 能量集中在 2fc, 4fc 附近)');

% ---- (c) 采样FFT误差 (只画显著谐波, 微小谐波的相对误差无意义) ----
nexttile; hold on; grid on; box on
if hasFFT
    sigF = R.Amp(sel) > 0.01*R.Vdc;      % |A_h| > 1% Vdc
    semilogy(fHz(sigF), Err(sigF), 'o-', 'MarkerSize', 4, 'LineWidth', 1, ...
        'Color', [0.55 0.15 0.65]);
    yline(1, '--', '1%', 'Color', [0.5 0.5 0.5], 'HandleVisibility','off');
    ylabel('相对误差 / %');
    title(sprintf(['(c) 1\\mus 等步长采样FFT 相对精确法的误差 ' ...
        '(仅显著谐波 |A_h|>1%%V_{dc}, 最大 %.2f%%)'], R.Amp_fft_maxerr*100));
else
    text(0.5, 0.5, '未计算采样FFT (UseSampledFFT=false)', 'Units','normalized', ...
        'HorizontalAlignment','center');
end
set(gca, 'YScale','log'); xlim([0 Fmax/1000]);
xlabel('频率 / kHz');
set(gca, 'YLimMode','auto');

title(tl, sprintf('H桥单极性SPWM 输出电压谐波频谱  (M=%g, fc=%g Hz)', R.M, R.fc), ...
    'FontWeight','bold');
end
