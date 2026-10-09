function f = plot_spectra_vs_M(fc, f0, Vdc, Mlist, Kmax)
%PLOT_SPECTRA_VS_M  图5: 不同调制比 M 下的谐波频谱对比
%
%   f = plot_spectra_vs_M(2400, 50, 100, [0.2 0.4 0.6 0.8 1.0], 600)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
arguments
    fc    (1,1) double = 2400
    f0    (1,1) double = 50
    Vdc   (1,1) double = 100
    Mlist (1,:) double = [0.2 0.4 0.6 0.8 1.0]
    Kmax  (1,1) double = 600
end
set_cjk_font();

n = numel(Mlist);
A = zeros(Kmax, n);
for i = 1:n
    R = hbridge_spwm_harmonics(fc, f0, Vdc, Mlist(i), Ncyc=2, Kmax=Kmax, ...
        UseSampledFFT=false, Verbose=false);
    A(:,i) = R.Amp;
end
freq = (1:Kmax).' * f0 / 1000;      % kHz

f = figure('Color','w','Position',[80 60 1200 800]);
tl = tiledlayout(f, 3, 1, 'TileSpacing','compact');

nexttile; hold on; grid on; box on
for i = 1:n
    stem(freq, A(:,i), 'filled', 'MarkerSize', 2, 'DisplayName', sprintf('M=%g', Mlist(i)));
end
for m = 2:2:6
    xline(m*fc/1000, ':', sprintf('%dfc', m), 'Color', [0.4 0.4 0.4], 'FontSize', 8, ...
        'HandleVisibility','off');
end
xlim([0 min(Kmax*f0, 8*fc)/1000]); ylim([0 max(A(:))*1.1]);
xlabel('频率 / kHz'); ylabel('谐波幅值 / V');
title('(a) 不同调制比下的谐波幅值谱 (绝对幅值)');
legend('Location','northeast');

nexttile; hold on; grid on; box on
for i = 1:n
    stem(freq, 100*A(:,i)/A(1,i), 'filled', 'MarkerSize', 2, 'DisplayName', sprintf('M=%g', Mlist(i)));
end
xlim([0 min(Kmax*f0, 8*fc)/1000]);
xlabel('频率 / kHz'); ylabel('相对基波 / %');
title('(b) 以各自基波归一化 (可见 2fc 谐波群随 M 的变化)');
legend('Location','northeast');

nexttile; hold on; grid on; box on
hc = 2*round(fc/f0);
win = max(1,hc-15):min(Kmax,hc+15);
for i = 1:n
    plot((win*f0)/1000, A(win,i), 'o-', 'LineWidth', 1.2, ...
        'MarkerSize', 5, 'DisplayName', sprintf('M=%g', Mlist(i)));
end
xlabel('频率 / kHz'); ylabel('谐波幅值 / V');
title('(c) 2f_c 谐波群细节 (边带间隔 2f_0)');
legend('Location','northeast');

title(tl, sprintf('调制比 M 对输出电压频谱的影响  (fc=%g Hz, f0=%g Hz, Vdc=%g V)', ...
    fc, f0, Vdc), 'FontWeight','bold');
end
