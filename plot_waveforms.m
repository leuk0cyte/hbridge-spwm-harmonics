function f = plot_waveforms(R, tWin)
%PLOT_WAVEFORMS  图1: 三角载波/调制波, 4个IGBT驱动信号, 输出电压(含开关细节)
%
%   f = plot_waveforms(R)             默认显示2个基波周期
%   f = plot_waveforms(R, [0 0.04])   指定时间窗口(秒)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    R struct
    tWin (1,2) double = [0 0]
end
T0 = 1/R.f0;
if tWin(2) <= tWin(1)
    tWin = [0, 2*T0];
end
set_cjk_font();

tA = linspace(tWin(1), tWin(2), 40000).';
cA = spwm_tri_carrier(tA, R.fc);
vA = R.M * sin(2*pi*R.f0*tA);

tZoom = linspace(tWin(1), tWin(1) + min(2e-3, tWin(2)-tWin(1)), 40000).';
cZ = spwm_tri_carrier(tZoom, R.fc);
vZ = R.M * sin(2*pi*R.f0*tZoom);
Q1 = double(vZ >  cZ);
Q2 = double(vZ <= cZ);
Q3 = double(-vZ >  cZ);
Q4 = double(-vZ <= cZ);

tV = linspace(tWin(1), tWin(1) + min(2/R.fc, tWin(2)-tWin(1)), 40000).';
VV = R.Vdc * spwm_out_level(tV, R.fc, R.f0, R.M);

f = figure('Color','w','Position',[60 40 1150 900]);
tl = tiledlayout(f, 4, 1, 'TileSpacing','compact');

nexttile; hold on; grid on; box on
plot(tA*1000, cA, 'k-', 'LineWidth', 0.5);
plot(tA*1000, vA, 'r-', 'LineWidth', 1.2);
plot(tA*1000, -vA, 'b--', 'LineWidth', 1.2);
yline(0, 'Color', [0.6 0.6 0.6]);
xlim(tWin*1000); ylim([-1.35 1.35]);
ylabel('幅值'); title('(a) 三角载波 c(t) 与调制波 v_m(t), -v_m(t)');
legend('三角载波', '调制波 v_m', '反相调制波 -v_m', 'Location','northeast');

nexttile; hold on; grid on; box on
OFF = 1.35;
plot(tZoom*1000, Q1,         'r-', 'LineWidth', 1.0);
plot(tZoom*1000, Q2 + OFF,   'b-', 'LineWidth', 1.0);
plot(tZoom*1000, Q3 + 2*OFF, 'g-', 'LineWidth', 1.0);
plot(tZoom*1000, Q4 + 3*OFF, 'm-', 'LineWidth', 1.0);
xlim(tZoom([1 end])*1000);
ylim([-0.3 3*OFF+0.5]);
yticks([0.5 0.5+OFF 0.5+2*OFF 0.5+3*OFF]);
yticklabels({'Q1','Q2','Q3','Q4'});
ylabel('驱动信号');
title(sprintf('(b) H桥4个IGBT的GE驱动信号 (放大 %g ms 显示开关细节)', (tZoom(end)-tZoom(1))*1000));
legend({'Q1 左上','Q2 左下','Q3 右上','Q4 右下'}, 'Location','northeast','NumColumns',2);

nexttile; hold on; grid on; box on
[sa, sb, sL] = local_segments(R, tWin(1), tWin(2));
stairs([sa; sb(end)]*1000, [sL; sL(end)], 'g-', 'LineWidth', 1.0);
xlim(tWin*1000); ylim([-R.Vdc*1.15 R.Vdc*1.15]);
yline(R.Vdc, ':', 'Color', [0.6 0.6 0.6]); yline(-R.Vdc, ':', 'Color', [0.6 0.6 0.6]);
ylabel('电压 / V'); title('(c) H桥输出电压 V_{out} = Vdc(Q1-Q3), 三电平 (精确开关时刻)');

nexttile; hold on; grid on; box on
stairs(tV*1000, VV, 'g-', 'LineWidth', 1.3);
xlim(tV([1 end])*1000); ylim([-R.Vdc*1.2 R.Vdc*1.2]);
xlabel('时间 / ms'); ylabel('电压 / V');
title(sprintf('(d) 输出电压局部放大 (约 %g 个开关周期, 可见 ±Vdc 与 0V 续流态)', ...
    (tV(end)-tV(1))*R.fc));

title(tl, sprintf('H桥单极性SPWM波形  (fc=%g Hz, f0=%g Hz, M=%g, Vdc=%g V)', ...
    R.fc, R.f0, R.M, R.Vdc), 'FontWeight','bold');
end

function [a, b, L] = local_segments(R, t1, t2)
a = R.exact.segA; b = R.exact.segB; L = R.exact.segLev;
k = (b > t1) & (a < t2);
if ~any(k)
    a = t1; b = t2; L = 0; return
end
a = max(a(k), t1); b = min(b(k), t2); L = L(k);
end
