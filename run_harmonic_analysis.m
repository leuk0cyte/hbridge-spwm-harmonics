%% ========================================================================
%  H桥单极性SPWM (变频器输出) 谐波分析  ——  主脚本
%  对应原 Python 脚本: h_bridge_svpwm.py
%  用法: 修改下面 "参数设置" 后按 F5 运行, 或命令行执行 run_harmonic_analysis
%  输出: matlab_harmonics\output\ 目录下的 CSV 数据 + PNG 图形
%% ========================================================================
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
clear; clc; close all;

here = fileparts(mfilename('fullpath'));
if isempty(here), here = pwd; end
addpath(here);
outdir = fullfile(here, 'output');
if ~exist(outdir, 'dir'), mkdir(outdir); end

%% ---------------------------- 1. 参数设置 ----------------------------
fc    = 2400;        % 开关(载波)频率 [Hz]        <<< 可调 (原Python为400)
f0    = 50;          % 基波(调制波)频率 [Hz]
Vdc   = 100;         % 直流母线电压 [V]
M     = 0.8;         % 调制比 0.1 ~ 1.0           <<< 可调

Ncyc  = 10;          % 参与FFT的基波周期数(整数)
DT    = 1e-6;        % 采样步长 [s], 与原Python一致 (1us)
Kmax  = 4000;        % 最高谐波次数 => 4000*50 = 200 kHz

Rload = 0.5;         % 负载等效电阻 [ohm] (设为0则只算电压谐波)
Lload = 2e-3;        % 负载等效电感 [H]

%% ---------------------------- 2. 谐波分析 ----------------------------
R = hbridge_spwm_harmonics(fc, f0, Vdc, M, ...
        Ncyc=Ncyc, DT=DT, Kmax=Kmax, Rload=Rload, Lload=Lload, Verbose=true);
S = harmonic_report(R, 20); %#ok<NASGU>

%% ---------------------------- 3. 导出数据 ----------------------------
csvname = fullfile(outdir, sprintf('spectrum_fc%g_M%.2f.csv', fc, M));
export_spectrum_csv(R, csvname);
fprintf('\n  [OK] 频谱数据已导出: %s\n', csvname);

%% ---------------------------- 4. 绘图 ----------------------------
f1 = plot_waveforms(R);
export_fig_png(f1, fullfile(outdir, 'fig1_waveform.png'));

f2 = plot_spectrum(R);
export_fig_png(f2, fullfile(outdir, 'fig2_spectrum.png'));

f3 = plot_spectrum_detail(R);
export_fig_png(f3, fullfile(outdir, 'fig3_spectrum_detail.png'));

%% ------------------- 5. 调制比 M 扫描 (0.1 ~ 1.0) -------------------
Mlist = 0.1:0.1:1.0;
Sw = sweep_modulation_index(fc, f0, Vdc, Mlist, ...
        Kmax=min(Kmax, 2000), Ncyc=2, Rload=Rload, Lload=Lload, Verbose=true);
writetable(Sw.tbl, fullfile(outdir, 'thd_vs_M.csv'));

f4 = plot_thd_vs_M(Sw);
export_fig_png(f4, fullfile(outdir, 'fig4_thd_vs_M.png'));

f5 = plot_spectra_vs_M(fc, f0, Vdc, [0.2 0.4 0.6 0.8 1.0], Kmax);
export_fig_png(f5, fullfile(outdir, 'fig5_spectra_vs_M.png'));

%% ---------------------------- 6. 汇总 ----------------------------
fprintf('\n');
fprintf('======================= 汇总 (ASCII) =======================\n');
fprintf(' fc=%.14g Hz  f0=%g Hz  Vdc=%g V  M=%g  Kmax=%d (%.0f kHz)\n', ...
    fc, f0, Vdc, M, Kmax, Kmax*f0/1000);
fprintf(' A1 = %.6f V (theory %.6f V)\n', R.Amp(1), R.A1_theory);
fprintf(' DC = %.3e V   Vrms(wave) = %.4f V   Vrms(spec) = %.4f V\n', ...
    R.DC, R.RMS_time, R.RMS_spec);
fprintf(' THD(<=50) = %.3f %%   THD(<=2fc) = %.3f %%   THD(<=%.0fkHz) = %.3f %%\n', ...
    100*R.THD_50, 100*R.THD_2fc, Kmax*f0/1000, 100*R.THD_all);
fprintf(' WTHD = %.3f %%   THDi(<=50) = %.3f %%\n', 100*R.WTHD, 100*R.THD_i_50);
fprintf(' 1us-FFT max rel err (significant harmonics) = %.3f %%\n', 100*R.Amp_fft_maxerr);
fprintf(' 输出文件:\n');
d = dir(fullfile(outdir, '*'));
for k = 1:numel(d)
    if ~d(k).isdir
        fprintf('   %-40s %8.1f kB\n', d(k).name, d(k).bytes/1024);
    end
end
fprintf('============================================================\n');
