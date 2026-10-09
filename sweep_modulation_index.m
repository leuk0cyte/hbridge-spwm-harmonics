function S = sweep_modulation_index(fc, f0, Vdc, Mlist, opts)
%SWEEP_MODULATION_INDEX  扫描调制比 M, 统计 THD / WTHD / 基波幅值
%
%   S = sweep_modulation_index(2400, 50, 100, 0.1:0.1:1.0, Kmax=2000, Rload=0.5, Lload=2e-3)
%
%   S.tbl   : 结果表格 (可直接 writetable 导出)
%   S.M, S.A1, S.THD_50, S.THD_2fc, S.THD_all, S.WTHD, S.THD_i_50 ...
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    fc    (1,1) double {mustBePositive} = 2400
    f0    (1,1) double {mustBePositive} = 50
    Vdc   (1,1) double {mustBePositive} = 100
    Mlist (1,:) double {mustBeNonnegative} = 0.1:0.1:1.0
    opts.Kmax  (1,1) double {mustBePositive, mustBeInteger} = 2000
    opts.Ncyc  (1,1) double {mustBePositive, mustBeInteger} = 2
    opts.Rload (1,1) double {mustBeNonnegative} = 0
    opts.Lload (1,1) double {mustBeNonnegative} = 0
    opts.Verbose (1,1) logical = true
end

n = numel(Mlist);
S = struct();
S.M        = Mlist(:);
S.A1       = nan(n,1);
S.A1_theory= nan(n,1);
S.THD_50   = nan(n,1);
S.THD_2fc  = nan(n,1);
S.THD_4fc  = nan(n,1);
S.THD_all  = nan(n,1);
S.WTHD     = nan(n,1);
S.RMS_total= nan(n,1);
S.THD_i_50 = nan(n,1);
S.THD_i_all= nan(n,1);
S.fc = fc; S.f0 = f0; S.Vdc = Vdc; S.Kmax = opts.Kmax;

if opts.Verbose
    fprintf('\n--- 调制比扫描 (fc=%g Hz, f0=%g Hz, Vdc=%g V, Kmax=%d => %g kHz) ---\n', ...
        fc, f0, Vdc, opts.Kmax, opts.Kmax*f0/1000);
    fprintf('   %5s %10s %10s %10s %10s %10s %10s %10s\n', ...
        'M', 'A1/V', 'THD2fc/%', 'THD4fc/%', 'THDmax/%', 'WTHD/%', 'THDi50/%', 'THDimax/%');
end

for i = 1:n
    R = hbridge_spwm_harmonics(fc, f0, Vdc, Mlist(i), ...
        Ncyc=opts.Ncyc, Kmax=opts.Kmax, Rload=opts.Rload, Lload=opts.Lload, ...
        UseSampledFFT=false, Verbose=false);
    S.A1(i)        = R.Amp(1);
    S.A1_theory(i) = R.A1_theory;
    S.THD_50(i)    = R.THD_50;
    S.THD_2fc(i)   = R.THD_2fc;
    S.THD_4fc(i)   = R.THD_4fc;
    S.THD_all(i)   = R.THD_all;
    S.WTHD(i)      = R.WTHD;
    S.RMS_total(i) = R.RMS_time;
    S.THD_i_50(i)  = R.THD_i_50;
    S.THD_i_all(i) = R.THD_i_all;
    if opts.Verbose
        fprintf('   %5.2f %10.4f %10.3f %10.3f %10.3f %10.3f %10.3f %10.3f\n', ...
            Mlist(i), S.A1(i), 100*S.THD_2fc(i), 100*S.THD_4fc(i), ...
            100*S.THD_all(i), 100*S.WTHD(i), 100*S.THD_i_50(i), 100*S.THD_i_all(i));
    end
end

S.tbl = table(S.M, S.A1, S.A1_theory, ...
    100*S.THD_50, 100*S.THD_2fc, 100*S.THD_4fc, 100*S.THD_all, 100*S.WTHD, ...
    100*S.THD_i_50, 100*S.THD_i_all, S.RMS_total, ...
    'VariableNames', {'M', 'A1_V', 'A1_theory_V', 'THD_50_percent', ...
    'THD_2fc_percent', 'THD_4fc_percent', 'THD_all_percent', 'WTHD_percent', ...
    'THD_current_50_percent', 'THD_current_all_percent', 'Vrms_V'});
end
