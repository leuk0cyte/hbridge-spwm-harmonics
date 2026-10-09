function T = export_spectrum_csv(R, fname)
%EXPORT_SPECTRUM_CSV  把谐波分析结果导出为 CSV 表格
%
%   T = export_spectrum_csv(R, 'spectrum.csv')
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

H   = R.harm;
pct = 100 * R.Amp / max(R.Amp(1), eps);

T = table(H, R.freq, R.Amp, pct, R.Ph_deg, ...
    'VariableNames', {'harmonic_order', 'frequency_Hz', 'amplitude_V', ...
                      'percent_of_fundamental', 'phase_deg'});

if ~isempty(R.Amp_fft)
    T.amplitude_fft_1us_V      = R.Amp_fft;
    T.rel_error_vs_exact       = R.err_rel_fft;
end
if ~isempty(R.Iamp)
    T.current_amplitude_A      = R.Iamp;
end

if nargin >= 2 && ~isempty(fname)
    writetable(T, fname);
end
end
