function c = spwm_tri_carrier(t, fc)
%SPWM_TRI_CARRIER  三角载波 (与原Python脚本逐点一致)
%   phase = frac(t*fc);  c = 2*|2*phase-1| - 1
%   峰值 +1 (t = k/fc), 谷值 -1 (t = (k+0.5)/fc), 范围 [-1, +1]
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
ph = t .* fc - floor(t .* fc);
c  = 2 * abs(2*ph - 1) - 1;
end
