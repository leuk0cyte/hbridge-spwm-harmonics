function r = spwm_level_residual(t, fc, f0, M)
%SPWM_LEVEL_RESIDUAL  跳变条件残差 |v(t)| - |c(t)|, 其零点即开关跳变时刻
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
r = abs(M*sin(2*pi*f0*t)) - abs(spwm_tri_carrier(t, fc));
end
