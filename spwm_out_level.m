function lv = spwm_out_level(t, fc, f0, M)
%SPWM_OUT_LEVEL  H桥单极性SPWM输出电压的归一化电平 {-1, 0, +1}
%   左桥臂 Q1 = (v > c), 右桥臂 Q3 = (-v > c), v = M*sin(2*pi*f0*t)
%   Vout/Vdc = Q1 - Q3 = sign(v) .* (|v| > |c|)
%   +1 -> +Vdc, -1 -> -Vdc, 0 -> 续流(零电平)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
v  = M * sin(2*pi*f0*t);
a  = abs(spwm_tri_carrier(t, fc));
lv = zeros(size(t));
lv(v >  a) =  1;
lv(v < -a) = -1;
end
