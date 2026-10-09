function [Q1, Q2, Q3, Q4, dL, dR] = spwm_gate_states(t, fc, f0, M, td)
%SPWM_GATE_STATES  H桥4个IGBT的门极信号 (含固定死区 td)
%
%   [Q1,Q2,Q3,Q4] = spwm_gate_states(t, fc, f0, M, td)
%   [Q1,Q2,Q3,Q4,dL,dR] = spwm_gate_states(...)   额外返回本桥臂是否处于死区
%
%   死区模型 (固定死区, 开通延迟 td / 关断立即):
%     理想门极  Q1=(v>c), Q2=1-Q1, Q3=(-v>c), Q4=1-Q3
%     实际门极  上管本桥臂有跳变后的 td 内不导通 -> 该桥臂两管全部关断
%   等价判据: 器件在 t 时刻导通  <=>  理想门极在区间 [t-td, t] 内恒为 1
%     (故区间 (t-td, t] 内有任何跳变 -> 该桥臂处于死区, 两管均关断)
%   死区期间桥臂电位由续流二极管决定(见 spwm_deadtime_out_level)。
%
%   td = 0 时退化为理想互补驱动, 与原 Python 脚本发电平完全一致。
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    t
    fc (1,1) double {mustBePositive}
    f0 (1,1) double {mustBePositive}
    M  (1,1) double {mustBeNonnegative}
    td (1,1) double {mustBeNonnegative} = 0
end

t  = t(:);
v  = M * sin(2*pi*f0*t);
c  = spwm_tri_carrier(t, fc);
q1 = v >  c;        % 左桥臂上管理想驱动
q3 = -v > c;        % 右桥臂上管理想驱动

if td > 0
    [~, eL, eR] = spwm_ideal_edges(fc, f0, M, min(t) - 2*td, max(t));
    dL = local_edge_in_window(eL, t - td, t) > 0;
    dR = local_edge_in_window(eR, t - td, t) > 0;
else
    dL = false(size(t));
    dR = false(size(t));
end

Q1 = double(~dL &  q1);   Q2 = double(~dL & ~q1);
Q3 = double(~dR &  q3);   Q4 = double(~dR & ~q3);
end

% ------------------------------------------------------------------
function n = local_edge_in_window(e, lo, hi)
% 统计每个 [lo(i), hi(i)] 半开区间 (lo, hi] 内的跳变个数
n = zeros(numel(hi), 1);
if isempty(e)
    return
end
e = sort(e(:));
n = local_cumcount(e, hi) - local_cumcount(e, lo);
end

function C = local_cumcount(e, x)
% C(i) = 跳变时刻中 <= x(i) 的个数
C = discretize(x(:), [-inf; e(:); inf]) - 1;
C(isnan(C)) = 0;
end