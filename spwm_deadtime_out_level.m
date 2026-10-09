function [lv, levL, levR, dL, dR] = spwm_deadtime_out_level(t, fc, f0, M, td, phi)
%SPWM_DEADTIME_OUT_LEVEL  含死区的H桥输出电压归一化电平 {-1, 0, +1}
%
%   lv = spwm_deadtime_out_level(t, fc, f0, M, td)
%   lv = spwm_deadtime_out_level(t, fc, f0, M, td, phi)
%   [lv, levL, levR, dL, dR] = spwm_deadtime_out_level(...)
%
%   td  : 固定死区时间 [s] (td=0 时与理想单极性SPWM完全一致)
%   phi : 负载电流滞后角 [rad]; 电流 i(t) ∝ sin(2*pi*f0*t - phi)
%         也可以直接给一个函数句柄 phi = @(t) i(t), 由它决定电流符号
%
%   死区期间的桥臂电位由续流二极管钳位(与电流方向有关):
%     左桥臂: i > 0 (电流流出左桥臂中点) -> 下管二极管续流 -> 电平 0
%             i < 0                        -> 上管二极管续流 -> 电平 1
%     右桥臂: i > 0 (电流流入右桥臂中点) -> 上管二极管续流 -> 电平 1
%             i < 0                        -> 下管二极管续流 -> 电平 0
%   由此每个桥臂每次开关在载波周期内的平均电压误差为
%     dV = -sign(i) * td * fc * Vdc
%   单极性调制下**两个桥臂都在斩波**(只是在不同时刻翻转), 两者误差叠加,
%   故整机每个载波周期的平均误差为
%     dV_out = -2 * sign(i) * td * fc * Vdc
%   (t_d=4us, fc=2400Hz, Vdc=100V 时为 ±1.92 V, 已数值验证)
%   另: 若某个脉冲宽度 < td, 该管永不导通, 电平被二极管钳位(窄脉冲丢失),
%   这正是死区造成低次谐波与基波幅值下降(本例 4us 时约 -1.9%)的根源。
%
%   输出电压 Vout = Vdc * lv = Vdc * (levL - levR)。
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    t
    fc  (1,1) double {mustBePositive}
    f0  (1,1) double {mustBePositive}
    M   (1,1) double {mustBeNonnegative}
    td  (1,1) double {mustBeNonnegative} = 0
    phi = 0
end

t = t(:);

% ---- 实际门极(含死区) ----
[Q1, ~, Q3, ~, dL, dR] = spwm_gate_states(t, fc, f0, M, td);
levL = double(Q1);
levR = double(Q3);

% ---- 死区期间由二极管钳位 ----
if td > 0
    if isa(phi, 'function_handle')
        ipos = phi(t) > 0;                  % 由给定的电流函数决定方向
    else
        ipos = sin(2*pi*f0*t - phi) > 0;    % 基波电流方向
    end
    levL(dL) = double(~ipos(dL));           % 左桥臂: i>0 -> 0, i<0 -> 1
    levR(dR) = double( ipos(dR));           % 右桥臂: i>0 -> 1, i<0 -> 0
end

lv = levL - levR;
end