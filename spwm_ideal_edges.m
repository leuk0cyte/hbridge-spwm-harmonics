function [eAll, eL, eR] = spwm_ideal_edges(fc, f0, M, t1, t2)
%SPWM_IDEAL_EDGES  理想开关跳变时刻 (|v|=|c|), 并按左/右桥臂分类
%
%   [eAll, eL, eR] = spwm_ideal_edges(fc, f0, M)          默认区间 [0, 1/f0]
%   [eAll, eL, eR] = spwm_ideal_edges(fc, f0, M, t1, t2)  指定扫描区间
%
%   输出:
%     eAll : 全部跳变时刻, 即输出电压电平发生变化(|v|=|c|)的时刻
%     eL   : 左桥臂上管 Q1 门极翻转时刻 (v = c 处)
%     eR   : 右桥臂上管 Q3 门极翻转时刻 (v = -c 处)
%
%   调制与载波定义 (与原 Python 脚本一致):
%     v(t) = M*sin(2*pi*f0*t),  c(t) = spwm_tri_carrier(t, fc)
%     Q1 = (v > c)   (左桥臂上管),   Q3 = (-v > c)  (右桥臂上管)
%   注意: 两个桥臂在**不同时刻**翻转 —— Q1 在 v=c 处, Q3 在 v=-c 处;
%         两者仅在 c=0 时才同时翻转, 因此 eL 与 eR 通常是 eAll 的不同子集。
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。

arguments
    fc (1,1) double {mustBePositive}
    f0 (1,1) double {mustBePositive}
    M  (1,1) double {mustBeNonnegative}
    t1 (1,1) double = 0
    t2 (1,1) double = 0
end

if t2 <= t1
    t2 = t1 + 1/f0;
end
eAll = zeros(0,1);  eL = zeros(0,1);  eR = zeros(0,1);
if M <= 0
    return
end

% ---- 1) 粗扫描: 找出"活动指示"发生变化的区间 ----
dtScan = max(1/(1000*fc), (t2 - t1)/4e6);
ts = (t1:dtScan:t2).';
if ts(end) < t2 - eps(t2)
    ts(end+1) = t2;
end
act = abs(M*sin(2*pi*f0*ts)) > abs(spwm_tri_carrier(ts, fc));
chg = find(act(1:end-1) ~= act(2:end));
if isempty(chg)
    return
end

% ---- 2) 二分法细化到机器精度 ----
eAll = local_refine_roots(@(t) abs(M*sin(2*pi*f0*t)) - abs(spwm_tri_carrier(t, fc)), ...
                          ts(chg), ts(chg+1));
eAll = sort(eAll(:));

% ---- 3) 按桥臂分类: 看该时刻前后哪个上管的驱动翻转了 ----
ep = max(1e-12, dtScan/100);
v1 = M*sin(2*pi*f0*(eAll-ep));  c1 = spwm_tri_carrier(eAll-ep, fc);
v2 = M*sin(2*pi*f0*(eAll+ep));  c2 = spwm_tri_carrier(eAll+ep, fc);
eL = eAll( (v1 > c1) ~= (v2 > c2) );
eR = eAll( (-v1 > c1) ~= (-v2 > c2) );
end

% ------------------------------------------------------------------
function tr = local_refine_roots(F, a, b)
% 向量化二分法 (要求两端点函数值异号); 端点同号时退化为局部极小定位
fa = F(a); fb = F(b);
ok = (fa .* fb) < 0;
tr = 0.5*(a + b);
if any(ok)
    lo = a(ok); hi = b(ok); flo = fa(ok);
    for it = 1:80
        mid  = 0.5*(lo + hi);
        fm   = F(mid);
        same = (flo .* fm) > 0;
        lo(same) = mid(same); flo(same) = fm(same);
        hi(~same) = mid(~same);
    end
    tr(ok) = 0.5*(lo + hi);
end
for k = find(~ok).'
    tr(k) = fminbnd(@(t) abs(F(t)), a(k), b(k));
end
end