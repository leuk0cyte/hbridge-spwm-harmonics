function export_fig_png(f, fname, res)
%EXPORT_FIG_PNG  保存图形为 PNG, 并先去坐标区工具栏(避免导出警告/多余图标)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
arguments
    f
    fname (1,:) char
    res (1,1) double = 200
end
ax = findall(f, 'Type', 'axes');
for k = 1:numel(ax)
    try
        ax(k).Toolbar = [];
    catch
    end
end
exportgraphics(f, fname, 'Resolution', res);
end
