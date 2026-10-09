function set_cjk_font()
%SET_CJK_FONT  让图形中的中文正常显示(Windows 下的微软雅黑)
%
% SPDX-License-Identifier: GPL-3.0-or-later
% Copyright (C) 2026 leuk0cyte
% 本文件是 hbridge-spwm-harmonics 的一部分, 以 GPLv3 或更新版本发布, 无任何担保, 详见 LICENSE。
try
    set(groot, 'defaultAxesFontName',   'Microsoft YaHei');
    set(groot, 'defaultTextFontName',   'Microsoft YaHei');
    set(groot, 'defaultLegendFontName', 'Microsoft YaHei');
    set(groot, 'defaultAxesFontSize',   10);
catch
end
end
