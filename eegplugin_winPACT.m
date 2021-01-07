% eegplugin_winPACT() - Analyze phase-amplitude coupling.

% Copyright (C) 2017 Makoto Miyakoshi, SCCN, INC, UCSD.
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.
%
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with this program; if not, write to the Free Software
% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1.07  USA

function vers = eegplugin_winPACT(fig, trystrs, catchstrs)
    
vers = 'beta';
if nargin < 3
    error('eegplugin_winPACT requires 3 arguments');
end

% Create a highLevelMenu.
highLevelManu = findobj(fig, 'tag', 'tools');
submenu       = uimenu(highLevelManu, 'label', 'winPACT','separator','on');

% Add submenus.
toolsmenu = findobj(fig, 'tag', 'tools');
uimenu( submenu, 'label', '(Optimize)', 'callback', 'winPACT_optimize');
uimenu( submenu, 'label', '(Simulate)', 'callback', 'winPACT_simulate');
uimenu( submenu, 'label', 'Precompute', 'callback', 'winPACT_precompute');
uimenu( submenu, 'label', 'Visualize',  'callback', 'winPACT_visualize');
