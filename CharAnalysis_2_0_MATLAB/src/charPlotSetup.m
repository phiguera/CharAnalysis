function S = charPlotSetup(Charcoal, figPosition)
% charPlotSetup   Return shared layout constants for CharAnalysis figures.
%
%   S = charPlotSetup(Charcoal)
%   S = charPlotSetup(Charcoal, figPosition)
%
%   Returns a struct S with fields used by all CharPlotFig functions:
%     S.wm          - width multiplier scaled to record length
%     S.FS          - base font size
%     S.FW          - font weight string
%     S.zoneText    - cell array of zone label strings
%     S.figPosition - normalized figure position vector [x y w h]
%
%   An optional figPosition argument overrides the default. CharPlotResults
%   passes a progressively offset position so stacked figures are visible;
%   individual CharPlotFig calls use the default.

%% Width multiplier: scales figure width to record length
if max(Charcoal.ybp) > 15000
    S.wm = 0.5  * 2.8e-5;
elseif max(Charcoal.ybp) > 5500
    S.wm = 1.25 * 2.8e-5;
elseif max(Charcoal.ybp) > 2500
    S.wm = 1.5  * 2.8e-5;
else
    S.wm = 5    * 2.8e-5;
end

%% Typography
S.FS = 8;
S.FW = 'bold';

%% Zone labels (supports up to 7 zones)
S.zoneText = {'Zone 1','Zone 2','Zone 3','Zone 4','Zone 5','Zone 6','Zone 7'};

%% Figure position
if nargin < 2 || isempty(figPosition)
    S.figPosition = [0.1013  0.1933  0.8648  0.6943];
else
    S.figPosition = figPosition;
end

end
