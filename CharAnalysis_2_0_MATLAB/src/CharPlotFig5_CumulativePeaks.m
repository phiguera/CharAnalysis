function CharPlotFig5_CumulativePeaks(r)
% CharPlotFig5_CumulativePeaks   Figure 5: Cumulative peaks through time.
%
%   CharPlotFig5_CumulativePeaks(results)
%
%   Cumulative sum of identified peaks as a function of time. The slope
%   at any point is the instantaneous fire frequency (fires yr^-1).
%   Areas with missing values are indicated by grey boxes.
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
site         = r.site;
gapIn        = r.gapIn;
Post         = r.Post;

zoneDiv = Pretreatment.zoneDiv;
nGaps   = size(gapIn, 1);

CharcoalCharPeaks = Post.CharcoalCharPeaks;

S = charPlotSetup(Charcoal);
wm          = S.wm;
FS          = S.FS;
figPosition = S.figPosition - [0.0165*2 0.0225*2 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 5 - Cumulative peaks through time
%% ------------------------------------------------------------------------
figure(5); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Cumulative peaks through time')

xPlot = Charcoal.ybpI(CharcoalCharPeaks(:,end) > 0);
yPlot = flipud(cumsum(logical(find(CharcoalCharPeaks(:,end)>0))));

if ~isempty(xPlot)
    plot(xPlot, yPlot, '.k')
    set(gca, 'XDir', 'reverse', 'TickDir', 'out', 'XMinorTick', 'on', ...
        'YMinorTick', 'on', 'FontSize', FS, ...
        'XTick', 0:1000:max(zoneDiv), 'XTickLabel', 0:max(zoneDiv)/1000, ...
        'Position', [.1 .1 0.8*wm*range(zoneDiv) max(yPlot)*0.005])
    if nGaps > 0
        charDrawGaps(gca, gapIn, Charcoal)
        legend('Cumulative peaks', 'missing values')
    end
    xlabel('time (cal. yr BP x 1000)', 'FontSize', FS)
    ylabel('cumulative number of peaks', 'FontSize', FS)
    xlim([min(Charcoal.ybp), max(zoneDiv)])
    title([char(site), ': Cumulative fires as a function of time'], ...
        'FontSize', FS, 'FontWeight', 'bold')
    box off; grid on
end

end
