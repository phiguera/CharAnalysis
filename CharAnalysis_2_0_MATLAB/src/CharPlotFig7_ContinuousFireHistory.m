function CharPlotFig7_ContinuousFireHistory(r)
% CharPlotFig7_ContinuousFireHistory   Figure 7: Continuous fire history.
%
%   CharPlotFig7_ContinuousFireHistory(results)
%
%   Top panel:    Peak magnitude as bars with identified peaks as + symbols.
%   Middle panel: Fire return intervals and smoothed FRI curve with CIs.
%   Bottom panel: Smoothed fire frequency.
%
%   In all panels, areas with missing values are indicated by grey boxes.
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
PeakAnalysis = r.PeakAnalysis;
site         = r.site;
gapIn        = r.gapIn;
Post         = r.Post;

zoneDiv = Pretreatment.zoneDiv;
nGaps   = size(gapIn, 1);

peakIn       = Post.peakIn;
peakScreenIn = Post.peakScreenIn;
peak_mag     = Post.peak_mag;
ff_sm        = Post.ff_sm;
FRIyr        = Post.FRIyr;
FRI          = Post.FRI;
smFRIyr      = Post.smFRIyr;
smFRI        = Post.smFRI;
smFRIci      = Post.smFRIci;
alpha        = Post.alpha;

S = charPlotSetup(Charcoal);
wm          = S.wm;
FS          = S.FS;
zoneText    = S.zoneText;
figPosition = S.figPosition - [0.0165*4 0.0225*4 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 7 - Continuous fire history
%% ------------------------------------------------------------------------
figure(7); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Continuous fire history: peak magnitude, FRIs through time, and smoothed fire frequency')

% -- Panel: Peak magnitude -------------------------------------------------
subplot(3,5,1:4); hold off
x = Charcoal.ybpI;
H = bar(x, peak_mag);
set(H, 'FaceColor', 'k', 'BarWidth', 1.0); hold on
plot(x(peakScreenIn), ones(size(peakScreenIn))*0.8*max(peak_mag), ...
    '.', 'color', [.75 .75 .75])
plot(x(peakIn), 0.8*max(peak_mag), '+r')
y_lim7 = [0, 1.1*max(peak_mag)];
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), 'XTickLabel', [], ...
    'XLim', [zoneDiv(1) zoneDiv(end)], 'Ylim', y_lim7, 'box', 'off', ...
    'position', [0.1 .75 wm*range(zoneDiv) .12])
if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)], [max(y_lim7)*0.90 max(y_lim7)], ...
                '-k', 'color', [.5 .5 .5], 'linewidth', 2)
        end
        text(mean(zoneDiv(z-1:z)), max(y_lim7)*0.95, char(zoneText(z-1)), ...
            'horizontalalignment', 'center', 'fontweight', 'normal', 'fontsize', FS)
    end
end
if nGaps > 0; charDrawGaps(gca, gapIn, Charcoal); end
ylabel({'peak magnitude', '(pieces cm^-^2 peak^-^1)'}, 'FontSize', FS)
title({'Peak magnitude, FRIs, and fire frequ.', '', char(site), ''}, ...
    'FontSize', FS, 'FontWeight', 'Bold', 'VerticalAlignment', 'Bottom')

% -- Panel: Smoothed fire frequency ----------------------------------------
subplot(3,5,9); hold off
plot(Charcoal.ybpI, ff_sm, 'k', 'linewidth', 1); hold on
y_lim9 = [0, max(1.1*max(ff_sm), eps)];
if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)], [max(y_lim9)*0.90 max(y_lim9)], ...
                '-k', 'color', [.5 .5 .5], 'linewidth', 2)
        end
    end
end
if nGaps > 0; charDrawGaps(gca, gapIn, Charcoal); end
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), ...
    'XTickLabel', 0:max(zoneDiv)/1000, 'XLim', [zoneDiv(1) zoneDiv(end)], ...
    'box', 'off', 'ylim', y_lim9, 'position', [0.1 .45 wm*range(zoneDiv) .12])
ylabel({'fire frequency', ['(fires ', num2str(PeakAnalysis.peakFrequ), ' yr^-^1)']})
xlabel('time (cal. yr BP x 1000)', 'FontSize', FS)

% -- Panel: FRIs through time ----------------------------------------------
subplot(3,5,11:14); hold off
if isnumeric(FRIyr) && ~isnan(FRIyr(1)) && length(smFRI) > 2
    x_lim7 = [0, 1.1*max(FRI)];
    X = [smFRIyr'; flipud(smFRIyr')];
    Y = [smFRIci(:,1); flipud(smFRIci(:,2))];
    plot(FRIyr, FRI, 'sk', 'MarkerFaceColor', [.75 .75 .75], 'MarkerSize', 4); hold on
    plot(smFRIyr, smFRI, '-k')
    Hf = fill(X, Y, [0.8 0.8 0.8]);
    set(Hf, 'edgecolor', [0.8 0.8 0.8])
    plot(FRIyr, FRI, 'sk', 'MarkerFaceColor', [.75 .75 .75], 'MarkerSize', 4)
    plot(smFRIyr, smFRI, '-k')
    set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
        'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), 'Ylim', x_lim7, ...
        'XTickLabel', [], 'XLim', [zoneDiv(1) zoneDiv(end)], 'box', 'off', ...
        'position', [0.1 .6 wm*range(zoneDiv) .12], 'yScale', 'linear')
    if nGaps > 0; charDrawGaps(gca, gapIn, Charcoal); end
    if length(zoneDiv) > 2
        for z = 2:length(zoneDiv)
            if z < length(zoneDiv)
                plot([zoneDiv(z) zoneDiv(z)], [max(x_lim7)*0.90 max(x_lim7)], ...
                    '-k', 'color', [.5 .5 .5], 'linewidth', 2)
            end
        end
    end
    ylabel({['FRI (yr fire^-^1)'], ...
        [num2str(PeakAnalysis.peakFrequ), '-yr mean'], ...
        [num2str((1-alpha)*100), '% CI']}, 'FontSize', FS)
    grid on
end

end
