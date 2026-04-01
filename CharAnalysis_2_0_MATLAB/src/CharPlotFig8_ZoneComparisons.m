function CharPlotFig8_ZoneComparisons(r)
% CharPlotFig8_ZoneComparisons   Figure 8: Between-zone comparisons of raw CHAR.
%
%   CharPlotFig8_ZoneComparisons(results)
%
%   Left panel:  Cumulative distribution functions (CDFs) of raw CHAR
%                values within each zone. Two-sample KS tests compare
%                zones pairwise; a table of p-values is displayed.
%   Right panel: Box plots of raw CHAR values by zone (10th, 25th, 50th,
%                75th, and 90th percentiles).
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
site         = r.site;

zoneDiv = Pretreatment.zoneDiv;
nZones  = length(zoneDiv) - 1;

S = charPlotSetup(Charcoal);
FS          = S.FS;
zoneText    = S.zoneText;
figPosition = S.figPosition - [0.0165*5 0.0225*5 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 8 - Between-zone comparisons of raw CHAR distributions
%% ------------------------------------------------------------------------
figure(8); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Between-zone comparisons of raw charcoal distributions')

col = [0 0 0; .5 .5 .5; 0 1 0; 0 0 1; 1 0 1; 1 0 0; 1 0.4 0; 0 0.5 0];

Charcoal.accZone = NaN(1000, nZones);
for i = 1:nZones
    CHAR_z = Charcoal.acc(Charcoal.ybp >= zoneDiv(i) & ...
                           Charcoal.ybp <  zoneDiv(i+1));
    Charcoal.accZone(1:length(CHAR_z), i) = CHAR_z;
    subplot(3,5,[6 7 8 11 12 13])
    [f, xc] = ecdf(CHAR_z);
    plot(xc, f, 'color', col(i,:), 'linewidth', 2); hold on
end
set(gca, 'tickdir', 'out')
xlabel('CHAR (pieces cm^-^2 yr^-^1)')
ylabel('cum. prop.')
title('CDFs of zone-specific raw CHAR, with KS test results')
box off
legend(char(zoneText))

if nZones > 1
    h_ks = NaN(nZones-1, nZones-1);
    pKS  = NaN(nZones-1, nZones-1);
    for i = 1:nZones-1
        for j = 2:nZones
            [h_ks(i,j-1), pKS(i,j-1)] = kstest2( ...
                Charcoal.accZone(:,i), Charcoal.accZone(:,j));
        end
    end
    pKSResults = zeros(size(pKS)+1);
    pKSResults(2:end,2:end) = pKS;
    pKSResults(1,2:end) = 2:length(pKS)+1;
    pKSResults(2:end,1) = (1:length(pKS))';
    text(0.5*max(Charcoal.acc), 0.63, 'KS p-value matrix:', 'FontWeight', 'bold')
    text(0.45*max(Charcoal.acc), 0.48, num2str(round(pKSResults*1000)/1000))
    text(0.45*max(Charcoal.acc), 0.58, 'Zone', 'FontWeight', 'Bold', ...
        'BackgroundColor', 'w', 'HorizontalAlignment', 'Center')
end
text(0, 1.5, [char(site), ': Between-zone comparisons of raw CHAR distributions'], ...
    'FontSize', FS+2, 'FontWeight', 'Bold')

subplot(3,5,[9 10 14 15])
boxplot(fliplr(Charcoal.accZone), 'colors', 'kkk', 'symbol', '.k')
set(gca, 'xticklabel', fliplr(1:nZones), 'tickdir', 'out', ...
    'ylim', [min(min(Charcoal.accZone)) max(max(Charcoal.accZone))])
ylabel('CHAR (pieces cm^-^2 yr^-^1)')
xlabel('zone number')
title('Box plots of raw CHAR per zone')
box off

end
