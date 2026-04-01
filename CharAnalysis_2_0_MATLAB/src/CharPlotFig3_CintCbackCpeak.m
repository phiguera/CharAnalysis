function CharPlotFig3_CintCbackCpeak(r)
% CharPlotFig3_CintCbackCpeak   Figure 3: C_interpolated, C_background, and C_peak.
%
%   CharPlotFig3_CintCbackCpeak(results)
%
%   Panel (a): C_interpolated with C_background overlaid.
%   Panel (b): C_peak with positive and negative threshold values, and
%              identified peaks marked as + symbols. Peaks failing the
%              minimum-count criterion are shown as grey dots.
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
PeakAnalysis = r.PeakAnalysis;
CharThresh   = r.CharThresh;
Smoothing    = r.Smoothing;
site         = r.site;
Results      = r.Results;
gapIn        = r.gapIn;
Post         = r.Post;

zoneDiv   = Pretreatment.zoneDiv;
yrInterp  = Pretreatment.yrInterp;
yr        = Smoothing.yr;
transform = Pretreatment.transform;
nGaps     = size(gapIn, 1);

peakIn       = Post.peakIn;
peakScreenIn = Post.peakScreenIn;

S = charPlotSetup(Charcoal);
wm          = S.wm;
FS          = S.FS;
FW          = S.FW;
zoneText    = S.zoneText;
figPosition = S.figPosition;

%% ------------------------------------------------------------------------
%% FIGURE 3 - C_interpolated, C_background, C_peak
%% ------------------------------------------------------------------------
figure(3); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Resampled charcoal (C_interpolated), low-frequency trends (C_background), and detrended series (C_peak)')

% -- Panel (a): C_int and C_background -------------------------------------
subplot(3,5,1:4)
x  = Charcoal.ybpI;
y  = Charcoal.accI;
y2 = Charcoal.accIS;
H  = bar(x, y); hold on
set(H, 'facecolor', [0 0 0], 'edgecolor', [0 0 0])
plot(x, y2, 'color', [.5 .5 .5], 'linewidth', 2)
xlim([zoneDiv(1), zoneDiv(end)])
ylim([0, 1.1*max(y)])
if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv], [max(y)*1.01 max(y)*1.1], '-k', ...
        'color', [.5 .5 .5], 'linewidth', 2)
end
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', ...
    'YMinorTick', 'on', 'TickDir', 'out', ...
    'XTick', 0:1000:max(zoneDiv), 'XTickLabel', [], ...
    'Position', [0.1 .55 wm*range(zoneDiv) .175])
box off
ylabel(charYLabel(transform), 'FontSize', FS)
for z = 2:length(zoneDiv)
    if z < length(zoneDiv)
        plot([zoneDiv(z) zoneDiv(z)], [max(y)*1.01 max(y)*1.1], ...
            '-k', 'color', [.5 .5 .5], 'linewidth', 2)
    end
    if length(zoneDiv) > 2
        text(mean(zoneDiv(z-1:z)), max(y)*1.05, zoneText(z-1), ...
            'horizontalalignment', 'center', 'fontweight', 'normal', 'fontsize', FS)
    end
end
title({['(a) ', char(site), ': C_i_n_t_e_r_p_o_l_a_t_e_d (', num2str(yrInterp), ' yr)', ...
    ' and C_b_a_c_k_g_r_o_u_n_d defined by ', ...
    num2str(yr(end)), '-yr trends']}, 'FontSize', FS, 'FontWeight', FW)

% -- Panel (b): C_peak and thresholds -------------------------------------
subplot(3,5,10)
y  = Charcoal.peak;
y2 = CharThresh.pos(:,end);
y3 = CharThresh.neg(:,end);
if PeakAnalysis.cPeak == 1
    H1 = bar(x, y, 1);
else
    H1 = bar(x, y, 1, 'BaseValue', 1);
end
set(H1, 'facecolor', [0 0 0], 'edgecolor', [0 0 0])
hold on
plot(x, y2, 'r'); plot(x, y3, 'r')
if PeakAnalysis.cPeak == 1
    plot([min(x) max(x)], [0 0], 'k')
elseif PeakAnalysis.cPeak == 2
    plot([min(x) max(x)], [1 1], 'k')
end
plot(x(peakScreenIn), ones(size(peakScreenIn))*0.8*max(y), '.', ...
    'color', [.75 .75 .75])
plot(x(peakIn), 0.8*max(y), '+k')
xlim([zoneDiv(1), zoneDiv(end)])
ylim([min(y3), 1.1*max(y)])
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Position', [0.1 .32 wm*range(zoneDiv) .175], 'TickDir', 'out', ...
    'XTick', 0:1000:max(zoneDiv), 'XTickLabel', 0:max(zoneDiv)/1000)
if PeakAnalysis.cPeak == 1
    title({['(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d - C_b_a_c_k_g_r_o_u_n_d), ' ...
        'thresholds defining C_n_o_i_s_e, and peaks'], ...
        'identified (gray dots fail to pass peak-magnitude test)'}, ...
        'VerticalAlignment', 'middle', 'FontSize', FS, 'FontWeight', 'bold')
elseif PeakAnalysis.cPeak == 2
    title(['(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d / C_b_a_c_k_g_r_o_u_n_d), ' ...
        'thresholds defining C_n_o_i_s_e, and peaks identified ' ...
        '(gray peaks fail to pass peak-magnitude test)'], ...
        'VerticalAlignment', 'middle', 'FontSize', FS, 'FontWeight', 'bold')
else
    title(['(b) C_i_n_t_e_r_p_o_l_a_t_e_d, thresholds defining C_n_o_i_s_e, ' ...
        'and peaks identified (gray peaks fail to pass peak-magnitude test)'], ...
        'VerticalAlignment', 'middle', 'FontSize', FS, 'FontWeight', 'bold')
end
xlabel('time (cal. yr BP x 1000)', 'FontSize', FS)
ylabel(charYLabel(transform), 'FontSize', FS)
box off

if nGaps > 0
    charDrawGaps(gca, gapIn, Charcoal)
end

end
