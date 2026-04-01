function CharPlotFig9_ThresholdDetails(r)
% CharPlotFig9_ThresholdDetails   Figure 9: Alternative displays of threshold values.
%
%   CharPlotFig9_ThresholdDetails(results)
%
%   Panel (a): C_int with C_back and C_thresh, and identified peaks.
%   Panel (b): C_peak and C_thresh displayed in the ratio domain.
%   Panel (c): C_peak and C_thresh displayed in the residual domain.
%
%   This figure illustrates how the selected threshold would appear under
%   both C_peak definitions. Produced only when Results.allFigures == 1,
%   but can also be called directly at any time.
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
PeakAnalysis = r.PeakAnalysis;
CharThresh   = r.CharThresh;
Post         = r.Post;

zoneDiv = Pretreatment.zoneDiv;

peakScreenIn      = Post.peakScreenIn;
CharcoalCharPeaks = Post.CharcoalCharPeaks;

S = charPlotSetup(Charcoal);
wm          = S.wm;
FS          = S.FS;
figPosition = S.figPosition - [0.0165*6 0.0225*6 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 9 - Alternative displays of threshold value(s)
%% ------------------------------------------------------------------------
figure(9); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Alternative displays of threshold value(s)')

CHAR_ratio     = Charcoal.accI ./ Charcoal.accIS;
CHAR_residuals = Charcoal.accI  - Charcoal.accIS;
y_peaks = CharcoalCharPeaks(:,end);  y_peaks(y_peaks==0) = -99;

if PeakAnalysis.cPeak ~= -99
    if PeakAnalysis.cPeak == 1
        y1 = Charcoal.accIS + CharThresh.pos(:,end);
        y2 = y1 ./ Charcoal.accIS;
        y3 = CharThresh.pos(:,end);
    else  % ratio
        y1 = Charcoal.accIS .* CharThresh.pos(:,end);
        y2 = CharThresh.pos(:,end);
        y3 = Charcoal.accIS - (Charcoal.accIS ./ CharThresh.pos(:,end));
    end
else
    y1 = CharThresh.pos(:,end);
    y2 = CharThresh.pos(:,end) ./ Charcoal.accIS;
    y3 = CharThresh.pos(:,end) - Charcoal.accIS;
end

% -- Panel (a): C_int with C_back and threshold ---------------------------
subplot(3,1,1)
stairs(Charcoal.ybpI, Charcoal.accI, 'color', [.75 .75 .75]); hold on
plot(Charcoal.ybpI, Charcoal.accIS, 'linewidth', 1.5, 'color', 'k')
plot(Charcoal.ybpI, y1, 'r', 'linewidth', 1.5)
plot(Charcoal.ybpI, max(Charcoal.accI)*0.80*y_peaks, '+k', 'color', 'r')
plot(Charcoal.ybpI(peakScreenIn), ...
    max(Charcoal.accI)*0.80*logical(peakScreenIn), '.', 'color', [.75 .75 .75])
xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)]); ylim([0 max(Charcoal.accI)])
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), ...
    'XTickLabel', 0:max(zoneDiv)/1000, 'position', [0.1 .75 wm*range(zoneDiv) .175])
box off; grid on
ylabel('CHAR (pieces cm^-^2 yr^-^1)')
title('(a) Relationship between C_b_a_c_k_g_r_o_u_n_d and peak threshold')
legend('C_i_n_t_e_r_p_o_l_a_t_e_d', 'C_b_a_c_k_g_r_o_u_n_d', ...
    'C_t_h_r_e_s_h_o_l_d', 'peaks', 'Location', 'Eastoutside')

% -- Panel (b): C_peak as ratio -------------------------------------------
subplot(3,1,2)
stairs(Charcoal.ybpI, CHAR_ratio, 'color', [.75 .75 .75]); hold on
plot(Charcoal.ybpI, y2, 'r', 'linewidth', 1.5)
xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), ...
    'XTickLabel', 0:max(zoneDiv)/1000, 'position', [0.1 .50 wm*range(zoneDiv) .175])
box off; grid on
ylabel({'threshold', '(threshold / low-frequency CHAR)'})
title('(b) Peak threshold as ratios of C_b_a_c_k_g_r_o_u_n_d')
legend('C_p_e_a_k', 'C_t_h_r_e_s_h_o_l_d', 'peaks identified', 'Location', 'Eastoutside')

% -- Panel (c): C_peak as residual ----------------------------------------
subplot(3,1,3)
stairs(Charcoal.ybpI, CHAR_residuals, 'color', [.75 .75 .75]); hold on
plot(Charcoal.ybpI, y3, 'r', 'linewidth', 1.5)
xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
if min(y3) <= 0; ylim([0 abs(max(y3)*2)]); else; ylim([0 max(y3)*2]); end
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), ...
    'XTickLabel', 0:max(zoneDiv)/1000, 'position', [0.1 .20 wm*range(zoneDiv) .175])
box off; grid on
xlabel('time (cal. yr BP x 1000)')
ylabel({'threshold', '(threshold - low-frequency CHAR)'})
title('(c) Peak threshold as residuals of C_b_a_c_k_g_r_o_u_n_d')
legend('C_p_e_a_k', 'C_t_h_r_e_s_h_o_l_d', 'peaks identified', 'Location', 'Eastoutside')

end
