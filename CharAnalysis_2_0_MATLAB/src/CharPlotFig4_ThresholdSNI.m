function CharPlotFig4_ThresholdSNI(r)
% CharPlotFig4_ThresholdSNI   Figure 4: Sensitivity to alternative thresholds and SNI.
%
%   CharPlotFig4_ThresholdSNI(results)
%
%   Panel (a): C_int with C_back and all threshold values; peaks identified
%              with the final threshold marked as +.
%   Panel (b): For global threshold -- C_peak distribution with noise PDF
%              and peak count curve. For local threshold -- mean FRI and
%              95% CI per zone for each of the three threshold values.
%   Panel (c): SNI time series.
%   Panel (d): Global SNI boxplot.
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
transform = Pretreatment.transform;

CharcoalCharPeaks = Post.CharcoalCharPeaks;

S = charPlotSetup(Charcoal);
wm          = S.wm;
FS          = S.FS;
FW          = S.FW;
zoneText    = S.zoneText;
figPosition = S.figPosition - [0.0165 0.0225 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 4 - Sensitivity to alternative thresholds and quality of record
%% ------------------------------------------------------------------------
figure(4); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Sensitivity to alternative thresholds and quality of record')

% -- Panel (a): C_int with C_back and threshold ----------------------------
subplot(3,5,1:4)
x  = Charcoal.ybpI;
y  = Charcoal.accI;
y2 = Charcoal.accIS;
y3 = CharThresh.pos(:,end);
y4 = CharThresh.neg(:,end);
y5 = CharcoalCharPeaks;  y5(y5==0) = -99;
H  = bar(x, y, 1); hold on
set(H, 'facecolor', [0 0 0], 'edgecolor', [0 0 0])
plot(x, y2, 'color', [.5 .5 .5], 'linewidth', 2)
if PeakAnalysis.cPeak == 1
    plot(x, y2+y3, 'r');  plot(x, y2+y4, 'r')
else
    plot(x, y2.*y3, 'r'); plot(x, y2.*y4, 'r')
end
xlim([zoneDiv(1), zoneDiv(end)]); ylim([0, 1.15*max(y)])
if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv], [max(y)*1.01 max(y)*1.1], '-k', 'color', [.5 .5 .5])
end
for z = 2:length(zoneDiv)
    if z < length(zoneDiv)
        plot([zoneDiv(z) zoneDiv(z)], [max(y)*1.01 max(y)*1.1], ...
            '-k', 'color', [.5 .5 .5], 'linewidth', 2)
    end
    if length(zoneDiv) > 2
        text(mean(zoneDiv(z-1:z)), max(y)*1.05, char(zoneText(z-1)), ...
            'horizontalalignment', 'center', 'fontweight', 'normal', 'fontsize', FS)
    end
end
plot(x, max(y)*0.78*y5(:,1), '.k', 'color', [.5 .5 .5])
plot(x, max(y)*0.85*y5(:,2), '.k', 'color', [.5 .5 .5])
plot(x, max(y)*0.92*y5(:,3), '.k', 'color', [.5 .5 .5])
if     sum(y5(:,1) - y5(:,end)) == 0;  mIndex = 0.78;
elseif sum(y5(:,3) - y5(:,end)) == 0;  mIndex = 0.92;
else;                                   mIndex = 0.85;
end
plot(x, max(y)*mIndex*y5(:,end), '.k', 'color', [1 1 1])
plot(x, max(y)*mIndex*y5(:,end), '+k')
ylabel(charYLabel(transform), 'FontSize', FS)
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'TickDir', 'out', 'XTick', 0:1000:max(zoneDiv), 'XTickLabel', [], ...
    'Position', [0.1 .55 wm*range(zoneDiv) .175])
title([{char(site)}, {''}, {'(a) C_i_n_t_e_r_p_o_l_a_t_e_d, C_b_a_c_k_g_r_o_u_n_d, and peak ID (+)'}], ...
    'FontSize', FS, 'FontWeight', FW)
box off

% -- Panel (b): Threshold sensitivity -------------------------------------
subplot(3,5,5)
if PeakAnalysis.threshType == 1        % Global threshold
    cpk = Charcoal.peak;
    y3tot = zeros(size(CharThresh.possible));
    for bi = 1:length(CharThresh.possible)
        y3tot(bi) = sum(Charcoal.peak > CharThresh.possible(bi));
    end
    if PeakAnalysis.cPeak == 1
        xPoss = CharThresh.possible;
        xPos2 = xPoss(xPoss>0); %#ok<NASGU>
    else
        xPoss = CharThresh.possible;
        xPos2 = xPoss(xPoss>1); %#ok<NASGU>
    end
    [n, xh] = histcounts(cpk, xPoss, 'Normalization', 'probability');
    xh_c = xh(1:end-1) + diff(xh)/2;
    H1 = bar(xh_c, n, 'BarWidth', 1.0);
    set(H1, 'FaceColor', [.75 .75 .75], 'EdgeColor', [.75 .75 .75])
    hold on
    yyaxis left
    if PeakAnalysis.threshMethod > 1
        y2pdf = CharThresh.noisePDF;
        y2pdf_plot = interp1(CharThresh.possible, y2pdf, xh_c, 'linear', 'extrap');
        plot(xh_c, y2pdf_plot * mean(diff(xh_c)), 'k--', 'linewidth', 1.5)
    else
        plot(-99, -99, 'k--')
    end
    ylabel('relative frequency', 'fontsize', FS)
    ylim([0, 1.01*max(n)])
    set(gca, 'Ycolor', 'k')
    yyaxis right
    if PeakAnalysis.cPeak == 1
        xPlotThresh = xh_c(xh_c > 0);
        yPlotThresh = interp1(CharThresh.possible, y3tot, xPlotThresh, ...
                              'linear', 'extrap');
    else
        xPlotThresh = xh_c(xh_c >= 1);
        yPlotThresh = interp1(CharThresh.possible, y3tot, xPlotThresh, ...
                              'linear', 'extrap');
    end
    plot(xPlotThresh, yPlotThresh, 'k-', 'linewidth', 1.5)
    ylabel('# of peaks identified', 'fontsize', FS, 'rotation', 270, ...
        'verticalalignment', 'bottom', 'color', 'k')
    ylim([0, 0.99*max(y3tot)])
    set(gca, 'Ycolor', 'k')
    yyaxis left
    for t = 1:4
        plot([CharThresh.pos(1,t) CharThresh.pos(1,t)], [0 max(n)], '--', ...
            'color', [0.75 0.75 0.75])
        if t == 4
            plot([CharThresh.pos(1,t) CharThresh.pos(1,t)], [0 max(n)], '-k')
            text(CharThresh.pos(1,t), 0, '<', ...
                'FontSize', FS+2, 'Rotation', 90, 'FontWeight', 'bold')
        end
    end
    if PeakAnalysis.threshMethod > 1
        text(CharThresh.pos(1,end)*1.25, 0.9*max(n), ...
            {[num2str(PeakAnalysis.threshValues(end)*100), '^t^h percentile'], ...
             ['= ', num2str(round(max(xh_c(xh_c<=CharThresh.pos(1,end)))*1e4)/1e4)]}, ...
            'FontSize', FS, 'HorizontalAlignment', 'left', 'backgroundcolor', 'w')
    else
        text(CharThresh.pos(1,end)*1.25, 0.9*max(n), ...
            {'threshold value', ...
             ['= ', num2str(round(max(xh_c(xh_c<=CharThresh.pos(1,end)))*1e4)/1e4)]}, ...
            'FontSize', FS, 'HorizontalAlignment', 'left', 'backgroundcolor', 'w')
    end
    set(gca, 'TickDir', 'out', 'FontSize', FS, 'XMinorTick', 'on', ...
        'XLim', [min(Charcoal.peak) 0.75*max(Charcoal.peak)], ...
        'position', [0.17+wm*range(zoneDiv) .55 .16 .175])
    box off
    if PeakAnalysis.cPeak == 1
        xlabel('residual CHAR value (pieces cm^-^2 yr^-^1)', 'FontSize', FS)
    else
        xlabel('CHAR ratio (pieces cm^-^2 yr^-^1)', 'FontSize', FS)
    end
    if PeakAnalysis.threshMethod > 1
        title({'(b) High-frequency CHAR distribution,', 'est. noise dist., and threshold values'}, 'FontWeight', FW)
    else
        title({'(b) High-frequency CHAR distribution', 'and threshold values'}, 'FontWeight', FW)
    end

else    % Local threshold: mFRI sensitivity by zone
    plot_in = fliplr(1:length(zoneDiv)-1);
    for i = 1:length(plot_in)
        xPlot1 = Charcoal.ybpI(CharcoalCharPeaks(:,1)>0);
        xPlot2 = Charcoal.ybpI(CharcoalCharPeaks(:,2)>0);
        xPlot3 = Charcoal.ybpI(CharcoalCharPeaks(:,3)>0);
        xPlot1 = xPlot1(xPlot1>=zoneDiv(i) & xPlot1<zoneDiv(i+1));
        xPlot2 = xPlot2(xPlot2>=zoneDiv(i) & xPlot2<zoneDiv(i+1));
        xPlot3 = xPlot3(xPlot3>=zoneDiv(i) & xPlot3<zoneDiv(i+1));
        FRI_thresh = NaN*ones(200,3);
        FRI_thresh(1:length(diff(xPlot1)),1) = diff(xPlot1)';
        FRI_thresh(1:length(diff(xPlot2)),2) = diff(xPlot2)';
        FRI_thresh(1:length(diff(xPlot3)),3) = diff(xPlot3)';
        mFRI_thresh = NaN(3,1);  mFRI_ci = NaN(3,2);
        for j = 1:3
            vals = FRI_thresh(FRI_thresh(:,j)>0, j);
            mFRI_thresh(j) = mean(vals);
            if length(vals) > 1
                bm = NaN(1000,1);
                for b = 1:1000
                    bm(b) = mean(vals(randi(length(vals),length(vals),1)));
                end
                mFRI_ci(j,:) = prctile(bm, [2.5 97.5]);
            end
        end
        xp = [plot_in(i)-0.25, plot_in(i), plot_in(i)+0.25];
        errorbar(xp, mFRI_thresh, mFRI_thresh-mFRI_ci(:,1), ...
            mFRI_ci(:,2)-mFRI_thresh, '.k', 'color', [.5 .5 .5])
        hold on
        in2 = find(PeakAnalysis.threshValues == PeakAnalysis.threshValues(end), 1);
        errorbar(xp(in2), mFRI_thresh(in2), ...
            mFRI_thresh(in2)-mFRI_ci(in2,1), mFRI_ci(in2,2)-mFRI_thresh(in2), '+k')
    end
    set(gca, 'TickDir', 'out', 'FontSize', FS, ...
        'xtick', 1:length(zoneDiv)-1, 'xticklabel', fliplr(1:length(zoneDiv)-1), ...
        'position', [0.17+wm*range(zoneDiv) .55 .10 .175])
    box off
    xlabel('zone')
    ylabel({'zone-specific mean FRI', '+/- 95% ci (years fire^-^1)'})
    title({'(b) Sensitivity to', 'alternative thresholds'}, 'fontweight', 'bold')
end

% -- Panel (c): Local SNI time series -------------------------------------
subplot(3,5,11:13)
if PeakAnalysis.threshType == 1
    SNI_plot = CharThresh.SNI .* ones(size(Charcoal.peak));
else
    SNI_plot = CharThresh.SNI;
end
plot(Charcoal.ybpI, SNI_plot, 'k')
y_lim = [0 10];
xlim([zoneDiv(1), zoneDiv(end)]); ylim(y_lim)
if y_lim(2) < 20;    y_tick = 0:2:y_lim(2);
elseif y_lim(2) < 50; y_tick = 0:5:y_lim(2);
else;                 y_tick = 0:10:y_lim(2);
end
set(gca, 'FontSize', FS, 'XDir', 'reverse', 'XMinorTick', 'on', 'YMinorTick', 'on', ...
    'Position', [0.1 .35 wm*range(zoneDiv) .12], 'TickDir', 'out', ...
    'XTick', 0:1000:max(zoneDiv), 'XTickLabel', 0:max(zoneDiv)/1000, 'YTick', y_tick)
hold on; plot([zoneDiv(1) zoneDiv(end)], [3 3], 'k--'); box off
xlabel('time (cal. yr BP x 1000)', 'FontSize', FS)
ylabel('signal-to-noise index', 'FontSize', FS)
title('(c) Local signal-to-noise index', 'FontSize', FS, 'FontWeight', FW)

% -- Panel (d): Global SNI boxplot ----------------------------------------
subplot(3,5,14)
boxplot(CharThresh.SNI, 'colors', 'kkk', 'symbol', '.k')
set(gca, 'TickDir', 'out', 'FontSize', FS, 'YMinorTick', 'off', ...
    'position', [0.17+wm*range(zoneDiv) .35 .10 .12], ...
    'xtick', 1, 'xticklabel', '', 'ylim', y_lim, 'YTick', y_tick, 'yscale', 'linear')
box off; grid off; hold on
plot([zoneDiv(1) zoneDiv(end)], [3 3], 'k--')
title('(d) Global signal-to-noise index', 'FontSize', FS, 'FontWeight', FW)
xlabel({'global signal-to-', 'noise distribution'}, 'FontSize', FS)
ylabel('')
text(1.2, median(CharThresh.SNI), num2str(nanmedian(CharThresh.SNI)), ...
    'FontSize', FS, 'BackgroundColor', 'w')

end
