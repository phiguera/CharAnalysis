function [z, GOF_i, SNI_i] = bkgCharSensitivity(Charcoal, CharThresh, ...
    PeakAnalysis, Pretreatment, Smoothing, Results, site)
% bkgCharSensitivity    Run CharAnalysis with multiple background windows.
%   [z, GOF_i, SNI_i] = bkgCharSensitivity(Charcoal, CharThresh, ...
%       PeakAnalysis, Pretreatment, Smoothing, Results, site)
%
%   Analyses a charcoal record using multiple smoothing window widths (all
%   with the same smoothing method) and plots graphs illustrating the
%   sensitivity of results to the choice of background smoothing window.
%
%   v2.0 changes vs v1.1
%     - global plotData  removed; set locally and passed to CharThreshGlobal.
%     - global bkgSensIn removed; set locally and passed to CharThreshGlobal.
%     - hist()    replaced by histcounts().
%     - plotyy()  replaced by yyaxis in subplot (2,2,4).
%     - print string form replaced by function form for consistency.

%% ── Local flags (replaced globals) ──────────────────────────────────────
plotData = 0;       % Do not redraw Figure 2 during the sensitivity loop.

% bkgSensIn tells CharThreshGlobal to reuse Figure 2 axis limits for the
% threshold bin vector instead of recomputing from each smoothed series.
if PeakAnalysis.threshType == 1    % Global threshold only.
    bkgSensIn = 1;
else
    bkgSensIn = 0;                 % Not used for local threshold path.
end

%% ── Background window range ──────────────────────────────────────────────
% Maximum window must be shorter than the record length.
if max(Charcoal.ybpI) >= 1000
    bkgMax = 1000;
else
    bkgBin   = 100:100:900;
    [n, edges] = histcounts(max(Charcoal.ybpI), [bkgBin, bkgBin(end)+100]);
    bkgMax   = bkgBin(n > 0);
end

if PeakAnalysis.threshType == 1
    bkgMin = 100;
else
    % Smallest window such that each local threshold has >= 30 samples.
    bkgMin   = 30 * Pretreatment.yrInterp;
    bkgBin   = 100:100:500;
    [n, edges] = histcounts(bkgMin, [bkgBin, bkgBin(end)+100]);
    bkgMin   = bkgBin(n > 0);
end

%% ── Allocate output arrays ───────────────────────────────────────────────
bkSmooth = bkgMin:100:bkgMax;      % [yr] Windows to evaluate.

if PeakAnalysis.threshType == 1    % Global threshold
    figure(2); xlimIn = get(gca, 'xlim');
    x = xlimIn(1) : range(xlimIn)/250 : xlimIn(2);
    x = x(x > 0);
    z     = NaN * ones(length(bkSmooth), length(x));
    SNI_i = NaN * ones(length(bkSmooth), 1);
    GOF_i = NaN * ones(length(bkSmooth), 1);
else                               % Local threshold
    x     = 1:3;
    z     = NaN * ones(length(bkSmooth), 4);
    mFRI  = NaN * ones(length(bkSmooth), 4, 2);
    SNI_i = NaN * ones(length(bkSmooth), length(CharThresh.SNI));
    GOF_i = NaN * ones(length(bkSmooth), length(CharThresh.GOF));
end

SmoothingI = Smoothing;

%% ── Sensitivity loop ─────────────────────────────────────────────────────
for i = 1:length(bkSmooth)
    SmoothingI.yr = bkSmooth(i);
    disp(['    C_background sensitivity iteration ', num2str(i), ...
          ' of ', num2str(length(bkSmooth)), ...
          ': window = ', num2str(bkSmooth(i)), ' yrs.'])

    % Smooth charcoal record with this window width.
    [Charcoal] = CharSmooth(Charcoal, Pretreatment, SmoothingI, Results, plotData);

    % Calculate peak CHAR component.
    if PeakAnalysis.cPeak == 1
        Charcoal.peak = Charcoal.accI - Charcoal.accIS;
    else
        Charcoal.peak = Charcoal.accI ./ Charcoal.accIS;
    end

    % Define thresholds — pass plotData and bkgSensIn as explicit args.
    if PeakAnalysis.threshType == 1
        [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment, ...
            PeakAnalysis, site, Results, plotData, bkgSensIn);
    else
        [CharThresh] = CharThreshLocal(Charcoal, SmoothingI, ...
            PeakAnalysis, site, Results, plotData);
    end

    % Identify peaks.
    [Charcoal, CharThresh] = CharPeakID(Charcoal, Pretreatment, ...
        PeakAnalysis, CharThresh);

    % Store results.
    if PeakAnalysis.threshType == 1
        z(i,:)    = sum(Charcoal.charPeaks);
        GOF_i(i)  = CharThresh.GOF(1);
        SNI_i(i)  = CharThresh.SNI;
    else
        z(i,:)      = sum(Charcoal.charPeaks);
        mFRI(i,:,1) = mean(diff( ...
            Charcoal.ybpI(Charcoal.charPeaks(:,end) > 0)));
        mFRI(i,:,2) = std(diff( ...
            Charcoal.ybpI(Charcoal.charPeaks(:,end) > 0)));
        GOF_i(i,:)  = CharThresh.GOF;
        SNI_i(i,:)  = CharThresh.SNI;
    end

    clear CharThresh
end

%% ── Plot results (Figure 10) ─────────────────────────────────────────────
figure(10); clf
set(gcf, 'color', 'w', ...
    'name',     'Sensitivity to different background windows', ...
    'units',    'normalized', ...
    'position', [0.0435  0.0571  0.8649  0.6943])

if PeakAnalysis.threshType == 1   % ── Global threshold: contour plot ─────

    [X, Y] = meshgrid(x, bkSmooth);
    contourf(X, Y, z)
    colormap(1 - gray); brighten(0.4); colorbar
    xlim([x(1), prctile(x, 50)])
    set(gca, 'TickDir', 'out')
    xlabel('threshold (C_p_e_a_k units)')
    ylabel('smoothing window width (yr)')
    text(1.2*max(x), min(bkSmooth) + range(bkSmooth)/2, ...
        'peaks identified (#)', 'Rotation', 270, 'HorizontalAlignment', 'Center')
    title('Peaks identified with varying threshold and C_b_a_c_k_g_r_o_u_n_d criteria')
    box off

else                               % ── Local threshold: diagnostic panels ─

    x_lim = [min(bkSmooth) - 0.5*mean(diff(bkSmooth)), ...
              max(bkSmooth) + 0.5*mean(diff(bkSmooth))];

    subplot(2,2,1)
    boxplot(GOF_i')
    set(gca, 'XTick', 1:length(bkSmooth), 'XTickLabel', bkSmooth, ...
        'TickDir', 'out')
    ylim([0 1])
    ylabel('KS-test results (p-value)')
    title('(a) Noise distribution goodness of fit')
    box off

    subplot(2,2,2)
    boxplot(SNI_i')
    set(gca, 'XTick', 1:length(bkSmooth), 'XTickLabel', bkSmooth, ...
        'yscale', 'linear', 'TickDir', 'out')
    ylim([0, min([max(max(SNI_i)), 10])])
    hold on
    plot(ones(length(bkSmooth), 1) * 3, 'k--')
    ylabel('signal-to-noise index')
    title('(b) Signal-to-noise index')
    box off

    subplot(2,2,3)
    plot(bkSmooth, nanmedian(GOF_i, 2) + nanmedian(SNI_i, 2), '-ok', ...
        'linewidth', 2)
    xlim(x_lim)
    xlabel('smoothing window width (yr)')
    ylabel('sum of median SNI & GOF')
    set(gca, 'TickDir', 'out')
    title('(c) Sum of median SNI and GOF value')
    box off

    subplot(2,2,4)
    % ── yyaxis replaces plotyy ────────────────────────────────────────────
    yyaxis left
    errorbar(bkSmooth, mFRI(:,end,1), mFRI(:,end,2), 'ko-', 'linewidth', 2)
    ylabel('fire-return interval (yr)', 'color', 'k')
    set(gca, 'ycolor', 'k', 'xlim', x_lim, 'TickDir', 'out', ...
        'yscale', 'linear')

    yyaxis right
    plot(bkSmooth, z(:,end), 'bo-', 'linewidth', 2)
    ylabel('# of peaks identified', 'color', 'b', ...
        'rotation', 270, 'verticalalignment', 'bottom')
    set(gca, 'ycolor', 'b', 'xlim', x_lim)

    xlabel('smoothing window width (yr)')
    title('(d) Mean fire-return interval +/- standard deviation')
    box off

end

%% ── Save figure ──────────────────────────────────────────────────────────
if Results.saveFigures == 1
    set(gcf, 'PaperPositionMode', 'auto', 'PaperType', 'uslegal')
    orient(gcf, 'landscape')
    print('-dpdf', '-r300', '10_sensitivity_to_C_background.pdf')
    print('-dtiff', '-r300', '10_sensitivity_to_C_background.tif')
end