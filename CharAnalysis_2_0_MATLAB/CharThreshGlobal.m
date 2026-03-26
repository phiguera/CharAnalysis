function [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment, ...
        PeakAnalysis, site, Results, plotData, bkgSensIn)
% CharThreshGlobal    Calculate global threshold value for charcoal peak ID.
%   [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment,
%       PeakAnalysis, site, Results, plotData, bkgSensIn)
%
%   Determines a threshold value based on the distribution of Cpeak across
%   the entire record.  The noise component of Cpeak is modelled as either
%   a 0-mean Gaussian (cPeak = 1, residuals) or a 1-mean Gaussian
%   (cPeak = 2, ratios), or by a Gaussian mixture model.
%
%   INPUTS
%     Charcoal     : struct containing peak (Cpeak series)
%     Pretreatment : struct  (zoneDiv)
%     PeakAnalysis : struct  (cPeak, threshMethod, threshValues)
%     site         : site name string (for plot titles)
%     Results      : struct  (allFigures)
%     plotData     : scalar flag  1 = draw Figure 2
%                    (explicit argument in v2.0; was global in v1.1)
%     bkgSensIn    : scalar flag  0 = normal run,  1 = called from
%                    bkgCharSensitivity (reuse Figure 2 axis limits for
%                    posThreshBins instead of recomputing from data)
%                    (explicit argument in v2.0; was global in v1.1)
%
%   OUTPUT
%     CharThresh : struct with fields:
%                    possible  : vector of candidate threshold values
%                    pos       : [N x 4] positive threshold matrix
%                    neg       : [N x 1] (or [N x 4]) negative threshold
%                    noisePDF  : noise probability density function
%                    SNI       : signal-to-noise index (scalar)
%                    GOF       : goodness-of-fit placeholder (-999 vector)
%
%   v2.0 changes vs v1.1
%     - plotData   passed as explicit argument; global removed.
%     - bkgSensIn  passed as explicit argument; global removed.
%       In v1.1 the global declaration was commented out but the variable
%       was still being written by bkgCharSensitivity and read here via
%       side effects.  It is now a clean argument in both directions.
%     - hist() replaced by charHistCounts() in the plot block.
%     - Bug fix (silent): in the data-defined threshold bin-lookup loop
%       (threshMethod > 1) the second operand of the closest-bin comparison
%       was PeakAnalysis.threshValues(i) instead of thresh(i), meaning the
%       tiebreak compared a percentile value against a CHAR value.  Both
%       operands now consistently use thresh(i).
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── LOCAL VARIABLES ───────────────────────────────────────────────────────
zoneDiv     = Pretreatment.zoneDiv;
figPosition = [0.1178  0.2158  0.8648  0.6943];

%% ── DEFINE CANDIDATE THRESHOLD BINS ─────────────────────────────────────
%
%   On a normal run (bkgSensIn == 0) the bins span the full range of Cpeak.
%   During the background sensitivity loop (bkgSensIn == 1) the bins are
%   fixed to the x-axis limits already shown on Figure 2, so that the
%   contour plot in bkgCharSensitivity uses a consistent x-axis across all
%   smoothing windows.

if bkgSensIn == 0
    posThreshBins = min(Charcoal.peak) : ...
                    range(Charcoal.peak) / 250 : ...
                    max(Charcoal.peak);
else
    figure(2);
    xlimIn        = get(gca, 'xlim');
    posThreshBins = xlimIn(1) : range(xlimIn) / 250 : xlimIn(2);
end

CharThresh.possible = posThreshBins;

%% ── USER-DEFINED THRESHOLD (threshMethod == 1) ───────────────────────────

if PeakAnalysis.threshMethod == 1

    CharThresh.pos = ones(length(Charcoal.peak), ...
                          length(PeakAnalysis.threshValues));

    for i = 1:length(PeakAnalysis.threshValues)
        tv   = PeakAnalysis.threshValues(i);
        in1  = find(CharThresh.possible >= tv, 1,        'first');
        in2  = find(CharThresh.possible <= tv, 1,        'last');

        % Closest bin to the requested value
        if abs(CharThresh.possible(in1) - tv) <= ...
           abs(CharThresh.possible(in2) - tv)
            inFinal = in1;
        else
            inFinal = in2;
        end

        CharThresh.pos(:,i) = CharThresh.possible(inFinal);
    end

    CharThresh.neg      = -99 * ones(length(Charcoal.peak), ...
                                     length(PeakAnalysis.threshValues));
    CharThresh.noisePDF = -99;

end

%% ── DATA-DEFINED THRESHOLD (threshMethod == 2 or 3) ─────────────────────

if PeakAnalysis.threshMethod > 1

    %% Estimate noise distribution parameters (muHat, sigmaHat)

    if PeakAnalysis.threshMethod == 2

        if PeakAnalysis.cPeak == 1
            % Residuals: noise assumed to be a zero-mean Gaussian.
            % Estimate sigma from the negative half of Cpeak (mirrored).
            negVals  = Charcoal.peak(Charcoal.peak <= 0);
            sigmaHat = std([negVals; abs(negVals)]);
            muHat    = 0;
        else
            % Ratios: noise assumed to be a one-mean Gaussian.
            % Shift values so the distribution is centred on 0, mirror,
            % then shift back.
            subVals  = Charcoal.peak(Charcoal.peak <= 1) - 1;
            sigmaHat = std([subVals; abs(subVals)] + 1);
            muHat    = 1;
        end

        CharThresh.noisePDF = normpdf(posThreshBins, muHat, sigmaHat);

    end

    if PeakAnalysis.threshMethod == 3

        % Gaussian mixture model (two-component)
        [mu, sig] = GaussianMixture(Charcoal.peak, 2, 2, false);

        if mu(1) == mu(2)
            beep
            disp('WARNING: poor fit of Gaussian mixture model.')
            disp('         Re-fitting starting with three classes.')
            [mu, sig] = GaussianMixture(Charcoal.peak, 3, 2, false);
        end

        % Noise component = the Gaussian with the smaller mean
        noiseIdx = find(mu == min(mu), 1);
        sigmaHat = sig(noiseIdx);
        muHat    = mu(noiseIdx);

        CharThresh.noisePDF = normpdf(posThreshBins, muHat, sigmaHat);

    end

    %% Map percentile thresholds to the nearest bin in posThreshBins

    thresh = norminv(PeakAnalysis.threshValues, muHat, sigmaHat);

    CharThresh.pos = ones(length(Charcoal.peak), ...
                          length(PeakAnalysis.threshValues));

    for i = 1:length(PeakAnalysis.threshValues)
        in1 = find(CharThresh.possible >= thresh(i), 1,        'first');
        in2 = find(CharThresh.possible <= thresh(i), 1,        'last');

        % v2.0 bug fix: both operands now use thresh(i).
        % v1.1 used PeakAnalysis.threshValues(i) on the right-hand side,
        % which compared a percentile (e.g. 0.95) against a CHAR value.
        if abs(CharThresh.possible(in1) - thresh(i)) <= ...
           abs(CharThresh.possible(in2) - thresh(i))
            inFinal = in1;
        else
            inFinal = in2;
        end

        CharThresh.pos(:,i) = CharThresh.possible(inFinal);
    end

    % Negative threshold: mirror of the lowest positive percentile
    threshNeg      = norminv(1 - PeakAnalysis.threshValues(3), ...
                             muHat, sigmaHat);
    CharThresh.neg = threshNeg * ones(length(Charcoal.peak), 1);

end

%% ── SIGNAL-TO-NOISE INDEX (Kelly et al.) ────────────────────────────────

signal = Charcoal.peak(Charcoal.peak >  CharThresh.pos(1,4));
noise  = Charcoal.peak(Charcoal.peak <= CharThresh.pos(1,4));

if ~isempty(signal) && length(noise) > 2
    CharThresh.SNI = (1 / length(signal)) .* ...
        sum((signal - mean(noise)) ./ std(noise)) .* ...
        ((length(noise) - 2) / length(noise));
else
    CharThresh.SNI = 0;
end

%% ── GOODNESS-OF-FIT PLACEHOLDER ─────────────────────────────────────────
%   GOF is computed sample-by-sample in CharThreshLocal; for the global
%   threshold a single KS test is not applied, so this field is filled with
%   a sentinel value for downstream compatibility.

CharThresh.GOF = -999 * ones(size(Charcoal.peak));

%% ── PLOT FIGURE 2 ────────────────────────────────────────────────────────

if plotData == 1 && Results.allFigures == 1

    figure(2); clf;
    set(gcf, 'color', 'white', ...
             'name',  ['Peak CHAR distribution, threshold values, ' ...
                       'and noise distribution (if selected)'], ...
             'units', 'normalized', ...
             'position', figPosition)

    % Histogram of Cpeak using charHistCounts (replaces hist())
    [n, x] = charHistCounts(Charcoal.peak, CharThresh.possible);
    h = bar(x, n/sum(n), 1);
    set(h, 'faceColor', [.5 .5 .5], 'edgeColor', [.5 .5 .5])
    hold on

    % Noise PDF overlay (only when data-defined threshold is used)
    if PeakAnalysis.threshMethod > 1
        plot(CharThresh.possible, ...
             CharThresh.noisePDF * mean(diff(CharThresh.possible)), ...
             'k', 'lineWidth', 2)
    end

    % All candidate threshold levels (dashed), selected level (solid red)
    for t = 1:4
        plot([CharThresh.pos(1,t)  CharThresh.pos(1,t)], ...
             [0  max(n/sum(n))], '--k')
    end
    plot([CharThresh.pos(1,4)  CharThresh.pos(1,4)], ...
         [0  max(n/sum(n))], 'r')

    set(gca, 'tickdir', 'out', 'box', 'off')
    ylabel('proportion or scaled density')
    xlabel('peak CHAR (# cm^-^2 yr^-^1)')
    title([char(site) ':  ' num2str(zoneDiv(1)) ' to ' ...
           num2str(zoneDiv(end)) ' cal. yr BP'])

    % Annotation: threshold value and SNI
    xAnnot = CharThresh.pos(1,4);
    yAnnot = max(n/sum(n));

    if PeakAnalysis.threshMethod == 1
        legend('peak CHAR distribution', ...
               'possible threshold values', ...
               'selected threshold value')
        text(1.5*xAnnot, 0.75*yAnnot, ...
             ['Threshold = ' num2str(xAnnot)], ...
             'backgroundColor', 'w')
        text(1.75*xAnnot, 0.60*yAnnot, ...
             ['signal-to-noise index = ' num2str(CharThresh.SNI)], ...
             'backgroundColor', 'w')

    elseif PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak == 1
        legend('peak CHAR distribution', ...
               'estimated noise PDF, via 0-mean Gaussian', ...
               'possible threshold values', ...
               'selected threshold value')
        text(1.5*xAnnot, 0.75*yAnnot, ...
             ['Threshold = ' num2str(xAnnot) ';  ' ...
              num2str(round(PeakAnalysis.threshValues(4)*100)) ...
              '^t^h percentile'], ...
             'backgroundColor', 'w')
        text(1.75*xAnnot, 0.60*yAnnot, ...
             ['signal-to-noise index = ' num2str(CharThresh.SNI)], ...
             'backgroundColor', 'w')

    elseif PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak == 2
        legend('peak CHAR distribution', ...
               'estimated noise PDF, via 1-mean Gaussian', ...
               'possible threshold values', ...
               'selected threshold value')
        text(1.5*xAnnot, 0.75*yAnnot, ...
             ['Threshold = ' num2str(xAnnot) ';  ' ...
              num2str(round(PeakAnalysis.threshValues(4)*100)) ...
              '^t^h percentile'], ...
             'backgroundColor', 'w')
        text(1.75*xAnnot, 0.60*yAnnot, ...
             ['signal-to-noise index = ' num2str(CharThresh.SNI)], ...
             'backgroundColor', 'w')

    elseif PeakAnalysis.threshMethod == 3
        legend('peak CHAR distribution', ...
               'estimated noise PDF, via Gaussian mixture model', ...
               'possible threshold values', ...
               'selected threshold value')
        text(1.5*xAnnot, 0.75*yAnnot, ...
             ['Threshold = ' num2str(xAnnot) ';  ' ...
              num2str(round(PeakAnalysis.threshValues(4)*100)) ...
              '^t^h percentile'], ...
             'backgroundColor', 'w')
        text(1.75*xAnnot, 0.60*yAnnot, ...
             ['signal-to-noise index = ' num2str(CharThresh.SNI)], ...
             'backgroundColor', 'w')
    end

end

end
