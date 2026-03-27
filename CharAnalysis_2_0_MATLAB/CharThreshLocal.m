function [CharThresh] = CharThreshLocal(Charcoal, Smoothing, ...
        PeakAnalysis, site, Results, plotData)
% CharThreshLocal    Calculate local (sliding-window) threshold for peak ID.
%   [CharThresh] = CharThreshLocal(Charcoal, Smoothing,
%       PeakAnalysis, site, Results, plotData)
%
%   Determines a per-sample threshold value based on the distribution of
%   Cpeak within a sliding window centred on each sample.  The noise
%   component within each window is modelled as a 0-mean Gaussian (cPeak=1,
%   residuals), a 1-mean Gaussian (cPeak=2, ratios), or a two-component
%   Gaussian mixture model (threshMethod=3).
%
%   INPUTS
%     Charcoal     : struct containing peak (Cpeak), ybpI, accI
%     Smoothing    : struct  (yr [window width in years])
%     PeakAnalysis : struct  (cPeak, threshMethod, threshValues, minCountP)
%     site         : site name string (for plot titles)
%     Results      : struct  (allFigures)
%     plotData     : scalar flag  1 = draw Figure 2 diagnostic plots
%                    (explicit argument in v2.0; was global in v1.1)
%
%   OUTPUT
%     CharThresh : struct with fields:
%                    pos  : [N x nThreshValues] positive threshold matrix
%                    neg  : [N x nThreshValues] negative threshold matrix
%                    SNI  : [N x 1] signal-to-noise index time series
%                    GOF  : [N x 1] KS goodness-of-fit p-values
%
%   v2.0 changes vs v1.1
%     - plotData passed as explicit argument; global removed.
%     - All smooth() calls (Curve Fitting Toolbox) replaced by charLowess().
%     - hist() replaced by charHistCounts() in the diagnostic plot block.
%     - Window-selection logic restructured using a single half-window
%       variable hw for clarity and symmetry.
%     - Bug fix: GMM poor-fit fallback now passes X (local window) to
%       GaussianMixture instead of Charcoal.peak (entire record), which
%       caused severe slowdown in v2.0 prior to this fix.
%     - Bug fix: subplot title index corrected from r/Smoothing.yr to
%       Smoothing.yr/r for dimensional consistency.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ?? LOCAL VARIABLES ???????????????????????????????????????????????????????
threshYr    = Smoothing.yr;
figPosition = [0.1178  0.2158  0.8648  0.6943];

r   = mean(diff(Charcoal.ybpI));   % record resolution (yr per sample)
hw  = round(0.5 * threshYr / r);   % half-window width in samples
P   = PeakAnalysis.threshValues(4);% percentile used for threshold + plots
nPk = length(Charcoal.peak);
nTV = length(PeakAnalysis.threshValues);
FS  = 8;

%% ?? PREALLOCATE OUTPUT ARRAYS ????????????????????????????????????????????
CharThresh.pos = NaN(nPk, nTV);
CharThresh.neg = NaN(nPk, nTV);
CharThresh.SNI = NaN(nPk, 1);
CharThresh.GOF = NaN(nPk, 1);

muHat    = NaN(nPk, 2);
sigmaHat = NaN(nPk, 2);
propN    = NaN(nPk, 2);

%% ?? OPEN FIGURE 2 IF PLOTTING ????????????????????????????????????????????
if plotData ~= 0 && Results.allFigures == 1
    plotIn = 1;
    figure(2); clf;
    set(gcf, 'color', 'w', ...
             'name',  'Local distributions of C_peak values', ...
             'units', 'normalized', ...
             'position', figPosition);
end

%% ?? PER-SAMPLE LOOP: FIT NOISE DISTRIBUTION AND COMPUTE THRESHOLD ????????
for i = 1:nPk

    % Progress report every 100 samples
    if mod(i, 100) == 0
        disp(['    Calculating local threshold ' num2str(i) ...
              ' of ' num2str(nPk)]);
    end

    %% Select window samples
    %
    %   Three boundary conditions expressed symmetrically via hw:
    %     Start  : not enough samples to the left  -> extend right
    %     End    : not enough samples to the right -> extend left
    %     Middle : full symmetric window available

    if i <= hw
        X = Charcoal.peak(1 : hw + i);
    elseif i > nPk - hw
        X = Charcoal.peak(i - hw : end);
    else
        X = Charcoal.peak(i - hw : i + hw);
    end

    % Remove NaN values before any distribution fitting.
    X = X(~isnan(X));

    % If fewer than 4 valid samples remain, fall back to simple Gaussian
    % using whatever valid data is available rather than skipping entirely.
    % Skipping leaves NaN thresholds which propagates through peak ID.
    if length(X) < 4
        if PeakAnalysis.cPeak == 1
            negVals = X(X <= 0);
            if isempty(negVals)
                sigmaHat(i,1) = std(X);
            else
                sigmaHat(i,1) = std([negVals; abs(negVals)]);
            end
            muHat(i,1) = 0;
        else
            subVals       = X(X <= 1) - 1;
            sigmaHat(i,1) = std([subVals; abs(subVals)] + 1);
            muHat(i,1)    = 1;
        end
        muHat(i,2)    = muHat(i,1);
        sigmaHat(i,2) = sigmaHat(i,1);
        propN(i,:)    = [1 0];
        % Compute threshold, SNI and GOF from fallback params then continue
        % to next sample — no further distribution fitting needed.
        CharThresh.pos(i,:) = norminv(PeakAnalysis.threshValues', ...
                                      muHat(i,1), sigmaHat(i,1));
        CharThresh.neg(i,:) = norminv(1 - PeakAnalysis.threshValues, ...
                                      muHat(i,1), sigmaHat(i,1));
        tPos = norminv(P, muHat(i,1), sigmaHat(i,1));
        sig_i   = X(X >  tPos);
        noise_i = X(X <= tPos);
        if ~isempty(sig_i) && length(noise_i) > 2
            CharThresh.SNI(i) = (1/length(sig_i)) .* ...
                sum((sig_i - mean(noise_i))./std(noise_i)) .* ...
                ((length(noise_i)-2)/length(noise_i));
        else
            CharThresh.SNI(i) = 0;
        end
        continue
    end

    %% Estimate local noise distribution

    if PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak == 1
        % Residuals: zero-mean Gaussian noise
        negVals     = X(X <= 0);
        sigmaHat(i) = std([negVals; abs(negVals)]);
        muHat(i)    = 0;

    elseif PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak ~= 1
        % Ratios: one-mean Gaussian noise
        subVals     = X(X <= 1) - 1;
        sigmaHat(i) = std([subVals; abs(subVals)] + 1);
        muHat(i)    = 1;

    elseif PeakAnalysis.threshMethod == 3
        % Gaussian mixture model
        if sum(X) == 0
            % Degenerate window — all values zero, GMM cannot fit.
            % Fall back to zero-mean (residuals) or one-mean (ratios)
            % Gaussian assumption for this window only.
            if PeakAnalysis.cPeak == 1
                muHat(i,1) = 0;
                negVals = X(X <= 0);
                if isempty(negVals)
                    sigmaHat(i,1) = std(X);
                else
                    sigmaHat(i,1) = std([negVals; abs(negVals)]);
                end
            else
                muHat(i,1)    = 1;
                subVals       = X(X <= 1) - 1;
                sigmaHat(i,1) = std([subVals; abs(subVals)] + 1);
            end
            muHat(i,2)    = muHat(i,1);
            sigmaHat(i,2) = sigmaHat(i,1);
            propN(i,:)    = [1 0];
        else
            [muHat(i,:), sigmaHat(i,:), ~, propN(i,:)] = ...
                GaussianMixture(X, 2, 2, false);

            if muHat(i,1) == muHat(i,2)
                beep
                disp('WARNING: poor fit of Gaussian mixture model;')
                disp('         re-fitting starting with three classes.')
                % Fix: pass X (local window), not Charcoal.peak (full record)
                [mu_tmp, sig_tmp] = GaussianMixture(X, 3, 2, false);
                muHat(i,1:2)    = mu_tmp(1:2);
                sigmaHat(i,1:2) = sig_tmp(1:2);
            end

            % Ensure noise component is first, signal second
            noiseIdx  = find(muHat(i,:) == min(muHat(i,:)), 1);
            signalIdx = find(muHat(i,:) == max(muHat(i,:)), 1);
            muHat(i,:)    = [muHat(i,noiseIdx)    muHat(i,signalIdx)];
            sigmaHat(i,:) = [sigmaHat(i,noiseIdx) sigmaHat(i,signalIdx)];
            propN(i,:)    = [propN(i,noiseIdx)    propN(i,signalIdx)];
        end
    end

    %% Compute local threshold values from inverse normal CDF
    CharThresh.pos(i,:) = norminv(PeakAnalysis.threshValues', ...
                                  muHat(i,1), sigmaHat(i,1));
    CharThresh.neg(i,:) = norminv(1 - PeakAnalysis.threshValues, ...
                                  muHat(i,1), sigmaHat(i,1));

    %% Signal-to-noise index (Kelly et al.)
    tPos    = norminv(P, muHat(i,1), sigmaHat(i,1));
    sig_i   = X(X >  tPos);
    noise_i = X(X <= tPos);

    if ~isempty(sig_i) && length(noise_i) > 2
        CharThresh.SNI(i) = (1 / length(sig_i)) .* ...
            sum((sig_i - mean(noise_i)) ./ std(noise_i)) .* ...
            ((length(noise_i) - 2) / length(noise_i));
    else
        CharThresh.SNI(i) = 0;
    end

    %% Goodness-of-fit: KS test of noise component vs fitted normal
    noiseForKS = X(X <= tPos);

    if length(noiseForKS) > 3
        ksBin = linspace(min(noiseForKS), max(noiseForKS), 101);
        ksCdf = normcdf(ksBin, muHat(i,1), sigmaHat(i,1));
        [~, ksP] = kstest(noiseForKS, [ksBin'  ksCdf']);
        CharThresh.GOF(i) = ksP;
    end

    %% Diagnostic plot: sample distributions (Figure 2 subplots)
    if plotData ~= 0 && Results.allFigures == 1

        plotStep  = max(round(length(Charcoal.peak) / 25), 1);
        inPlotSet = (hw * 2) : plotStep : nPk;

        if any(i == inPlotSet) && plotIn <= 25

            subplot(5, 5, plotIn)

            [n_h, x_h] = charHistCounts(X, 50);
            H1 = bar(x_h, n_h / sum(n_h), 1);
            set(H1, 'facecolor', [.75 .75 .75], ...
                    'edgecolor', [.75 .75 .75]);
            hold on

            if PeakAnalysis.threshMethod == 3
                pdf1 = normpdf(x_h, muHat(i,1), sigmaHat(i,1)) ...
                       * mean(diff(x_h)) * propN(i,1);
                pdf2 = normpdf(x_h, muHat(i,2), sigmaHat(i,2)) ...
                       * mean(diff(x_h)) * propN(i,2);
                plot(x_h, pdf1,        'k', 'linewidth', 1)
                plot(x_h, pdf2,        'k', 'linewidth', 1)
                plot(x_h, pdf1 + pdf2, 'b', 'linewidth', 1)
            else
                plot(x_h, normpdf(x_h, muHat(i,1), sigmaHat(i,1)) ...
                          * mean(diff(x_h)), 'k', 'linewidth', 2)
            end

            plot([tPos  tPos], [0  0.5], 'r')

            set(gca, 'tickdir', 'out', 'ylim', [0  0.25], 'fontsize', 10)
            box off

            xl = get(gca, 'xlim');
            yl = get(gca, 'ylim');

            titleYr = round(Charcoal.ybpI( ...
                round(mean(max(1, i - round(0.5*(Smoothing.yr/r))) : ...
                           min(nPk, i + round(0.5*(Smoothing.yr/r)))))));

            if plotIn == 1
                title([char(site) ':  ' num2str(titleYr) ' yr BP'], ...
                      'fontsize', FS)
            else
                title([num2str(titleYr) ' yr BP'], 'fontsize', FS)
            end

            text(xl(2), yl(2), ...
                 [{['SNI_i = '    num2str(round(CharThresh.SNI(i)*100)/100)]}, ...
                  {['KS p-val = ' num2str(round(CharThresh.GOF(i)*100)/100)]}, ...
                  {['t_i = '      num2str(round(tPos*1000)/1000)]}], ...
                 'horizontalalignment', 'right', ...
                 'verticalalignment',   'top', ...
                 'FontSize', FS);

            if plotIn == 11
                ylabel('proportion OR density (scaled)', 'FontSize', FS)
            end
            if plotIn == 23
                xlabel('CHAR (pieces cm^-^2 yr^-^1)', 'FontSize', FS)
            end

            plotIn = plotIn + 1;

        end
    end

end % per-sample loop

%% ?? SMOOTH THRESHOLD AND SNI SERIES WITH LOWESS ??????????????????????????
span = Smoothing.yr / r;

CharThresh.SNI(:,1) = charLowess(CharThresh.SNI(:,1), span, 'lowess');
CharThresh.SNI(CharThresh.SNI < 0, 1) = 0;

for i = 1:nTV
    CharThresh.pos(:,i) = charLowess(CharThresh.pos(:,i), span, 'lowess');
    CharThresh.neg(:,i) = charLowess(CharThresh.neg(:,i), span, 'lowess');
end

end