function [FRIyr, FRI, smFRIyr, smFRI, smFRIci] = smoothFRI(yr, peaks, ...
    winWidth, alpha, nBoot, FRI_param, makeFigure)
% smoothFRI   Compute and optionally plot smoothed fire-return intervals.
%
%   [FRIyr, FRI, smFRIyr, smFRI, smFRIci] = smoothFRI(yr, peaks, ...
%       winWidth, alpha, nBoot, FRI_param, makeFigure)
%
%   yr         - year for each sample (column vector)
%   peaks      - binary series, same size as yr; 1 identifies a peak
%   winWidth   - [yr] window width for smoothing
%   alpha      - confidence level for CIs (e.g. 0.05 for 95 %)
%   nBoot      - number of bootstrap samples for CIs
%   FRI_param  - 0 = std-dev bands; 1 = bootstrapped (1-alpha)*100 % CIs
%   makeFigure - 1 = produce figure; 0 = suppress figure
%
%   PEH, July 2010.
%   v2.0: smooth() -> charLowess(); zeros*ones bug fixed to NaN*ones;
%         kstest2 three-output form corrected to two-output form.
%         AnDarksamtest retained (bundled).

%% ── Input defaults ───────────────────────────────────────────────────────
if nargin < 7
    makeFigure = 0;
end

stepLength = 100;   % [yr] Time step to move window across record.
maxFRIs    = 100;   % Conservative upper bound on FRIs per window.

%% ── Derived parameters ───────────────────────────────────────────────────
steps   = -70:stepLength:max(max(yr));          % [cal yr BP] time steps
nSteps  = length(steps);
peakYr  = yr(peaks == 1);                       % years with peaks
spaceIn = round(winWidth/stepLength - 1);
x_lim   = [-60, ceil(max(peakYr)/1000)*1000];  % x-axis limits for figure

%% ── Allocate output arrays ───────────────────────────────────────────────
% v1 had "zeros*ones" (evaluates to 0, not NaN). Fixed to NaN*ones.
FRImatrix.FRI    = NaN * ones(maxFRIs, nSteps);
FRImatrix.params = NaN * ones(3, nSteps);
FFmatrix         = NaN * ones(nSteps, 1);
p                = NaN * ones(nSteps, 1);   % Anderson-Darling p-values
p2               = NaN * ones(nSteps, 1);   % KS two-sample p-values
p3               = NaN * ones(nSteps, 1);   % Wilcoxon rank-sum p-values

%% ── Main loop: FRI stats and between-window comparisons ─────────────────
for i = 1:nSteps
    startYr  = steps(i) - 0.5*winWidth;
    endYr    = steps(i) + 0.5*winWidth;
    lengthIn = yr >= startYr & yr < endYr;
    yrIn     = find(peakYr >= startYr & peakYr < endYr);
    fris     = diff(peakYr(yrIn));              % FRIs for this window

    if length(fris) > 1
        FRImatrix.FRI(1:length(fris), i) = fris;
        FFmatrix(i) = length(peakYr(yrIn)) / range(yr(lengthIn));

        if FRI_param == 0           % Standard-deviation bands
            FRImatrix.params(1, i) = mean(fris);
            FRImatrix.params(2, i) = mean(fris) + std(fris)/sqrt(length(fris));
            FRImatrix.params(3, i) = mean(fris) - std(fris)/sqrt(length(fris));
        else                        % Bootstrapped CIs
            FRImatrix.params(1, i) = mean(fris);
            boot_means = NaN(nBoot, 1);
            for b = 1:nBoot
                boot_means(b) = mean(fris(randi(length(fris), length(fris), 1)));
            end
            FRImatrix.params(2:3, i) = prctile(boot_means, ...
                [alpha/2*100, (1-alpha/2)*100]);
        end

        % Between-window distribution comparisons (non-overlapping windows)
        if i > spaceIn
            x1 = FRImatrix.FRI(:, i - spaceIn);
            x1 = x1(x1 > -999 & ~isnan(x1));
            if ~isempty(x1)
                x2  = fris;
                % Anderson-Darling test (bundled AnDarksamtest — retained)
                x12 = NaN * ones(length(x1) + length(x2), 2);
                x12(1:length(x1), 1) = x1;
                x12(1:length(x1), 2) = 1;
                x12(length(x1)+1:end, 1) = x2;
                x12(length(x1)+1:end, 2) = 2;
                p(i)  = AnDarksamtest(x12, alpha);
                % Two-sample KS test — two-output form (modern MATLAB)
                [~, p2(i)] = kstest2(x1, x2);
                % Wilcoxon rank-sum test
                p3(i) = ranksum(x1, x2);
            end
        end
    end
end

%% ── Smooth FRI parameters through time ───────────────────────────────────
in      = find(FRImatrix.params(1,:) > -999 & ~isnan(FRImatrix.params(1,:)));
x       = steps(in);
smoothIn = (winWidth/stepLength) / length(x);

% charLowess replaces smooth(...,'lowess')
y  = charLowess(FRImatrix.params(1, in)', smoothIn);   % smoothed mean FRI
y3 = charLowess(FRImatrix.params(2, in)', smoothIn);   % smoothed upper CI
y2 = charLowess(FRImatrix.params(3, in)', smoothIn);   % smoothed lower CI
y4 = FRImatrix.params(1, in)';                         % unsmoothed mean

X = [x';       flipud(x')];
Y = [y2;       flipud(y3)];

%% ── Collect output ───────────────────────────────────────────────────────
FRIyr   = peakYr(1:end-1);
FRI     = diff(peakYr);
smFRIyr = x;
smFRI   = y;
smFRIci = [y3, y2];

%% ── Optional figure ──────────────────────────────────────────────────────
if makeFigure == 1
    figure(10); clf; set(gcf, 'color', 'w')
    plot(peakYr(1:end-1), diff(peakYr), 's', ...
        'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', [0.6 0.6 0.6])
    hold on
    plot(x, y, '-k', 'linewidth', 2)
    h = fill(X, Y, [0.8 0.8 0.8]);
    set(h, 'edgecolor', [0.8 0.8 0.8])
    set(gca, 'xdir', 'rev', 'yscale', 'log', 'xlim', x_lim, 'ylim', [50 1000])
    if FRI_param == 0
        ylabel([num2str(winWidth) '-yr mean FRI'])
        legend('FRI', 'mean FRI +/- 1 stdev.')
    else
        ylabel([num2str(winWidth) '-yr mean FRI'])
        legend('FRI', 'mean FRI', [num2str(100*(1-alpha)) '% CI'])
    end
    % Redraw data on top of fill
    plot(peakYr(1:end-1), diff(peakYr), 's', ...
        'markerfacecolor', [0.6 0.6 0.6], 'markeredgecolor', [0.6 0.6 0.6])
    plot(x, y,  '-k', 'linewidth', 2)
    plot(x, y4, '-k', 'linewidth', 1)
    grid on
end