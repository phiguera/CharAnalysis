function [Charcoal] = CharSmooth(Charcoal, Pretreatment, Smoothing, ...
                                  Results, plotData)
% CharSmooth    Smooth charcoal record to estimate low-frequency Cbackground.
%   [Charcoal] = CharSmooth(Charcoal, Pretreatment, Smoothing,
%                            Results, plotData)
%
%   Smooths the interpolated charcoal series (Charcoal.accI) using one of
%   five methods to estimate background CHAR (Cbackground).  The selected
%   result is stored in Charcoal.accIS.
%
%   SMOOTHING METHODS  (Smoothing.method)
%     1  Lowess         locally-weighted linear regression
%     2  Robust Lowess  same, with iterative outlier down-weighting
%     3  Moving average simple arithmetic mean within window
%     4  Running median  median within window, then robust Lowess pass
%     5  Running mode    modal bin within window, then robust Lowess pass
%
%   INPUTS
%     Charcoal     : struct containing accI (resampled CHAR) and ybpI, ybp
%     Pretreatment : struct  (yrInterp, zoneDiv, transform)
%     Smoothing    : struct  (method [1-5], yr [window width in years])
%     Results      : struct  (allFigures)
%     plotData     : scalar flag  1 = draw Figure 1 subplot 2
%                    (explicit argument in v2.0; was global in v1.1)
%
%   OUTPUT
%     Charcoal.accIS : smoothed CHAR series chosen by Smoothing.method
%
%   v2.0 changes vs v1.1
%     - plotData passed as explicit argument; global removed.
%     - All smooth() calls (Curve Fitting Toolbox) replaced by charLowess()
%       which uses only base MATLAB.
%     - Secondary lowess pass in methods 4 and 5 remains plain 'lowess',
%       matching v1.1 behaviour.
%     - hist() call in plot block replaced by charHistCounts().
%
%   NOTE ON METHOD CHOICES
%     Methods 1 and 2 differ only in robustness: method 1 (lowess) weights
%     by distance alone; method 2 (rlowess) additionally down-weights
%     residual outliers.  Changing method 1 to rlowess would make it
%     identical to method 2 and remove a meaningful user choice, so method 1
%     intentionally remains plain lowess.  The secondary lowess pass in
%     methods 4 and 5 likewise uses plain lowess, matching v1.1 behaviour.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── LOCAL VARIABLES ───────────────────────────────────────────────────────
r = Pretreatment.yrInterp;
s = Smoothing.yr / r;           % window width in data-point units

N = length(Charcoal.accI);
charAccIS = NaN(N, 5);          % preallocate all five columns

%% ── METHOD 1: LOWESS ─────────────────────────────────────────────────────
%   Locally-weighted scatter-plot smooth using linear polynomial fitting.
%   Replaces: smooth(Charcoal.accI, s, 'lowess')

charAccIS(:,1) = charLowess(Charcoal.accI, s, 'lowess');

%% ── METHOD 2: ROBUST LOWESS ──────────────────────────────────────────────
%   Same as method 1 but with iterative bisquare re-weighting to suppress
%   the influence of outliers on the fitted curve.
%   Replaces: smooth(Charcoal.accI, s, 'rlowess')

charAccIS(:,2) = charLowess(Charcoal.accI, s, 'rlowess');

%% ── METHOD 3: MOVING AVERAGE ─────────────────────────────────────────────
%   Each output sample is the arithmetic mean of accI values within the
%   window.  At the record boundaries the window shrinks symmetrically.
%   Replaces: smooth(Charcoal.accI, s, 'moving')

charAccIS(:,3) = charLowess(Charcoal.accI, s, 'moving');

%% ── METHOD 4: RUNNING MEDIAN + ROBUST LOWESS ─────────────────────────────
%   Step 1: assign each sample the median accI value within the window.
%   Step 2: smooth the resulting series with robust Lowess.
%
%   The median suppresses large individual spikes within each window.
%   The lowess secondary pass then smooths across window boundaries.
%
%   Replaces: smooth(Charcoal.accI, s, 'moving') + smooth(..., s, 'lowess')

hw4 = round(s / 2);
for i = 1:N
    if i <= hw4
        % Near the start: expand window forward to maintain ~s points
        win = Charcoal.accI(1 : round(s));
    elseif i >= N - hw4
        % Near the end: expand window backward
        win = Charcoal.accI(N - hw4 : end);
    else
        win = Charcoal.accI(round(i - 0.5*s) : round(i + 0.5*s));
    end
    charAccIS(i, 4) = median(win);
end

charAccIS(:,4) = charLowess(charAccIS(:,4), s, 'lowess');

%% ── METHOD 5: RUNNING MODE + ROBUST LOWESS ───────────────────────────────
%   Step 1: divide accI values within the window into 100 equally-spaced
%           bins; assign each sample the centre of the most-populated bin.
%           If multiple bins share the maximum count, their centres are
%           averaged (median of modal centres).
%   Step 2: smooth the resulting series with robust Lowess.
%
%   The modal bin naturally ignores large outliers that fall in sparsely
%   populated bins.  The lowess secondary pass smooths across windows.
%
%   Replaces: hist() + smooth(..., s, 'lowess')
%   hist() replaced by charHistCounts().

hw5  = round(s / 2);
nBin = 100;

for i = 1:N
    if i <= hw5
        win = Charcoal.accI(1 : round(s));
    elseif i >= N - hw5
        win = Charcoal.accI(N - hw5 : end);
    else
        win = Charcoal.accI(round(i - 0.5*s) : round(i + 0.5*s));
    end

    [n_mode, x_mode] = charHistCounts(win, nBin);
    modal_centres     = x_mode(n_mode == max(n_mode));
    charAccIS(i, 5)   = median(modal_centres);   % median of ties
end

charAccIS(:,5) = charLowess(charAccIS(:,5), s, 'lowess');

%% ── PLOT (FIGURE 1, SUBPLOT 2) ───────────────────────────────────────────
%   Only drawn when plotData == 1 and Results.allFigures == 1.
%   Subplot 1 is drawn by CharPretreatment; both live in Figure 1.

if plotData == 1 && Results.allFigures == 1

    figure(1);
    subplot(2,1,2)

    h = bar(Charcoal.ybpI, Charcoal.accI, 1);
    set(h, 'facecolor', [.5 .5 .5], 'edgecolor', [.5 .5 .5])
    hold on
    plot(Charcoal.ybpI, charAccIS, 'linewidth', 1.5)
    ylim([0, max(Charcoal.accI)]);

    legend('C_i_n_t_e_r_p_o_l_a_t_e_d', ...
           'Lowess', ...
           'Robust Lowess', ...
           'Moving Average', ...
           'Moving Median', ...
           'Moving Mode')

    xlabel('time (cal. yr BP)')
    ylabel('CHAR (# cm^-^2 yr^-^1)')
    set(gca, 'xlim',    [min(Charcoal.ybp)-100  max(Charcoal.ybp)], ...
             'box',     'off', ...
             'tickdir', 'out', ...
             'xdir',    'rev')
    title(['(b) C_i_n_t_e_r_p_o_l_a_t_e_d and options for a ' ...
           num2str(Smoothing.yr) ' yr C_b_a_c_k_g_r_o_u_n_d'])

end

%% ── STORE SELECTED METHOD ────────────────────────────────────────────────
Charcoal.accIS = charAccIS(:, Smoothing.method);

end
