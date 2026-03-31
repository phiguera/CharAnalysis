function [Charcoal] = CharSmooth(Charcoal, Pretreatment, Smoothing, ...
                                  Results, plotData)
% CharSmooth    Smooth charcoal record to estimate low-frequency Cbackground.
%   [Charcoal] = CharSmooth(Charcoal, Pretreatment, Smoothing,
%                            Results, plotData)
%
%   Smooths the interpolated charcoal series (Charcoal.accI) using one of
%   five methods to estimate background CHAR (Cbackground). The selected
%   result is stored in Charcoal.accIS.
%
%   SMOOTHING METHODS  (Smoothing.method)
%     1  Lowess         locally-weighted linear regression
%     2  Robust Lowess  same, with iterative outlier down-weighting
%     3  Moving average simple arithmetic mean within window
%     4  Running median  median within window, then Lowess pass
%     5  Running mode    modal bin within window, then Lowess pass
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
%     - When the Curve Fitting Toolbox is available, smooth() is called
%       directly on Charcoal.accI (including any NaNs) for methods 1-3
%       and the Lowess passes in methods 4-5, guaranteeing results
%       identical to v1.1. charLowess() is used as the fallback when the
%       toolbox is absent; in that case NaN entries are bridged by linear
%       interpolation before smoothing and restored afterward.
%     - hist() call in plot block replaced by charHistCounts().
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% -- Local variables ------------------------------------------------------
r = Pretreatment.yrInterp;
s = Smoothing.yr / r;           % window width in data-point units

N = length(Charcoal.accI);
charAccIS = NaN(N, 5);          % preallocate all five columns

%% -- Curve Fitting Toolbox check ------------------------------------------
% When the toolbox is available, smooth() is called directly on the raw
% accI series (including any NaNs from record gaps), exactly as v1.1 did.
% When the toolbox is absent, charLowess() is used instead; it requires
% NaN-free input, so NaN entries are bridged and restored around the call.
v      = ver;
hasCFT = any(strcmp({v.Name}, 'Curve Fitting Toolbox'));

%% -- NaN bridging (charLowess path only) ----------------------------------
% Only needed when hasCFT is false. smooth() handles NaNs internally;
% charLowess() does not.
nanMask    = isnan(Charcoal.accI);
accI_clean = Charcoal.accI;
if ~hasCFT && any(nanMask)
    xAll = (1:N)';
    accI_clean(nanMask) = interp1(xAll(~nanMask), ...
                                   Charcoal.accI(~nanMask), ...
                                   xAll(nanMask), 'linear', 'extrap');
end

%% -- Method 1: Lowess -----------------------------------------------------
% Locally-weighted scatter-plot smooth using linear polynomial fitting.
% Replaces: smooth(Charcoal.accI, s, 'lowess')
if hasCFT
    charAccIS(:,1) = smooth(Charcoal.accI, s, 'lowess');
else
    charAccIS(:,1) = charLowess(accI_clean, s, 'lowess');
end

%% -- Method 2: Robust Lowess ----------------------------------------------
% Same as method 1 but with iterative bisquare re-weighting to suppress
% the influence of outliers on the fitted curve.
% Replaces: smooth(Charcoal.accI, s, 'rlowess')
if hasCFT
    charAccIS(:,2) = smooth(Charcoal.accI, s, 'rlowess');
else
    charAccIS(:,2) = charLowess(accI_clean, s, 'rlowess');
end

%% -- Method 3: Moving average ---------------------------------------------
% Each output sample is the arithmetic mean of accI values within the
% window. At the record boundaries the window shrinks symmetrically.
% Replaces: smooth(Charcoal.accI, s, 'moving')
if hasCFT
    charAccIS(:,3) = smooth(Charcoal.accI, s, 'moving');
else
    charAccIS(:,3) = charLowess(accI_clean, s, 'moving');
end

%% -- Method 4: Running median + Lowess ------------------------------------
% Step 1: assign each sample the median accI value within the window.
%         Uses accI_clean so the loop is always NaN-free.
% Step 2: smooth the resulting series with Lowess.
hw4 = round(s / 2);
for i = 1:N
    if i <= hw4
        win = accI_clean(1 : round(s));
    elseif i >= N - hw4
        win = accI_clean(max(1, N - round(s) + 1) : end);
    else
        win = accI_clean(round(i - 0.5*s) : round(i + 0.5*s));
    end
    charAccIS(i, 4) = median(win);
end
if hasCFT
    charAccIS(:,4) = smooth(charAccIS(:,4), s, 'lowess');
else
    charAccIS(:,4) = charLowess(charAccIS(:,4), s, 'lowess');
end

%% -- Method 5: Running mode + Lowess --------------------------------------
% Step 1: divide accI values within the window into 100 equally-spaced
%         bins; assign each sample the centre of the most-populated bin.
%         If multiple bins share the maximum count, their centres are
%         averaged (median of modal centres).
%         Uses accI_clean so the loop is always NaN-free.
% Step 2: smooth the resulting series with Lowess.
hw5  = round(s / 2);
nBin = 100;
for i = 1:N
    if i <= hw5
        win = accI_clean(1 : round(s));
    elseif i >= N - hw5
        win = accI_clean(max(1, N - round(s) + 1) : end);
    else
        win = accI_clean(round(i - 0.5*s) : round(i + 0.5*s));
    end
    [n_mode, x_mode] = charHistCounts(win, nBin);
    modal_centres     = x_mode(n_mode == max(n_mode));
    charAccIS(i, 5)   = median(modal_centres);
end
if hasCFT
    charAccIS(:,5) = smooth(charAccIS(:,5), s, 'lowess');
else
    charAccIS(:,5) = charLowess(charAccIS(:,5), s, 'lowess');
end

%% -- Restore NaN positions (charLowess path only) -------------------------
if ~hasCFT && any(nanMask)
    charAccIS(nanMask, :) = NaN;
end

%% -- Plot (Figure 1, subplot 2) -------------------------------------------
% Only drawn when plotData == 1 and Results.allFigures == 1.
% Subplot 1 is drawn by CharPretreatment; both live in Figure 1.
if plotData == 1 && Results.allFigures == 1

    figure(1);
    subplot(2,1,2)

    h = bar(Charcoal.ybpI, Charcoal.accI, 1);
    set(h, 'facecolor', [.5 .5 .5], 'edgecolor', [.5 .5 .5])
    hold on
    plot(Charcoal.ybpI, charAccIS, 'linewidth', 1.5)
    ylim([0, max(Charcoal.accI)]);

    legend('C_i_n_t_e_r_p_o_l_a_t_e_d', ...
           'Lowess', 'Robust Lowess', 'Moving Average', ...
           'Moving Median', 'Moving Mode')

    xlabel('time (cal. yr BP)')
    ylabel('CHAR (# cm^-^2 yr^-^1)')
    set(gca, 'xlim',    [min(Charcoal.ybp)-100  max(Charcoal.ybp)], ...
             'box',     'off', ...
             'tickdir', 'out', ...
             'xdir',    'rev')
    title(['(b) C_i_n_t_e_r_p_o_l_a_t_e_d and options for a ' ...
           num2str(Smoothing.yr) ' yr C_b_a_c_k_g_r_o_u_n_d'])

end

%% -- Store selected method ------------------------------------------------
Charcoal.accIS = charAccIS(:, Smoothing.method);

end
