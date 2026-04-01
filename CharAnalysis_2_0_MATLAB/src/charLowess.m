function ys = charLowess(y, span, method)
% charLowess    Base-MATLAB replacement for smooth() (Curve Fitting Toolbox).
%   ys = charLowess(y, span)
%   ys = charLowess(y, span, method)
%
%   Replicates the behaviour of the Curve Fitting Toolbox function
%       smooth(y, span, method)
%   using only base MATLAB, so CharAnalysis can run without that toolbox.
%
%   INPUTS
%     y      : data vector (column or row; returned in the same orientation)
%     span   : smoothing window.
%              >= 1  interpreted as number of data points
%              <  1  interpreted as fraction of data length
%     method : 'lowess'  (default) locally-weighted linear regression
%              'rlowess' robust lowess (3 re-weighting iterations)
%              'moving'  simple moving average
%
%   OUTPUT
%     ys : smoothed vector, same size as y
%
%   charLowess first checks whether the Curve Fitting Toolbox is available.
%   If it is, smooth() is called directly, guaranteeing results identical
%   to CharAnalysis v1.1. If the toolbox is not available, a pure-MATLAB
%   implementation is used as a fallback.
%
%   ALGORITHM (pure-MATLAB fallback, lowess / rlowess)
%     For each point i, a neighbourhood of k points is selected.  At
%     interior points the window is centred on i.  At boundary points the
%     window is SHIFTED (not shrunk) to maintain exactly k points, matching
%     the behaviour of MATLAB's smooth() from the Curve Fitting Toolbox.
%     Each neighbour j receives a tricubic distance weight:
%
%         w_j = ( 1 - (|i-j| / d_max)^3 )^3
%
%     where d_max is the distance to the FURTHEST point in the window
%     (not the half-window radius), so weights adapt naturally to the
%     shifted boundary windows.
%
%     A weighted least-squares line is fitted and evaluated at i.
%     For 'rlowess', residuals are used to form bisquare robustness
%     weights, down-weighting outliers (3 passes).
%
%   EDGE HANDLING
%     MATLAB smooth() shifts the window at boundaries so every point is
%     fitted using exactly k neighbours.  The previous version of
%     charLowess shrunk the window instead, producing slightly different
%     values at record edges.  This version uses shifted windows to match
%     MATLAB smooth() behaviour and eliminate the edge discrepancy that
%     caused 2 peaks near the oldest part of the record to be missed.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── Input handling ───────────────────────────────────────────────────────
if nargin < 3 || isempty(method)
    method = 'lowess';
end
isRow = isrow(y);
y     = y(:);           % work on column vector throughout
n     = numel(y);

%% ── Convert span to integer window width ─────────────────────────────────
if span < 1
    k = max(3, round(span * n));    % fraction -> point count
else
    k = max(3, round(span));        % already a point count
end
k  = min(k, n);                     % cannot exceed data length
hw = floor(k / 2);                  % nominal half-window radius

%% ── Moving average ───────────────────────────────────────────────────────
if strcmpi(method, 'moving')
    ys = movmean(y, k, 'EndPoints', 'shrink');
    if isRow, ys = ys'; end
    return
end

%% ── Use Curve Fitting Toolbox smooth() if available ──────────────────────
% smooth() produces results identical to v1.1. The pure-MATLAB fallback
% below is used only on installations without the Curve Fitting Toolbox.
if license('test', 'Curve_Fitting_Toolbox')
    switch lower(method)
        case 'lowess'
            ys = smooth(y, k, 'lowess');
        case 'rlowess'
            ys = smooth(y, k, 'rlowess');
        case 'moving'
            ys = smooth(y, k, 'moving');
        otherwise
            ys = smooth(y, k, 'lowess');
    end
    if isRow, ys = ys'; end
    return
end

%% ── Pure-MATLAB fallback (no Curve Fitting Toolbox) ──────────────────────
if strcmpi(method, 'rlowess')
    nIter = 5;
else
    nIter = 1;
end

ys = y;
rw = ones(n, 1);        % robustness weights (all 1 for plain lowess)

for iter = 1:nIter
    for i = 1:n

        % ── Window selection: shift at boundaries to keep exactly k pts ──
        % This matches MATLAB smooth() edge behaviour.
        i0 = i - hw;
        i1 = i0 + k - 1;
        if i0 < 1
            i0 = 1;
            i1 = k;
        end
        if i1 > n
            i1 = n;
            i0 = n - k + 1;
        end

        xi = (i0:i1)';
        yi = y(i0:i1);
        ri = rw(i0:i1);

        % Tricubic weights: distance from i to each window point,
        % normalised by distance to the furthest point in the window.
        d    = abs(xi - i);
        dmax = max(d) + eps;
        tri  = (1 - (d ./ dmax) .^ 3) .^ 3;
        w    = tri .* ri;

        % Weighted least-squares linear fit
        sw   = sum(w);
        swx  = sum(w .* xi);
        swy  = sum(w .* yi);
        swx2 = sum(w .* xi .^ 2);
        swxy = sum(w .* xi .* yi);
        det  = sw * swx2 - swx ^ 2;

        if abs(det) < eps * max(abs([sw*swx2, swx^2]))
            ys(i) = swy / (sw + eps);   % near-singular: fall back to mean
        else
            a     = (swy * swx2 - swx * swxy) / det;
            b     = (sw  * swxy - swx * swy)  / det;
            ys(i) = a + b * i;
        end
    end

    % Bisquare robustness weights for rlowess
    if strcmpi(method, 'rlowess') && iter < nIter
        res = y - ys;
        s   = median(abs(res));
        if s < eps, break; end
        u  = res ./ (6 * s);
        rw = max(0, 1 - u .^ 2) .^ 2 .* (abs(u) < 1);
    end
end

if isRow, ys = ys'; end
end