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
%   ALGORITHM (lowess / rlowess)
%     For each point i, a neighbourhood of k points centred on i is
%     selected (clipped at record boundaries).  Each neighbour j receives a
%     tricubic distance weight:
%
%         w_j = ( 1 - (|i-j| / d_max)^3 )^3
%
%     A weighted least-squares line is fitted and evaluated at i to give
%     ys(i).  For 'rlowess', residuals from the previous pass are used to
%     form bisquare robustness weights, down-weighting outliers (3 passes).
%
%   NOTE ON SPAN CONVENTION
%     MATLAB smooth() treats span >= 1 as a point count and span < 1 as a
%     fraction.  R lowess() always uses a fraction f in (0,1).  When
%     translating to R, a thin wrapper converting point-count to fraction
%     is recommended (see Modernisation Report, Section 5.1).
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% INPUT HANDLING
if nargin < 3 || isempty(method)
    method = 'lowess';
end

isRow = isrow(y);
y     = y(:);           % work on column vector throughout
n     = numel(y);

%% CONVERT SPAN TO INTEGER WINDOW WIDTH
if span < 1
    k = max(3, round(span * n));    % fraction -> point count
else
    k = max(3, round(span));        % already a point count
end
k  = min(k, n);                     % cannot exceed data length
hw = floor(k / 2);                  % half-window radius in samples

%% MOVING AVERAGE
if strcmpi(method, 'moving')
    ys = movmean(y, k, 'EndPoints', 'shrink');
    if isRow, ys = ys'; end
    return
end

%% LOWESS / RLOWESS: tricubic-weighted local linear regression
nIter = 1;
if strcmpi(method, 'rlowess')
    nIter = 3;
end

ys = y;
rw = ones(n, 1);        % robustness weights (all 1 for plain lowess)

for iter = 1:nIter
    for i = 1:n
        i0 = max(1, i - hw);
        i1 = min(n, i + hw);
        xi = (i0:i1)';
        yi = y(i0:i1);
        ri = rw(i0:i1);

        d    = abs(xi - i);
        dmax = max(d) + eps;
        tri  = (1 - (d ./ dmax) .^ 3) .^ 3;
        w    = tri .* ri;

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
