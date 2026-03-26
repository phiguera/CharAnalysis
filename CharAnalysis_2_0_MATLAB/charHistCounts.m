function [n, x] = charHistCounts(y, bins)
% charHistCounts    Modern replacement for the two-output form of hist().
%   [n, x] = charHistCounts(y, bins)
%
%   hist() was deprecated as a recommended function in MATLAB R2014b.
%   The recommended replacement, histcounts(), uses bin *edges* rather than
%   bin *centres*, so existing CharAnalysis code that relied on hist()
%   cannot simply swap in histcounts() without adjusting indices.
%   This wrapper preserves the hist() calling convention so each call site
%   in CharAnalysis requires only a one-word change.
%
%   INPUTS
%     y    : data vector
%     bins : bin specification — matches hist() behaviour exactly:
%              scalar  -> number of equally spaced bins (auto range)
%              vector  -> explicit bin centres (uniform spacing assumed)
%
%   OUTPUTS
%     n : count per bin  (row vector), same length as x
%     x : bin centres    (row vector)
%
%   EXAMPLES
%     % Scalar bins (like hist(y, 50))
%     [n, x] = charHistCounts(y, 50);
%
%     % Vector of centres (like hist(y, centers))
%     centers = linspace(0, 10, 30);
%     [n, x] = charHistCounts(y, centers);
%
%     % Normalised histogram (proportion per bin)
%     bar(x, n / sum(n));
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

y = y(:);

if isscalar(bins)
    %% Scalar: ask histcounts to auto-select edges over the data range
    [n, edges] = histcounts(y, bins);
    x = edges(1:end-1) + diff(edges) / 2;   % recover centres from edges

else
    %% Vector of centres: convert to edges, then call histcounts
    bins   = bins(:)';                       % ensure row vector
    d      = diff(bins);

    % Build edges halfway between adjacent centres; extend at both ends
    inner  = bins(1:end-1) + d / 2;
    left   = bins(1)   - d(1)   / 2;
    right  = bins(end) + d(end) / 2;
    edges  = [left, inner, right];

    n = histcounts(y, edges);
    x = bins;
end

end
