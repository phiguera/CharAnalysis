function [Charcoal, CharThresh] = CharPeakID(Charcoal, Pretreatment, ...
        PeakAnalysis, CharThresh)
% CharPeakID    Identify charcoal peaks based on threshold values.
%   [Charcoal, CharThresh] = CharPeakID(Charcoal, Pretreatment,
%       PeakAnalysis, CharThresh)
%
%   Identifies interpolated samples where Cpeak exceeds each threshold
%   value, removes consecutive exceedances so that only one sample per
%   event is retained, then screens the identified peaks against a
%   minimum-count criterion.
%
%   INPUTS
%     Charcoal     : struct containing peak, accI, countI, volI, ybpI
%     Pretreatment : struct  (yrInterp)
%     PeakAnalysis : struct  (threshType, threshValues, minCountP)
%     CharThresh   : struct  (possible, pos, neg, minCountP)
%
%   OUTPUTS  (fields added to input structs)
%     Charcoal.charPeaks       : [N x T] binary peak matrix (1 = peak)
%     Charcoal.charPeaksThresh : [N x T] threshold value at each peak
%     Charcoal.peaksTotal      : [1 x T] total peak count per threshold
%     Charcoal.threshFRI       : fire-return intervals per threshold
%     CharThresh.minCountP     : [N x T] p-values from minimum-count test
%
%   v2.0 changes vs v1.1
%     - Peak flagging: O(N x T) double for-loop replaced by a single
%       vectorised broadcast comparison.  Produces identical results.
%     - Consecutive-peak removal: two sequential O(N x T) for-loops
%       replaced by a diff-based column-wise operation.  Produces
%       identical results (see note on trailing-element behaviour below).
%     - charPeaksThresh construction: O(N x T) for-loop replaced by a
%       single broadcast multiply.
%     - Minimum-count analysis loop retained as-is.  It iterates over
%       identified peaks (~10-200 per record), not over all N x T samples,
%       so vectorising it would add complexity with negligible speed gain.
%
%   NOTE ON CONSECUTIVE-PEAK REMOVAL
%     The v1.1 algorithm (and v2.0) retains the LAST sample of each
%     consecutive run of above-threshold values (highest array index =
%     oldest age in cal yr BP).  For a run of flagged samples the loop
%     halved every element whose successor was also flagged, leaving the
%     trailing element as the only one still at value 2 after the pass.
%     The comment in v1.1 said "keep first peak as 2" but the implemented
%     behaviour was the reverse.  v2.0 preserves the implemented behaviour
%     exactly and documents it correctly.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── LOCAL VARIABLES ───────────────────────────────────────────────────────
r         = Pretreatment.yrInterp;
N         = length(Charcoal.peak);
alphaPeak = PeakAnalysis.minCountP;

%% ── BUILD THRESHOLD VALUE MATRIX ─────────────────────────────────────────
%
%   thresholdValues is [N x T]:
%     Global threshold  - constant columns (one value per column, same for
%                         all rows); built from CharThresh.possible.
%     Local threshold   - CharThresh.pos is already [N x T].

if PeakAnalysis.threshType == 1

    % Global: restrict to the positive-valued portion of the threshold
    % grid, then replicate into an [N x T] matrix.
    CharThreshPossible = CharThresh.possible(CharThresh.possible > 0);
    nThresholds        = length(CharThreshPossible);
    thresholdValues    = repmat(CharThreshPossible, N, 1);   % [N x T]

else

    % Local: threshold varies by sample; [N x T] already.
    thresholdValues = CharThresh.pos;
    nThresholds     = size(thresholdValues, 2);

end

%% ── VECTORISED PEAK FLAGGING ─────────────────────────────────────────────
%
%   v1.1 used a double for-loop:
%       for i = 1:N
%           for j = 1:nThresholds
%               if Charcoal.peak(i) > thresholdValues(i,j)
%                   Charcoal.charPeaks(i,j) = 2;
%               ...
%
%   v2.0: Charcoal.peak is [N x 1]; thresholdValues is [N x T].
%   MATLAB broadcasts the comparison across all T columns simultaneously.
%   Exceedances are marked as 2 (not 1) so the consecutive-removal step
%   can use value 1 as a "downgraded" marker distinct from "no peak" (0).

Charcoal.charPeaks = 2 * double(Charcoal.peak > thresholdValues);

%% ── VECTORISED CONSECUTIVE-PEAK REMOVAL ──────────────────────────────────
%
%   v1.1 walked forward through the array twice:
%     Pass 1: if sample i AND sample i+1 are both flagged, halve sample i
%             (2 -> 1).  The last sample of each run therefore stays at 2.
%     Pass 2: anything < 2 -> 0;  anything >= 2 -> 1.
%   Net result: the LAST sample of each consecutive run is kept as 1.
%
%   v2.0 uses diff to find run-end positions in a single pass per column:
%     diff(pk)(i) < 0  is true when pk(i) = 1 and pk(i+1) = 0,
%     i.e. sample i is the last flagged sample before a gap.
%     Appending pk(end) handles a run that reaches the final sample.

for j = 1:nThresholds
    pk    = Charcoal.charPeaks(:, j) > 0;          % logical [N x 1]
    isEnd = [diff(pk) < 0; pk(end)];               % true at run-end samples
    Charcoal.charPeaks(:, j) = double(isEnd);
end

%% ── charPeaksThresh: THRESHOLD VALUE AT EACH IDENTIFIED PEAK ─────────────
%
%   Used only for graphing.  Records the threshold value that each peak
%   was identified against, so downstream plot code can draw the appropriate
%   marker height.
%
%   v1.1 used an O(N x T) for-loop using thresholdValues(1,:) for all rows.
%   v2.0: single broadcast multiply.  thresholdValues(1,:) is [1 x T];
%   Charcoal.charPeaks is [N x T]; MATLAB broadcasts the row across N rows.

Charcoal.charPeaksThresh = Charcoal.charPeaks .* thresholdValues(1, :);

%% ── MINIMUM-COUNT ANALYSIS ───────────────────────────────────────────────
%
%   For each identified peak, tests whether the maximum charcoal count
%   in a window around the peak is statistically greater than the minimum
%   count using the Detre-White / Shuie-Bain test for unequal sample
%   volumes (Gavin 2006).  Peaks that do not pass are flagged in
%   CharThresh.minCountP and removed from Charcoal.charPeaks below.
%
%   This loop iterates over identified peaks (typically 10-200 per record),
%   not over all N samples, so it is not a performance bottleneck.

mcWindow = round(150 / r) * r;     % [yr] search window around each peak
d        = zeros(N, nThresholds);
CharThresh.minCountP = NaN(N, nThresholds);

for j = 1:nThresholds

    peakIndex = find(Charcoal.charPeaks(:, j));

    if length(peakIndex) <= 1
        continue        % need at least 2 peaks to define inter-peak windows
    end

    for i = 1:length(peakIndex)

        peakYr = Charcoal.ybpI(peakIndex);

        % Time window bounded by mcWindow on each side
        windowTime = [ max(Charcoal.ybpI(Charcoal.ybpI <= peakYr(i) + mcWindow)), ...
                       min(Charcoal.ybpI(Charcoal.ybpI >= peakYr(i) - mcWindow)) ];
        windowTime_in = [ find(Charcoal.ybpI == windowTime(1)), ...
                          find(Charcoal.ybpI == windowTime(2)) ];

        % Narrow window to adjacent peaks when they fall closer than mcWindow
        if i == 1
            windowPeak_in = [ find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i+1))), ...
                              find(Charcoal.ybpI == windowTime(2)) ];
        elseif i == length(peakIndex)
            windowPeak_in = [ find(Charcoal.ybpI == windowTime(1)), ...
                              find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i-1))) ];
        else
            windowPeak_in = [ find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i+1))), ...
                              find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i-1))) ];
        end

        if windowTime_in(1) > windowPeak_in(1)
            windowTime_in(1) = windowPeak_in(1);
        end
        if windowTime_in(2) < windowPeak_in(2)
            windowTime_in(2) = windowPeak_in(2);
        end

        windowSearch = windowTime_in;

        % Maximum count AFTER the peak (younger side, lower indices)
        countMax = round(max(Charcoal.countI(windowSearch(2) : peakIndex(i))));
        countMaxIn = windowSearch(2) - 1 + ...
            find(round(Charcoal.countI(windowSearch(2) : peakIndex(i))) ...
                 == countMax, 1);

        % Minimum count BEFORE the peak (older side, higher indices)
        countMin = round(min(Charcoal.countI(peakIndex(i) : windowSearch(1))));
        countMinIn = peakIndex(i) - 1 + ...
            find(round(Charcoal.countI(peakIndex(i) : windowSearch(1))) ...
                 == countMin, 1);

        volMax = Charcoal.volI(countMaxIn);
        volMin = Charcoal.volI(countMinIn);

        % Shuie-Bain (1982) expansion of Detre-White (1970) test statistic
        % for unequal sediment volumes (from Gavin 2006 / Charster)
        d(peakIndex(i), j) = ...
            ( abs(countMin - (countMin + countMax) * ...
                  (volMin / (volMin + volMax))) - 0.5 ) / ...
            sqrt( (countMin + countMax) * ...
                  (volMin  / (volMin + volMax)) * ...
                  (volMax  / (volMin + volMax)) );

        % p-value: t-distribution at 1e10 df is numerically the standard
        % normal (df -> inf limit); tcdf from Statistics Toolbox used here.
        CharThresh.minCountP(peakIndex(i), j) = 1 - tcdf(d(peakIndex(i), j), 1e10);

    end % each peak
end % each threshold

%% ── REMOVE PEAKS FAILING THE MINIMUM-COUNT CRITERION ────────────────────

for j = 1:nThresholds
    failIdx = find( Charcoal.charPeaks(:, j)    >  0 & ...
                    CharThresh.minCountP(:, j)   > alphaPeak );
    Charcoal.charPeaks(failIdx, j)       = 0;
    Charcoal.charPeaksThresh(failIdx, j) = 0;
end

%% ── FIRE-RETURN INTERVAL SENSITIVITY INDICES ────────────────────────────

for j = 1:nThresholds
    Charcoal.peaksTotal(j) = sum(Charcoal.charPeaks(:, j));

    inFRI = diff(Charcoal.ybpI(Charcoal.charPeaks(:, j) > 0));
    if ~isempty(inFRI)
        Charcoal.threshFRI(1:length(inFRI), j) = inFRI;
    end
end

end
