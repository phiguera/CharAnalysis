function CharValidateParams(charData, Pretreatment, Smoothing, ...
    PeakAnalysis, Results)
% CharValidateParams   Validate CharAnalysis input parameters before run.
%   CharValidateParams(charData, Pretreatment, Smoothing, PeakAnalysis,
%       Results)
%
%   Checks all user-supplied parameters for internal consistency and
%   plausible ranges. Throws a descriptive error and halts execution if
%   any check fails. Emits warnings for non-fatal but unusual settings.
%
%   Called from CharAnalysis.m immediately after CharParameters returns,
%   replacing the two inline if/error blocks that existed in v1.1.
%
%   Checks performed
%   ────────────────
%   charData
%     1. At least 6 columns present (cmTop cmBot ageTop ageBot charVol charCount)
%     2. No NaN values in the age columns (cols 3 and 4)
%     3. Ages are monotonically non-decreasing (ageTop col 3)
%
%   Pretreatment
%     4. zoneDiv has at least 2 values
%     5. zoneDiv is strictly ascending
%     6. yrInterp is non-negative (0 = auto)
%     7. zoneDiv(end) does not exceed the bottom age of the last raw sample
%        (warning only -- v2.0 marks overrun samples NaN, but results will
%        differ from v1.1 which silently filled them with zero CHAR)
%
%   Smoothing
%     8. Smoothing window (yr) is shorter than the record length
%     9. Smoothing method is in the range 1-5
%
%   PeakAnalysis
%     10. threshType is 1 or 2
%     11. threshMethod is 1, 2, or 3
%     12. Locally defined threshold cannot be user-defined
%         (threshMethod == 1 && threshType == 2)
%     13. threshValues all in (0, 1) when threshMethod > 1
%     14. cPeak is 1 (residuals) or 2 (ratios)
%     15. minCountP is in [0, 1]
%     16. peakFrequ is positive
%
%   CharAnalysis v2.0

%% ── 1. charData column count ─────────────────────────────────────────────
if size(charData, 2) < 6
    error('CharValidateParams: charData must have at least 6 columns (cmTop, cmBot, ageTop, ageBot, charVol, charCount). Found %d.', ...
        size(charData, 2))
end

%% ── 2. No NaN ages ───────────────────────────────────────────────────────
if any(isnan(charData(:,3))) || any(isnan(charData(:,4)))
    error('CharValidateParams: charData contains NaN values in the age columns (cols 3–4). All sample ages must be specified.')
end

%% ── 3. Ages monotonically non-decreasing ─────────────────────────────────
ageTops = charData(:,3);
if any(diff(ageTops) < 0)
    badRows = find(diff(ageTops) < 0);
    error('CharValidateParams: ageTop values (col 3) are not monotonically non-decreasing. First violation at row %d.', ...
        badRows(1))
end

%% ── 4. zoneDiv has at least 2 values ─────────────────────────────────────
if numel(Pretreatment.zoneDiv) < 2
    error('CharValidateParams: zoneDiv must contain at least 2 values (start and end of record).')
end

%% ── 5. zoneDiv strictly ascending ───────────────────────────────────────
if any(diff(Pretreatment.zoneDiv) <= 0)
    error('CharValidateParams: zoneDiv must be strictly ascending (youngest age first).')
end

%% ── 6. yrInterp non-negative ─────────────────────────────────────────────
if Pretreatment.yrInterp < 0
    error('CharValidateParams: yrInterp must be >= 0 (0 = use median sample resolution).')
end

%% ── 7. zoneDiv(end) within data extent (warning) ────────────────────────
% If zoneDiv(end) exceeds the bottom age of the last raw sample, the
% vectorized proportion matrix in CharPretreatment will find no overlapping
% raw data for those terminal intervals and mark them NaN. This causes
% charBkg to differ from v1.1, which silently filled them with zero CHAR.
% Set zoneDiv(end) to no greater than the bottom age of the last raw sample.
lastAgeBotInData = charData(end, 4);
if Pretreatment.zoneDiv(end) > lastAgeBotInData
    warning(['CharValidateParams: zoneDiv(end) (%g yr BP) exceeds the ' ...
             'bottom age of the last raw sample (%g yr BP). ' ...
             'Terminal interpolated samples beyond the data will be ' ...
             'set to NaN. Set zoneDiv(end) <= %g to avoid this.'], ...
        Pretreatment.zoneDiv(end), lastAgeBotInData, lastAgeBotInData)
end

%% ── 8. Smoothing window shorter than record ──────────────────────────────
recordLength = max(Pretreatment.zoneDiv) - min(Pretreatment.zoneDiv);
smoothYr     = Smoothing.yr(end);   % Use the longest window specified.
if smoothYr >= recordLength
    error('CharValidateParams: smoothing window (%d yr) must be shorter than the record length (%d yr).', ...
        smoothYr, recordLength)
end

%% ── 9. Smoothing method in range ─────────────────────────────────────────
if ~ismember(Smoothing.method, 1:5)
    error('CharValidateParams: smoothing method must be 1–5. Got %d.', ...
        Smoothing.method)
end

%% ── 10. threshType valid ──────────────────────────────────────────────────
if ~ismember(PeakAnalysis.threshType, [1 2])
    error('CharValidateParams: threshType must be 1 (global) or 2 (local). Got %d.', ...
        PeakAnalysis.threshType)
end

%% ── 11. threshMethod valid ───────────────────────────────────────────────
if ~ismember(PeakAnalysis.threshMethod, [1 2 3])
    error('CharValidateParams: threshMethod must be 1, 2, or 3. Got %d.', ...
        PeakAnalysis.threshMethod)
end

%% ── 12. Local threshold cannot be user-defined ───────────────────────────
if PeakAnalysis.threshMethod == 1 && PeakAnalysis.threshType == 2
    error('CharValidateParams: a locally defined threshold (threshType = 2) cannot be user-defined (threshMethod = 1). Change input parameters.')
end

%% ── 13. threshValues in (0,1) when threshMethod > 1 ─────────────────────
if PeakAnalysis.threshMethod > 1
    tv = PeakAnalysis.threshValues(1:4);
    if any(tv <= 0) || any(tv >= 1)
        error('CharValidateParams: all threshValues must be in the open interval (0, 1) when threshMethod = %d. Got [%s].', ...
            PeakAnalysis.threshMethod, num2str(tv))
    end
end

%% ── 14. cPeak valid ──────────────────────────────────────────────────────
if ~ismember(PeakAnalysis.cPeak, [1 2])
    error('CharValidateParams: cPeak must be 1 (residuals) or 2 (ratios). Got %d.', ...
        PeakAnalysis.cPeak)
end

%% ── 15. minCountP in [0,1] ───────────────────────────────────────────────
if PeakAnalysis.minCountP < 0 || PeakAnalysis.minCountP > 1
    error('CharValidateParams: minCountP must be in [0, 1]. Got %g.', ...
        PeakAnalysis.minCountP)
end

%% ── 16. peakFrequ positive ───────────────────────────────────────────────
if PeakAnalysis.peakFrequ <= 0
    error('CharValidateParams: peakFrequ must be positive. Got %g.', ...
        PeakAnalysis.peakFrequ)
end

%% ── Non-fatal warnings ───────────────────────────────────────────────────
% Warn if smoothing window is very short relative to record length.
if smoothYr < 0.05 * recordLength
    warning('CharValidateParams: smoothing window (%d yr) is less than 5%% of the record length (%d yr). Consider a longer window.', ...
        smoothYr, recordLength)
end

% Warn if local threshold window may contain fewer than 30 samples.
if PeakAnalysis.threshType == 2 && Pretreatment.yrInterp > 0
    samplesInWindow = smoothYr / Pretreatment.yrInterp;
    if samplesInWindow < 30
        warning('CharValidateParams: smoothing window (%d yr) at yrInterp = %d yr gives ~%d samples per local threshold window. At least 30 are recommended.', ...
            smoothYr, Pretreatment.yrInterp, floor(samplesInWindow))
    end
end

disp('      Parameter validation passed.')
end