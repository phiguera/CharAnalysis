function z_Compare_CharAnalysis_Results_V1_V2(paramsFile, referenceFile, varargin)
% CharRegressionTest   Regression test for CharAnalysis v2.0.
%   CharRegressionTest(paramsFile, referenceFile)
%   CharRegressionTest(paramsFile, referenceFile, 'Name', Value, ...)
%
%   Runs CharAnalysis on the specified input file and compares key output
%   columns against a validated reference CSV. Reports pass/fail for each
%   assertion and prints a summary.
%
%   REQUIRED INPUTS
%     paramsFile    : path to the CharAnalysis parameter file (.csv or .xls)
%     referenceFile : path to a validated v1.1 charResults CSV file
%
%   OPTIONAL NAME-VALUE PAIRS
%     'PeakTol'     : maximum allowed difference in peak counts (default 0)
%     'ThreshTol'   : maximum allowed mean threshold difference (default 0.001)
%     'MagTol'      : maximum allowed total magnitude difference % (default 2.0)
%     'Verbose'     : true/false, print per-sample details (default false)
%
%   EXAMPLE
%     CharRegressionTest('CO_charParams.csv', 'CO_charResults_V_1_1.csv')
%     CharRegressionTest('CO_charParams.csv', 'CO_charResults_V_1_1.csv', ...
%         'ThreshTol', 0.002, 'Verbose', true)
%
%   OUTPUT COLUMNS (referenceFile and charResults must share these columns)
%     Col  2 : age Top_i (yr BP)
%     Col  9 : thresh 1
%     Col 10 : thresh 2
%     Col 11 : thresh 3
%     Col 12 : thresh FinalPos
%     Col 16 : peaks 1
%     Col 17 : peaks 2
%     Col 18 : peaks 3
%     Col 19 : peaks Final
%     Col 20 : peaks Insig.
%     Col 21 : peak Mag
%     Col 22 : smPeak Frequ
%
%   CharAnalysis v2.0

%% ?? Parse inputs ?????????????????????????????????????????????????????????
p = inputParser;
addRequired(p,  'paramsFile',    @ischar);
addRequired(p,  'referenceFile', @ischar);
addParameter(p, 'PeakTol',   0,     @isnumeric);
addParameter(p, 'ThreshTol', 0.001, @isnumeric);
addParameter(p, 'MagTol',    2.0,   @isnumeric);
addParameter(p, 'Verbose',   false, @islogical);
parse(p, paramsFile, referenceFile, varargin{:});

peakTol   = p.Results.PeakTol;
threshTol = p.Results.ThreshTol;
magTol    = p.Results.MagTol;
verbose   = p.Results.Verbose;

%% ?? Header ???????????????????????????????????????????????????????????????
fprintf('\n')
fprintf('=================================================================\n')
fprintf('  CharAnalysis v2.0 Regression Test\n')
fprintf('  %s\n', datetime('now','Format','yyyy-MM-dd HH:mm:ss'))
fprintf('=================================================================\n')
fprintf('  Params file : %s\n', paramsFile)
fprintf('  Reference   : %s\n', referenceFile)
fprintf('  PeakTol     : %d\n', peakTol)
fprintf('  ThreshTol   : %.4f\n', threshTol)
fprintf('  MagTol      : %.1f%%\n', magTol)
fprintf('-----------------------------------------------------------------\n')

%% ?? Run CharAnalysis ?????????????????????????????????????????????????????
fprintf('\nRunning CharAnalysis...\n')
try
    results = CharAnalysis(paramsFile);
catch ME
    fprintf('FAIL: CharAnalysis threw an error:\n')
    fprintf('      %s\n', ME.message)
    fprintf('      (in %s, line %d)\n', ME.stack(1).name, ME.stack(1).line)
    printSummary(0, 0)
    return
end
fprintf('CharAnalysis completed successfully.\n\n')

%% ?? Locate v2.0 results file ?????????????????????????????????????????????
% CharAnalysis saves results to <site>_charResults.csv in the working dir.
[~, refName] = fileparts(referenceFile);
% Strip known v1.1 suffix patterns to recover site name
siteName = regexprep(refName, '(_V_1_1|_charResults.*)', '', 'ignorecase');
v20File  = [cd filesep siteName '_charResults.csv'];

if ~exist(v20File, 'file')
    % Fall back: look for any *_charResults.csv in current directory
    candidates = dir(fullfile(cd, '*_charResults.csv'));
    % Exclude the reference file itself
    candidates = candidates(~strcmpi({candidates.name}, ...
                                      [refName '.csv']));
    if isempty(candidates)
        fprintf('FAIL: Could not find v2.0 charResults CSV in %s\n', cd)
        printSummary(0, 0)
        return
    end
    v20File = fullfile(cd, candidates(1).name);
end
fprintf('v2.0 results file: %s\n\n', v20File)

%% ?? Load both result files ???????????????????????????????????????????????
try
    ref = readmatrix(referenceFile, 'NumHeaderLines', 1);
    v20 = readmatrix(v20File,       'NumHeaderLines', 1);
catch ME
    fprintf('FAIL: Could not read results files:\n      %s\n', ME.message)
    printSummary(0, 0)
    return
end

nRef = size(ref, 1);
nV20 = size(v20, 1);
if nRef ~= nV20
    fprintf('WARN: Row count differs — reference=%d, v2.0=%d\n', nRef, nV20)
end
nRows = min(nRef, nV20);

%% ?? Column definitions ???????????????????????????????????????????????????
peakCols   = [16 17 18 19];
peakNames  = {'peaks1','peaks2','peaks3','peaksFinal'};
threshCols = [9 10 11 12];
threshNames = {'thresh1','thresh2','thresh3','threshFinalPos'};
magCol     = 21;
freqCol    = 22;

%% ?? Run assertions ???????????????????????????????????????????????????????
nPass = 0;
nFail = 0;

%% ?? Assertion 1: Peak counts ?????????????????????????????????????????????
fprintf('--- Peak count comparison ---\n')
fprintf('%-16s  %6s  %6s  %6s  %s\n', 'Column','v1.1','v2.0','Diff','Result')
for c = 1:4
    col  = peakCols(c);
    n11  = nansum(ref(1:nRows, col));
    n20  = nansum(v20(1:nRows, col));
    diff = abs(n20 - n11);
    pass = diff <= peakTol;
    if pass; nPass = nPass+1; else; nFail = nFail+1; end
    fprintf('%-16s  %6d  %6d  %6d  %s\n', peakNames{c}, n11, n20, ...
        n20-n11, passStr(pass))
end

% Per-row peak final differences
finalDiffRows = find(ref(1:nRows, 19) ~= v20(1:nRows, 19));
fprintf('\nRows where peaksFinal differs: %d of %d\n', ...
    length(finalDiffRows), nRows)
if ~isempty(finalDiffRows) && verbose
    fprintf('Ages (yr BP) at differing rows:\n')
    disp(ref(finalDiffRows, 2))
end
fprintf('\n')

%% ?? Assertion 2: Threshold mean values ???????????????????????????????????
fprintf('--- Threshold comparison (mean values) ---\n')
fprintf('%-18s  %10s  %10s  %10s  %s\n', ...
    'Column','v1.1 mean','v2.0 mean','Diff','Result')
for c = 1:4
    col  = threshCols(c);
    m11  = nanmean(ref(1:nRows, col));
    m20  = nanmean(v20(1:nRows, col));
    diff = abs(m20 - m11);
    pass = diff <= threshTol;
    if pass; nPass = nPass+1; else; nFail = nFail+1; end
    fprintf('%-18s  %10.4f  %10.4f  %10.4f  %s\n', ...
        threshNames{c}, m11, m20, diff, passStr(pass))
end

% Sample-by-sample threshold spread
threshDiff = abs(ref(1:nRows,12) - v20(1:nRows,12));
fprintf('\nSample-by-sample threshFinalPos:\n')
fprintf('  Samples differing (> 1e-6): %d of %d\n', ...
    sum(threshDiff > 1e-6), nRows)
fprintf('  Mean absolute difference  : %.4f\n', nanmean(threshDiff))
fprintf('  Max absolute difference   : %.4f\n', nanmax(threshDiff))
fprintf('\n')

%% ?? Assertion 3: Peak magnitude ??????????????????????????????????????????
fprintf('--- Peak magnitude comparison ---\n')
tot11  = nansum(ref(1:nRows, magCol));
tot20  = nansum(v20(1:nRows, magCol));
pctDiff = 100 * abs(tot20 - tot11) / (tot11 + eps);
pass   = pctDiff <= magTol;
if pass; nPass = nPass+1; else; nFail = nFail+1; end
fprintf('Total magnitude:  v1.1 = %.2f,  v2.0 = %.2f,  diff = %.1f%%  %s\n', ...
    tot11, tot20, pctDiff, passStr(pass))

mn11   = nanmean(ref(ref(1:nRows,magCol)>0, magCol));
mn20   = nanmean(v20(v20(1:nRows,magCol)>0, magCol));
fprintf('Mean per peak:    v1.1 = %.4f,  v2.0 = %.4f\n', mn11, mn20)
fprintf('\n')

%% ?? Assertion 4: Smoothed fire frequency ????????????????????????????????
fprintf('--- Smoothed fire frequency comparison ---\n')
ff11   = ref(1:nRows, freqCol);
ff20   = v20(1:nRows, freqCol);
ffDiff = abs(ff11 - ff20);
maxFF  = nanmax(ffDiff);
meanFF = nanmean(ffDiff);
pass   = meanFF <= threshTol;
if pass; nPass = nPass+1; else; nFail = nFail+1; end
fprintf('Mean absolute difference: %.4f  %s\n', meanFF, passStr(pass))
fprintf('Max  absolute difference: %.4f\n', maxFF)
fprintf('\n')

%% ?? Assertion 5: No NaN in key output columns ????????????????????????????
fprintf('--- NaN check (v2.0 output) ---\n')
checkCols = [peakCols, threshCols, magCol];
checkNames = [peakNames, threshNames, {'peakMag'}];
allClean = true;
for c = 1:length(checkCols)
    nNaN = sum(isnan(v20(1:nRows, checkCols(c))));
    if nNaN > 0
        fprintf('%-18s  %d NaN values  FAIL\n', checkNames{c}, nNaN)
        nFail = nFail + 1;
        allClean = false;
    end
end
if allClean
    fprintf('All key output columns are NaN-free  %s\n', passStr(true))
    nPass = nPass + 1;
end
fprintf('\n')

%% ?? Summary ??????????????????????????????????????????????????????????????
printSummary(nPass, nFail)

end % CharRegressionTest

%% ?? Local helpers ????????????????????????????????????????????????????????
function s = passStr(tf)
    if tf
        s = 'PASS';
    else
        s = 'FAIL ***';
    end
end

function printSummary(nPass, nFail)
    nTotal = nPass + nFail;
    fprintf('=================================================================\n')
    fprintf('  SUMMARY:  %d of %d assertions passed', nPass, nTotal)
    if nFail == 0
        fprintf('  --  ALL PASS\n')
    else
        fprintf('  --  %d FAILED\n', nFail)
    end
    fprintf('=================================================================\n\n')
end