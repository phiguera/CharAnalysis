function z_test_params_writer()
% z_test_params_writer   Tests for the parameter-file reader and writers.
%   z_test_params_writer()
%
%   MATLAB counterpart of CharAnalysis_2_0_R/tests/test_params_writer.R.
%   Run from any folder; all files are written to a temporary folder.
%
%   1. Template: CharParamsReadTable reads the 26-row template.
%   2. Round trip: for each bundled site, read -> write a new file ->
%      re-read gives identical parameters (CO also: identical results).
%   3. A new file is never written over an existing one.
%   4. Revision: unchanged values do not rewrite the file; previous values
%      move to a time-stamped log column; hand-entered log columns are
%      preserved; the 10-log limit is enforced; a full 16-column file is
%      read correctly by CharParameters and runs.
%   5. Files not named <site>_charParams.csv are rejected.
%   6. CharAnalysis called with a full path writes results to the folder
%      holding the parameter file and returns to the original folder.
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

matDir   = fileparts(mfilename('fullpath'));
addpath(fullfile(matDir, 'src'));
template = fullfile(matDir, 'templates', 'template_charParams.csv');
startDir = pwd;
wState   = warning;
cleanup  = onCleanup(@() restoreState(startDir, wState)); %#ok<NASGU>

work = tempname;
mkdir(work);
fprintf('Working folder: %s\n\n', work)
nPass = 0;

%% 1. Template --------------------------------------------------------------
tmpl = CharParamsReadTable(template);
assert(size(tmpl, 1) == 26 && size(tmpl, 2) == 7, 'Template should be 26 x 7.')
assert(strcmp(tmpl{2, 2}, 'zoneDiv') && strcmp(tmpl{26, 2}, 'allFigures'), ...
    'Template labels not read as expected.')
nPass = pass(nPass, 'template read: 26 rows x 7 columns, labels in place');

%% 2. Round trip ------------------------------------------------------------
sites = {'CO', 'CH10', 'SI17', 'RA07', 'TL06'};
for k = 1:numel(sites)
    s    = sites{k};
    orig = fullfile(matDir, [s '_charParams.csv']);
    [~, P1, S1, A1, R1] = CharParameters(orig);
    v   = paramsVector(P1, S1, A1, R1);
    raw = str2double(CharParamsReadTable(orig));
    raw = raw(2:26, 3);
    assert(isequal(v, raw), '%s: parameter vector differs from file column 3.', s)

    d = fullfile(work, s); mkdir(d);
    copyfile(fullfile(matDir, [s '_charData.csv']), d);
    newFile = fullfile(d, [s '_charParams.csv']);
    CharParamsNew(v, newFile, template);
    [~, P2, S2, A2, R2] = CharParameters(newFile);
    assert(isequal(P1, P2) && isequal(S1, S2) && isequal(A1, A2) && ...
           isequal(R1, R2), '%s: parameters differ after round trip.', s)
    nPass = pass(nPass, ['round trip, identical parameters: ' s]);
end

% Identical analysis results (CO): original file vs. rewritten file
dA = fullfile(work, 'CO_orig'); mkdir(dA);
copyfile(fullfile(matDir, 'CO_charData.csv'), dA);
copyfile(fullfile(matDir, 'CO_charParams.csv'), dA);
rng(1); o1 = evalcRun(fullfile(dA, 'CO_charParams.csv'));
rng(1); o2 = evalcRun(fullfile(work, 'CO', 'CO_charParams.csv'));
close all
assert(isequaln(o1.Post.charResults, o2.Post.charResults), ...
    'CO: charResults differ between original and rewritten parameter files.')
nPass = pass(nPass, 'round trip, identical CharAnalysis results: CO');

%% 3. No overwrite ----------------------------------------------------------
coNew  = fullfile(work, 'CO', 'CO_charParams.csv');
before = fileread(coNew);
expectError(@() CharParamsNew(v0of(coNew), coNew, template), 'will not be overwritten');
assert(strcmp(before, fileread(coNew)), 'File changed after refused overwrite.')
nPass = pass(nPass, 'new file refuses to overwrite an existing file');

%% 4. Revision --------------------------------------------------------------
v0 = v0of(coNew);
t1 = datenum(2026, 11, 3, 14, 32, 0);

% 4a. Unchanged values: no rewrite
evalc('changed = CharParamsRevise(coNew, v0, t1, ''test'');');
assert(~changed && strcmp(before, fileread(coNew)), 'Unchanged revision rewrote the file.')
nPass = pass(nPass, 'unchanged values do not rewrite the file');

% 4b. One change: old values logged in column 7, new values in column 3
v1 = v0; v1(9) = v0(9) + 5;                   % yrInterpolate
assert(CharParamsRevise(coNew, v1, t1, 'test'), 'Revision reported no change.')
tab = CharParamsReadTable(coNew);
[~, P] = CharParameters(coNew);
assert(size(tab, 2) == 7, 'Expected 7 columns after first revision.')
assert(strcmp(tab{1, 7}, 'Replaced 2026-11-03 14:32, CharAnalysis test'), ...
    'Log header not as expected: %s', tab{1, 7})
assert(isequal(str2double(tab(2:26, 7)), v0) && isequal(str2double(tab(2:26, 3)), v1), ...
    'Logged or new values not as expected.')
assert(P.yrInterp == v1(9), 'CharParameters did not read the revised value.')
nPass = pass(nPass, 'revision logs previous values with time stamp; new values read correctly');

% 4c. Hand-entered log column is preserved and shifts right
hand = repmat({''}, size(tab, 1), 1);
hand{1} = 'my notes, run A'; hand{2} = '-51'; hand{10} = '15';
CharParamsWriteTable([tab, hand], coNew);
v2 = v1; v2(12) = v1(12) + 100;               % smoothing window
CharParamsRevise(coNew, v2, t1 + 1 / 1440, 'test');
tab = CharParamsReadTable(coNew);
assert(size(tab, 2) == 9, 'Expected 9 columns after second revision.')
assert(isequal(str2double(tab(2:26, 7)), v1) && isequal(str2double(tab(2:26, 8)), v0), ...
    'Log columns not ordered newest first.')
assert(strcmp(tab{1, 9}, 'my notes, run A') && strcmp(tab{2, 9}, '-51') && ...
    strcmp(tab{10, 9}, '15'), 'Hand-entered log column not preserved.')
assert(isequal(v0of(coNew), v2), 'CharParameters did not read the revised values.')
nPass = pass(nPass, 'hand-entered log column preserved; logs ordered newest first');

% 4d. Limit of 10 log columns (16 columns in total)
v = v2;
while size(CharParamsReadTable(coNew), 2) < 16
    v(12) = v(12) + 10;
    CharParamsRevise(coNew, v, [], 'test');
end
assert(size(CharParamsReadTable(coNew), 2) == 16, 'Expected 16 columns.')
assert(isequal(v0of(coNew), v), 'CharParameters misread the 16-column file.')
before = fileread(coNew);
vx = v; vx(12) = vx(12) + 10;
expectError(@() CharParamsRevise(coNew, vx, [], 'test'), 'maximum (10)');
assert(strcmp(before, fileread(coNew)), 'File changed after refused revision.')
o3 = evalcRun(coNew);
close all
assert(size(o3.Post.charResults, 1) > 0, '16-column file did not run.')
nPass = pass(nPass, '10-log limit enforced; a full 16-column file is read correctly and runs');

%% 5. Naming ----------------------------------------------------------------
expectError(@() CharParamsNew(v0, fullfile(work, 'CO_params.csv'), template), ...
    '_charParams.csv');
nPass = pass(nPass, 'files not named <site>_charParams.csv are rejected');

%% 6. Output folder and working folder --------------------------------------
dS = fullfile(work, 'CO_save'); mkdir(dS);
copyfile(fullfile(matDir, 'CO_charData.csv'), dS);
vs = v0; vs(24) = 1; vs(23) = 0;              % saveData = 1, saveFigures = 0
CharParamsNew(vs, fullfile(dS, 'CO_charParams.csv'), template);
here = pwd;
evalcRun(fullfile(dS, 'CO_charParams.csv'));
close all
assert(strcmp(pwd, here), 'CharAnalysis did not return to the original folder.')
assert(exist(fullfile(dS, 'CO_charResults.csv'), 'file') == 2, ...
    'Results were not written to the folder holding the parameter file.')
nPass = pass(nPass, 'results written next to the parameter file; working folder restored');

fprintf('\nAll %d parameter-writer tests passed.\n', nPass)

end

%% -- Local functions --------------------------------------------------------

function v = paramsVector(P, S, A, R)
% Inverse of the unpacking in CharParameters (zoneDiv padded with -9999).
zd = P.zoneDiv(:);
v  = [zd; -9999 * ones(8 - numel(zd), 1); P.yrInterp; P.transform; ...
      S.method; S.yr; A.cPeak; A.threshType; A.threshMethod; ...
      A.threshValues(:); A.minCountP; A.peakFrequ; A.bkgSens; ...
      R.saveFigures; R.save; R.allFigures];
end

function v = v0of(f)
[~, P, S, A, R] = CharParameters(f);
v = paramsVector(P, S, A, R);
end

function out = evalcRun(f)
% Run CharAnalysis quietly (command-window output suppressed).
evalc('out = CharAnalysis(f);');
end

function expectError(fn, pattern)
try
    fn();
catch ME
    assert(~isempty(strfind(ME.message, pattern)), ...      %#ok<STREMP>
        'Error message did not contain ''%s'':\n%s', pattern, ME.message)
    return
end
error('Expected an error containing ''%s''; none raised.', pattern)
end

function n = pass(n, label)
n = n + 1;
fprintf('PASS: %s\n', label)
end

function restoreState(d, w)
cd(d);
warning(w);
end
