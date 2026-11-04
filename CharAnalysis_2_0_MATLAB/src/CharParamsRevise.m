function changed = CharParamsRevise(fileName, values, timeStamp, versionStr)
% CharParamsRevise   Revise a *_charParams.csv file, logging previous values.
%   changed = CharParamsRevise(fileName, values)
%   changed = CharParamsRevise(fileName, values, timeStamp, versionStr)
%
%   fileName   : existing <site>_charParams.csv file.
%   values     : numeric vector of the 25 new parameter values.
%   timeStamp  : optional time for the log header (datenum; default now).
%   versionStr : optional version recorded in the log header
%                (default '2.0 MATLAB').
%   changed    : true if the file was revised; false if no value changed
%                (the file is then not rewritten).
%
%   If any value differs, the previous column-3 values move to a new
%   Parameter Log column (column 7), headed with a time stamp and version,
%   e.g. "Replaced 2026-11-04 14:32, CharAnalysis 2.0 MATLAB", and the new
%   values are written to column 3. Existing log columns, including
%   hand-entered ones, shift right unchanged. Log columns with no content in
%   the parameter rows (e.g., the template's empty, labelled column) are
%   dropped. Rows below the parameter block, if any, are kept.
%
%   CharParameters reads a fixed 16 fields per row (textscan), so a file
%   may hold at most 10 log columns (columns 7-16). When a file already
%   holds 10, the revision stops with an error and the file is unchanged.
%
%   Equivalent to the R function char_params_revise().
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

nParams = 25;
logCol  = 7;
maxCol  = 16;
maxLog  = maxCol - logCol + 1;     % 10 log columns

if nargin < 3 || isempty(timeStamp)
    timeStamp = now;
end
if nargin < 4 || isempty(versionStr)
    versionStr = '2.0 MATLAB';
end

charParamsCheck(values, fileName)
values = values(:);
tab    = CharParamsReadTable(fileName);
pRow   = 2:(nParams + 1);

oldStr = tab(pRow, 3);
old    = str2double(oldStr);
if ~any(isnan(old)) && all(old == values)
    [~, n, e] = fileparts(fileName);
    fprintf('No parameter values changed; %s was not rewritten.\n', [n e])
    changed = false;
    return
end

if size(tab, 2) < 6
    tab = [tab, repmat({''}, size(tab, 1), 6 - size(tab, 2))];
end
keep = [];
for j = logCol:size(tab, 2)
    if any(~cellfun(@isempty, strtrim(tab(pRow, j))))
        keep(end+1) = j;                                   %#ok<AGROW>
    end
end
if numel(keep) >= maxLog
    [~, n, e] = fileparts(fileName);
    error(['%s already holds %d logged parameter sets, the maximum (%d) ' ...
           'that keeps the file readable by CharParameters.\nSave a copy of ' ...
           'the file to keep its history, then delete older log columns ' ...
           '(rightmost) before revising.'], [n e], numel(keep), maxLog)
end

newCol      = repmat({''}, size(tab, 1), 1);
newCol{1}   = sprintf('Replaced %s, CharAnalysis %s', ...
                      datestr(timeStamp, 'yyyy-mm-dd HH:MM'), versionStr);
newCol(pRow) = oldStr;

out = [tab(:, 1:6), newCol, tab(:, keep)];
out(pRow, 3) = charParamsFmtNum(values);
CharParamsWriteTable(out, fileName)
changed = true;

end
