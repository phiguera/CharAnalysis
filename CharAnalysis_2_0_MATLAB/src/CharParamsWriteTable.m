function CharParamsWriteTable(tab, fileName)
% CharParamsWriteTable   Write a cell table of text to a parameter file.
%   CharParamsWriteTable(tab, fileName)
%
%   Writes tab (cell array of char) as comma-separated lines with LF line
%   endings. Fields are quoted only when they contain a comma, a quote, or
%   a line break, so column 3 (numeric values) is never quoted and both the
%   MATLAB reader (textscan in CharParameters) and the R reader parse it.
%
%   The file is first written to a temporary file in the same folder, read
%   back and verified, and then moved over the destination, so an
%   interrupted save cannot leave a partial parameter file. Content must be
%   plain ASCII.
%
%   Internal helper for CharParamsNew and CharParamsRevise; equivalent to
%   the R function .cp_write_atomic().
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

% Build lines
nRow  = size(tab, 1);
lines = cell(nRow, 1);
for r = 1:nRow
    f = tab(r, :);
    for c = 1:numel(f)
        s = f{c};
        if isempty(s)
            s = '';
        end
        s = char(s);
        if any(s > 127)
            error(['CharParamsWriteTable: parameter file content must be ' ...
                   'plain ASCII (row %d, column %d).'], r, c)
        end
        if any(s == ',' | s == '"' | s == sprintf('\n') | s == sprintf('\r'))
            s = ['"' strrep(s, '"', '""') '"'];
        end
        f{c} = s;
    end
    lines{r} = strjoin(f, ',');
end
content = [strjoin(lines', sprintf('\n')) sprintf('\n')];

% Write to a temporary file in the destination folder
outDir = fileparts(fileName);
if isempty(outDir)
    outDir = pwd;
end
[~, tmpName] = fileparts(tempname);
tmpFile = fullfile(outDir, ['.charParams_' tmpName '.tmp']);

fid = fopen(tmpFile, 'w');           % binary mode: LF is written as LF
if fid < 0
    error('CharParamsWriteTable: could not write to folder: %s', outDir)
end
fwrite(fid, content, 'char');
fclose(fid);

% Verify, then swap in
check = fileread(tmpFile);
if ~strcmp(check, content)
    delete(tmpFile)
    error('CharParamsWriteTable: verification of the temporary parameter file failed; nothing was written.')
end
[ok, msg] = movefile(tmpFile, fileName, 'f');
if ~ok
    if exist(tmpFile, 'file') == 2
        delete(tmpFile)
    end
    error('CharParamsWriteTable: could not write parameter file %s (%s)', ...
          fileName, msg)
end

end
