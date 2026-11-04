function tab = CharParamsReadTable(fileName)
% CharParamsReadTable   Read a *_charParams.csv file as a cell table of text.
%   tab = CharParamsReadTable(fileName)
%
%   Returns every row and column of the parameter file as a cell array of
%   char row vectors (empty cells are ''). Row 1 is the header; rows 2-26
%   hold the 25 parameters. Rows are padded to the widest row. Quoted
%   fields (which may contain commas or doubled quotes) are unquoted.
%
%   Used by the parameter interface (CharParamsGUI) and the parameter-file
%   writers (CharParamsNew, CharParamsRevise). Equivalent to the R function
%   char_params_read_table().
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

nParams = 25;

if exist(fileName, 'file') ~= 2
    error('CharParamsReadTable: parameter file not found: %s', fileName)
end

txt = fileread(fileName);
txt = strrep(txt, sprintf('\r\n'), sprintf('\n'));   % CRLF -> LF
txt = strrep(txt, sprintf('\r'), sprintf('\n'));     % lone CR -> LF

rows   = {};          % cell array of rows; each row a cell array of fields
fields = {};
buf    = blanks(0);
inQ    = false;
i      = 1;
n      = length(txt);
while i <= n
    c = txt(i);
    if inQ
        if c == '"'
            if i < n && txt(i+1) == '"'
                buf(end+1) = '"';                  %#ok<AGROW>
                i = i + 1;
            else
                inQ = false;
            end
        else
            buf(end+1) = c;                        %#ok<AGROW>
        end
    else
        if c == '"'
            inQ = true;
        elseif c == ','
            fields{end+1} = buf;                   %#ok<AGROW>
            buf = blanks(0);
        elseif c == sprintf('\n')
            fields{end+1} = buf;                   %#ok<AGROW>
            rows{end+1}   = fields;                %#ok<AGROW>
            fields = {};
            buf    = blanks(0);
        else
            buf(end+1) = c;                        %#ok<AGROW>
        end
    end
    i = i + 1;
end
if ~isempty(buf) || ~isempty(fields)
    fields{end+1} = buf;
    rows{end+1}   = fields;
end

% Drop blank lines (a row with one empty field)
keep = true(1, numel(rows));
for r = 1:numel(rows)
    keep(r) = ~(numel(rows{r}) == 1 && isempty(strtrim(rows{r}{1})));
end
rows = rows(keep);

if numel(rows) < nParams + 1
    error(['CharParamsReadTable: parameter file has %d rows; expected a ' ...
           'header plus %d parameter rows: %s'], numel(rows), nParams, fileName)
end

nCol = max(5, max(cellfun(@numel, rows)));
tab  = repmat({''}, numel(rows), nCol);
for r = 1:numel(rows)
    tab(r, 1:numel(rows{r})) = rows{r};
end

end
