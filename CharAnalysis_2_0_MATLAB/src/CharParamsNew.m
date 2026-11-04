function CharParamsNew(values, fileName, templateFile)
% CharParamsNew   Write a new *_charParams.csv file from the template.
%   CharParamsNew(values, fileName)
%   CharParamsNew(values, fileName, templateFile)
%
%   values       : numeric vector of the 25 parameter values, in file row
%                  order (zoneDiv padded to 8 values with -9999).
%   fileName     : destination, named <site>_charParams.csv.
%   templateFile : optional; default templates/template_charParams.csv in
%                  the CharAnalysis_2_0_MATLAB folder.
%
%   The file reproduces the template layout exactly (header row, 25
%   parameter rows, and the template's Stage, Variable, Units, and
%   Description text), because CharParameters reads parameters by row
%   position. Parameter Log columns start empty. An existing file is never
%   overwritten; use CharParamsRevise instead.
%
%   Equivalent to the R function char_params_new().
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

nParams = 25;
logCol  = 7;

charParamsCheck(values, fileName)
if exist(fileName, 'file') == 2
    error(['CharParamsNew: a parameter file already exists and will not be ' ...
           'overwritten: %s\nRevise it instead (CharParamsRevise).'], fileName)
end
if nargin < 3 || isempty(templateFile)
    templateFile = charParamsTemplatePath();
end

tmpl = CharParamsReadTable(templateFile);
tab  = repmat({''}, nParams + 1, logCol);
tab(:, 1:5) = tmpl(1:nParams + 1, 1:5);
% Keep the template's "Parameter Log" label in the header of column 7;
% log columns themselves start empty.
if size(tmpl, 2) >= logCol
    tab{1, logCol} = tmpl{1, logCol};
end
tab(2:end, 3) = charParamsFmtNum(values);

CharParamsWriteTable(tab, fileName)

end
