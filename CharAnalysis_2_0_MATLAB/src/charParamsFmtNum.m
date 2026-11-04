function s = charParamsFmtNum(x)
% charParamsFmtNum   Format parameter values as plain text for a CSV file.
%   s = charParamsFmtNum(x)
%
%   Returns a cell array (same number of elements as x) of char strings,
%   using up to 15 significant digits and no trailing zeros (e.g., 0.95,
%   500, -9999), readable by textscan %n in CharParameters and by R.
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

s = cell(numel(x), 1);
for i = 1:numel(x)
    if isnan(x(i))
        s{i} = '';
    else
        s{i} = sprintf('%.15g', x(i));
    end
end

end
