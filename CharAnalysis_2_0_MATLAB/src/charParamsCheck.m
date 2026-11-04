function charParamsCheck(values, fileName)
% charParamsCheck   Check inputs to CharParamsNew and CharParamsRevise.
%   charParamsCheck(values, fileName)
%
%   values must be 25 finite numbers; the file must be named
%   <site>_charParams.csv.
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

[~, name, ext] = fileparts(fileName);
if isempty(regexpi([name ext], '_charParams\.csv$', 'once'))
    error('Parameter file must be named <site>_charParams.csv: %s', [name ext])
end
if ~isnumeric(values) || numel(values) ~= 25
    error('Parameter values must be a numeric vector of length 25.')
end
bad = find(~isfinite(values));
if ~isempty(bad)
    error('Parameter values contain NaN or non-finite entries at position(s): %s.', ...
          num2str(bad(:)'))
end

end
