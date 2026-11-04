function f = charParamsTemplatePath()
% charParamsTemplatePath   Path to the bundled parameter-file template.
%   f = charParamsTemplatePath()
%
%   Returns <CharAnalysis_2_0_MATLAB>/templates/template_charParams.csv
%   (a byte copy of DataTemplates_and_Examples/template_charParams.csv).
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

srcDir = fileparts(mfilename('fullpath'));
f = fullfile(fileparts(srcDir), 'templates', 'template_charParams.csv');
if exist(f, 'file') ~= 2
    error('charParamsTemplatePath: template not found: %s', f)
end

end
