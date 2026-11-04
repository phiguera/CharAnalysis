function [CharAnalysisResults] = CharAnalysis(fileName, runMode, figSelection, saveFigs)
% CharAnalysis   Analyse charcoal time series based on selected parameters.
%
%    *****************************************************************
%                          CharAnalysis 2.0
%                     (c) 2004 - 2026, P.E. Higuera
%             MATLAB (r). (c) 1984 - 2024 The MathWorks, Inc.
%                       philip.higuera@umontana.edu
%                      Please read documentation at:
%                 https://phiguera.github.io/CharAnalysis/
%    *****************************************************************
%
%   CharAnalysis requires an input file in .xls or .csv format. This file
%   includes the selected parameters for peak analysis, and it either
%   includes (.xls) or references (.csv) the input charcoal dataset. See
%   the template file templateChar.xls for details.
%
%   v2.0 changes vs v1.1
%     - global plotData  removed; passed as explicit argument to
%       CharThreshGlobal (and set to 1 here for the main run).
%     - global bkgSensIn removed; passed as explicit argument to
%       CharThreshGlobal (set to 0 here; bkgCharSensitivity sets it to 1).
%     - CharPeakAnalysisResults replaced by CharPostProcess + CharPlotResults.
%     - Version banner updated to 2.0.
%
%   Parameter interface (November 2026; AI-assisted, Claude/Anthropic)
%     - Called with no input, CharAnalysis offers three options: run an
%       existing parameter file (file picker), revise an existing
%       parameter file, or create a new one from a <site>_charData.csv
%       file. The last two open CharParamsGUI and run the analysis after
%       the file is saved. Previous values are kept in time-stamped
%       Parameter Log columns (CharParamsRevise).
%     - When fileName includes a folder, the analysis runs in that folder
%       and returns to the original folder afterwards, so saved results
%       and figures go to the folder holding the parameter file.

%% Banner
disp(' ')
disp('*****************************************************************')
disp('                      CharAnalysis 2.0                           ')
disp('                  (c) 2004 - 2026, P.E. Higuera                  ')
disp('         MATLAB(r). (c) 1984 - 2024 The MathWorks, Inc.          ')
disp('                  philip.higuera@umontana.edu                    ')
disp('                   Please read documentation:                    ')
disp('             https://phiguera.github.io/CharAnalysis/            ')
disp('*****************************************************************')
disp(' ')

%% Default runMode and figSelection
if nargin < 2
    runMode = 'standard';
end
if nargin < 3
    figSelection = [];
end
if nargin < 4
    saveFigs = [];
end
runMode = lower(strtrim(runMode));
if ~ismember(runMode, {'standard', 'modular'})
    error('CharAnalysis: runMode must be ''standard'' or ''modular''.')
end

% Add src/ subfolder to path so all program functions are accessible.
baseDir = fileparts(mfilename('fullpath'));
addpath(fullfile(baseDir, 'src'));

%% No file supplied: run, revise, or create a parameter file
if nargin == 0
    disp('  CharAnalysis requires a parameter file (<site>_charParams.csv  ')
    disp('   or .xls). A .csv parameter file references the charcoal data  ')
    disp('   file <site>_charData.csv, which must be in the same folder.   ')
    disp(' ')
    disp('  If you choose to save the results (figures and/or data), they ')
    disp('   will be saved in the folder containing the parameter file.    ')
    disp(' ')
    disp('*****************************************************************')
    disp(' ')
    disp('  No parameter file specified. What would you like to do?')
    disp('    1: Run an existing parameter file (<site>_charParams.csv or .xls)')
    disp('    2: Revise an existing parameter file in the parameter interface')
    disp('    3: Create a new parameter file from a data file (<site>_charData.csv)')
    choice = str2double(input('  Selection (0 to cancel): ', 's'));
    if ~ismember(choice, [1 2 3])
        disp('CharAnalysis cancelled.')
        CharAnalysisResults = [];
        return
    end
    if choice == 1
        if usejava('desktop')
            [f, p] = uigetfile({'*_charParams.csv;*.xls;*.xlsx', ...
                'CharAnalysis parameter files'}, 'Select the parameter file');
            if isequal(f, 0)
                disp('No file selected. CharAnalysis cancelled.')
                CharAnalysisResults = [];
                return
            end
            fileName = fullfile(p, f);
        else
            fileName = input(['\n Input the file name OR the full path to the ' ...
                'site directory,\n including the file extension (if file name): '], 's');
        end
    else
        % Parameter interface (src/CharParamsGUI.m); returns the saved path.
        modes = {'revise', 'new'};
        fileName = CharParamsGUI('', modes{choice - 1}, 'Save and run');
        if isempty(fileName)
            disp('CharAnalysis cancelled.')
            CharAnalysisResults = [];
            return
        end
    end
end

warning('off', 'all')
clear charData

%% Strip any surrounding single quotes entered interactively
if length(fileName) > 1 && fileName(1) == '''' && fileName(end) == ''''
    fileName = fileName(2:end-1);
end

%% Handle directory-name input
if (~min(fileName(end-2:end) == 'xls') && ...
    ~min(fileName(end-2:end) == 'csv'))
    cd(fileName)
    i1   = regexp(fileName, '\');
    i2   = regexp(fileName, '/');
    i    = max([i1, i2]);
    site = fileName((i(end)+1):end);
    fileName = [site '_CharParams.csv'];
end

%% Work in the folder holding the parameter file
% Results and figures are written to the current folder, so run from the
% folder holding the parameter file and return to the original folder
% afterwards (including on error).
[fDir, fBase, fExt] = fileparts(fileName);
if ~isempty(fDir)
    origDir = pwd;
    cd(fDir)
    restoreDir = onCleanup(@() cd(origDir));   %#ok<NASGU>
    fileName = [fBase fExt];
end

%% (1) Read input file
disp('(1) Reading input file...')
[charData, Pretreatment, Smoothing, PeakAnalysis, Results, site] = ...
    CharParameters(fileName);
disp('      ...done.')

%% (1b) Validate parameters
disp('(1b) Validating input parameters...')
CharValidateParams(charData, Pretreatment, Smoothing, PeakAnalysis, Results);

%% Run-wide flags (set once, used by multiple downstream functions)
plotData  = 1;   % Draw diagnostic figures (Figures 1 and 2) on main run.
bkgSensIn = 0;   % Normal run - not inside the sensitivity loop.

%% (2) Pretreatment
disp('(2) Pretreating charcoal data...')
[Charcoal, Pretreatment, gapIn] = CharPretreatment(charData, site, ...
    Pretreatment, Results, plotData);
disp('      ...done.')

%% (3) Smooth to estimate low-frequency trends
disp('(3) Smoothing resampled CHAR to estimate low-frequency trends')
disp('    and calculating peak CHAR...')
[Charcoal] = CharSmooth(Charcoal, Pretreatment, Smoothing, Results, plotData);

if min(Charcoal.accIS) == 0 && PeakAnalysis.cPeak == 2
    error('Cannot calculate C_peak when C_background = 0; change parameters.')
end

%% (4) Calculate peak CHAR
if PeakAnalysis.cPeak == 1
    Charcoal.peak = Charcoal.accI - Charcoal.accIS;   % Residuals.
else
    Charcoal.peak = Charcoal.accI ./ Charcoal.accIS;  % Ratios.
end
disp('      ...done.')

%% (4) Define thresholds
disp('(4) Defining possible thresholds for peak identification...')

if PeakAnalysis.threshType == 1
    [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment, ...
        PeakAnalysis, site, Results, plotData, bkgSensIn);
else
    [CharThresh] = CharThreshLocal(Charcoal, Smoothing, ...
        PeakAnalysis, site, Results, plotData);
end
disp('      ...done.')

%% (5) Identify peaks
disp('(5) Identifying peaks based on possible thresholds...')
[Charcoal, CharThresh] = CharPeakID(Charcoal, Pretreatment, ...
    PeakAnalysis, CharThresh);
disp('      ...done.')

%% (6) Plot results and save data
if Results.save == 1
    disp('(6) Plotting results and saving data...')
else
    disp('(6) Plotting results...')
end

[Charcoal, Post] = CharPostProcess(Charcoal, Pretreatment, ...
    PeakAnalysis, CharThresh, Smoothing, Results, fileName, gapIn, site);

%% Assemble results struct (used by plotting functions and returned to workspace)
CharAnalysisResults.Charcoal     = Charcoal;
CharAnalysisResults.CharThresh   = CharThresh;
CharAnalysisResults.Pretreatment = Pretreatment;
CharAnalysisResults.Smoothing    = Smoothing;
CharAnalysisResults.PeakAnalysis = PeakAnalysis;
CharAnalysisResults.Results      = Results;
CharAnalysisResults.Post         = Post;
CharAnalysisResults.gapIn        = gapIn;
CharAnalysisResults.site         = site;

if strcmp(runMode, 'modular')
    CharFigureMenu(CharAnalysisResults, figSelection, saveFigs);
else
    CharPlotResults(CharAnalysisResults);
end

disp('      ...CharAnalysis finished.')

if Results.save == 1
    disp(' ')
    if min(fileName(end-2:end) == 'xls') > 0
        disp('   Results saved to input file:')
        disp(['     ', cd, '\', fileName])
    else
        disp('   Results saved to file:')
        disp(['     ', cd, '\', site, '_CharResults.csv'])
    end
    disp(' ')
end

%% (7) Background sensitivity analysis (optional)
if PeakAnalysis.bkgSens == 1
    disp('(7) Running C_background sensitivity analysis:')
    [z, GOF_i, SNI_i] = bkgCharSensitivity(Charcoal, CharThresh, ...
        PeakAnalysis, Pretreatment, Smoothing, Results, site);
    disp('C_background sensitivity analysis finished.')
end

%% Bring all figures to the front in order (standard path only)
if strcmp(runMode, 'standard')
    if PeakAnalysis.bkgSens == 1
        figIn = 1:10;
    else
        figIn = 1:9;
    end
    if Results.allFigures ~= 1
        figIn(ismember(figIn, [1, 2, 9])) = [];
    end
    for i = 1:length(figIn)
        figure(figIn(i))
    end
end