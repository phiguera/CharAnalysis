function paramsFile = CharParamsGUI(fileName, mode, buttonLabel)
% CharParamsGUI   Create or revise a parameter file with a graphical interface.
%   paramsFile = CharParamsGUI
%   paramsFile = CharParamsGUI(fileName)
%   paramsFile = CharParamsGUI(fileName, mode, buttonLabel)
%
%   Opens a window containing every parameter in a <site>_charParams.csv
%   file, with a summary of the charcoal record, plain-language definitions,
%   and advisory notes. On saving, the parameter file is written next to the
%   data file. When revising an existing file, the previous values are kept
%   in the file's Parameter Log columns (see CharParamsRevise).
%
%   INPUTS (all optional)
%     fileName    : a data file named <site>_charData.csv (create a new
%                   parameter file) or an existing <site>_charParams.csv
%                   (revise it). If empty, you are asked which to do and a
%                   file picker opens.
%     mode        : 'new' or 'revise'; used only when fileName is empty.
%     buttonLabel : text of the save button (default 'Save').
%
%   OUTPUT
%     paramsFile  : full path to the saved parameter file, or '' if the
%                   interface was cancelled. Run the analysis with
%                   CharAnalysis(paramsFile).
%
%   Data files must follow the <site>_charData.csv naming convention,
%   because CharParameters finds the data file from the parameter file
%   name. Errors (values CharValidateParams rejects) block saving; notes
%   advise. A file can hold up to 10 logged parameter sets.
%
%   Requires MATLAB R2018b or later (uifigure, uigridlayout). Uses no
%   features newer than R2019a, the documented minimum for CharAnalysis
%   2.0 (so no label word wrap, HTML labels, or tooltips). Equivalent to the R function char_params_gui(); the R
%   interface hides allFigures (no effect in R), whereas this interface
%   shows it because it controls Figures 1, 2, and 9 in MATLAB.
%
%   CharAnalysis v2.0 parameter interface.
%   AI-assisted development (Claude/Anthropic).

paramsFile = '';
if nargin < 1, fileName    = ''; end
if nargin < 2, mode        = ''; end
if nargin < 3 || isempty(buttonLabel), buttonLabel = 'Save'; end

%% -- Choose a file ---------------------------------------------------------
if isempty(fileName)
    if isempty(mode)
        optNew = 'Create a new parameter file';
        optRev = 'Revise an existing parameter file';
        ans1 = questdlg('What would you like to do?', ...
            'CharAnalysis parameters', optNew, optRev, 'Cancel', optNew);
        if strcmp(ans1, optNew)
            mode = 'new';
        elseif strcmp(ans1, optRev)
            mode = 'revise';
        else
            return
        end
    end
    if strcmp(mode, 'new')
        [f, p] = uigetfile('*_charData.csv', ...
            'Select the data file (<site>_charData.csv)');
    else
        [f, p] = uigetfile('*_charParams.csv', ...
            'Select the parameter file (<site>_charParams.csv)');
    end
    if isequal(f, 0)
        disp('No file selected; parameter interface cancelled.')
        return
    end
    fileName = fullfile(p, f);
end

fileName = charAbsPath(fileName);
if exist(fileName, 'file') ~= 2
    error('CharParamsGUI: file not found: %s', fileName)
end
[fDir, fName, fExt] = fileparts(fileName);
base = [fName fExt];

if ~isempty(regexpi(base, '_charData\.csv$', 'once'))
    mode       = 'new';
    site       = regexprep(base, '_charData\.csv$', '', 'ignorecase');
    dataFile   = fileName;
    paramsPath = fullfile(fDir, [site '_charParams.csv']);
    if exist(paramsPath, 'file') == 2
        ans2 = questdlg(sprintf(['A parameter file already exists for this ' ...
            'site:\n  %s\nRevise it instead?'], paramsPath), ...
            'CharAnalysis parameters', 'Revise', 'Cancel', 'Revise');
        if ~strcmp(ans2, 'Revise')
            return
        end
        mode = 'revise';
    end
elseif ~isempty(regexpi(base, '_charParams\.csv$', 'once'))
    mode       = 'revise';
    site       = regexprep(base, '_charParams\.csv$', '', 'ignorecase');
    paramsPath = fileName;
    dataFile   = fullfile(fDir, [site '_charData.csv']);
    if exist(dataFile, 'file') ~= 2
        error(['CharParamsGUI: the companion data file was not found:\n  %s\n' ...
               'CharAnalysis requires the data file and parameter file to sit ' ...
               'in the same folder, named <site>_charData.csv and ' ...
               '<site>_charParams.csv.'], dataFile)
    end
else
    error(['CharParamsGUI: ''%s'' does not follow the CharAnalysis naming ' ...
           'convention.\nData files must be named <site>_charData.csv (e.g., ' ...
           'CO_charData.csv), and parameter files <site>_charParams.csv. ' ...
           'CharAnalysis finds the data file from the parameter file name, so ' ...
           'this convention is required. Rename the file and try again.'], base)
end

%% -- Data and initial values -----------------------------------------------
charData = csvread(dataFile, 1, 0);
summ     = charDataSummary(charData);

if strcmp(mode, 'new')
    tmpl   = CharParamsReadTable(charParamsTemplatePath());
    init   = str2double(tmpl(2:26, 3));
    init(1:8) = [summ.ageTop; summ.ageBot; -9999 * ones(6, 1)];  % start/end only
    init(9)   = summ.resSuggest;
    old    = [];
    nLogs  = 0;
else
    tab    = CharParamsReadTable(paramsPath);
    old    = str2double(tab(2:26, 3));
    init   = old;
    nLogs  = 0;
    for j = 7:size(tab, 2)
        nLogs = nLogs + any(~cellfun(@isempty, strtrim(tab(2:26, j))));
    end
end
init = init(:);
if isnan(init(25))
    init(25) = 1;      % allFigures default
end

%% -- Labels and keys (one per parameter position) --------------------------
posKey = [repmat({'zoneDiv'}, 1, 8), {'yrInterp', 'transform', 'method', ...
    'smoYr', 'cPeak', 'threshType', 'threshMethod'}, ...
    repmat({'threshValues'}, 1, 4), {'minCountP', 'peakFrequ', 'bkgSens', ...
    'saveFigures', 'saveData', 'allFigures'}];
posLabel = [{'Zone start'}, arrayfun(@(k) sprintf('Zone value %d', k), 2:8, ...
    'UniformOutput', false), {'yrInterpolate', 'transform', ...
    'Smoothing method', 'Smoothing window (yr)', 'cPeak', 'threshType', ...
    'threshMethod', 'Peaks 1 threshold', 'Peaks 2 threshold', ...
    'Peaks 3 threshold', 'Peaks Final threshold', 'minCountP', 'peakFrequ', ...
    'C_background sensitivity', 'saveFigures', 'saveData', 'allFigures'}];
keyTab = struct('data', 'Data', 'zoneDiv', 'Pretreatment', ...
    'yrInterp', 'Pretreatment', 'transform', 'Pretreatment', ...
    'method', 'Smoothing', 'smoYr', 'Smoothing', 'cPeak', 'Peaks', ...
    'threshType', 'Peaks', 'threshMethod', 'Peaks', 'threshValues', 'Peaks', ...
    'minCountP', 'Peaks', 'peakFrequ', 'Results', 'bkgSens', 'Results', ...
    'saveFigures', 'Results', 'saveData', 'Results', 'allFigures', 'Results');

%% -- Build the window ------------------------------------------------------
if strcmp(mode, 'new')
    titleStr = sprintf('CharAnalysis parameters: %s (new file)', site);
else
    titleStr = sprintf('CharAnalysis parameters: %s (revising)', site);
end
scr = get(0, 'ScreenSize');
figW = 1060; figH = min(840, scr(4) - 80);
fig = uifigure('Name', titleStr, ...
    'Position', [max(20, (scr(3) - figW) / 2), max(40, (scr(4) - figH) / 2), figW, figH]);
fig.CloseRequestFcn = @(~, ~) onCancel();

colErr  = [1.00 0.88 0.88];
colNote = [1.00 0.95 0.84];
wrapDef  = 34;     % characters per line, column 1
wrapDesc = 66;     % characters per line, column 3
lineH    = 16;

gMain = uigridlayout(fig, [3 1]);
gMain.RowHeight   = {24, '1x', 36};
gMain.ColumnWidth = {'1x'};
lblTitle = uilabel(gMain, 'Text', titleStr, 'FontWeight', 'bold', 'FontSize', 13);
lblTitle.Layout.Row = 1; lblTitle.Layout.Column = 1;

tg = uitabgroup(gMain);
tg.Layout.Row = 2; tg.Layout.Column = 1;

gBtn = uigridlayout(gMain, [1 3]);
gBtn.Layout.Row = 3; gBtn.Layout.Column = 1;
gBtn.ColumnWidth = {'1x', 100, 140};
gBtn.RowHeight   = {'1x'};
gBtn.Padding     = [0 0 0 0];
btnCancel = uibutton(gBtn, 'Text', 'Cancel', 'ButtonPushedFcn', @(~, ~) onCancel());
btnCancel.Layout.Row = 1; btnCancel.Layout.Column = 2;
btnSave = uibutton(gBtn, 'Text', buttonLabel, 'FontWeight', 'bold', ...
    'ButtonPushedFcn', @(~, ~) onSave());
btnSave.Layout.Row = 1; btnSave.Layout.Column = 3;

ctl   = struct();     % input controls, by id
ctlOf = struct();     % key -> cell array of controls to highlight
notes = struct();     % tab name -> notes text area
tabs  = struct();     % tab name -> uitab

% ---- Data tab --------------------------------------------------------------
tabData = uitab(tg, 'Title', 'Data');
gData = uigridlayout(tabData, [4 1]);
gData.RowHeight = {230, 22, 22, '1x'};
summaryCell = {
    'Data file',                         [site '_charData.csv'];
    'Samples',                           sprintf('%d', summ.n);
    'Depth range (cm)',                  sprintf('%g to %g', summ.cmTop, summ.cmBot);
    'Age range (cal. yr BP)',            sprintf('%g to %g', summ.ageTop, summ.ageBot);
    'Record length (yr)',                sprintf('%g', summ.ageBot - summ.ageTop);
    'Sample resolution (yr per sample)', sprintf('median %.1f (range %.1f to %.1f)', ...
                                         summ.resMedian, summ.resMin, summ.resMax);
    'Samples with missing values',       sprintf('%d', summ.nMissing);
    'Breaks in sample depths',           sprintf('%d', summ.nGaps)};
tbl = uitable(gData, 'Data', summaryCell, 'ColumnName', {'Measure', 'Value'}, ...
    'RowName', {}, 'ColumnWidth', {300, 400});
tbl.Layout.Row = 1; tbl.Layout.Column = 1;
if strcmp(mode, 'revise')
    s1 = sprintf(['Revising %s_charParams.csv. The file holds %d of a maximum ' ...
        'of 10 logged parameter sets.'], site, nLogs);
else
    s1 = 'Creating a new parameter file from the template.';
end
l1 = uilabel(gData, 'Text', s1); l1.Layout.Row = 2; l1.Layout.Column = 1;
l2 = uilabel(gData, 'Text', sprintf(['Parameters will be saved to ' ...
    '%s_charParams.csv, in the same folder as the data file.'], site));
l2.Layout.Row = 3; l2.Layout.Column = 1;
notes.Data = uitextarea(gData, 'Editable', 'off', 'Value', {''});
notes.Data.Layout.Row = 4; notes.Data.Layout.Column = 1;
ctlOf.data = {};

% ---- Pretreatment tab ------------------------------------------------------
zinit  = init(1:8); zinit = zinit(~isnan(zinit) & zinit ~= -9999);
zStart = NaN; zEnd = NaN; zInt = [];
if numel(zinit) >= 1, zStart = zinit(1); end
if numel(zinit) >= 2, zEnd = zinit(end); zInt = zinit(2:end-1); end
zIntStr = strjoin(charParamsFmtNum(zInt)', ', ');

if strcmp(mode, 'new')
    zDesc1 = ['Start and end are pre-filled from the data (top age of the ' ...
        'youngest sample, bottom age of the oldest).'];
else
    zDesc1 = 'Values from the existing parameter file.';
end
gPre = newTab('Pretreatment', 3);
addRow(gPre, 1, 'zoneDiv', 'Analysis period and zone boundaries (cal. yr BP).', ...
    @makeZoneInputs, ...
    [zDesc1 ' Enter one or more zone boundaries, youngest first, separated ' ...
     'by commas (e.g., 1000, 3500, 5500). Leave blank for a single zone.'], ...
    {'Zones divide the record for plotting and for comparing fire return intervals between periods.', ...
     'Up to six boundaries can be entered; they must fall between the start and end.'}, 168);
addRow(gPre, 2, 'yrInterpolate', 'Time step the record is resampled to (yr).', ...
    @(g) makeNum(g, 'yrInterp', init(9), 'yrInterp'), ...
    sprintf(['Counts, volumes, and depths are interpolated before CHAR is ' ...
     'calculated. Median sample resolution is %.1f yr (suggested: %d). ' ...
     'Enter 0 to use the median automatically.'], summ.resMedian, summ.resSuggest), ...
    {'Resampling to equal time steps puts every sample on the same time scale, which the smoothing and threshold steps require.', ...
     'An interval finer than the sample resolution creates values that are not independent measurements.'}, 0);
addRow(gPre, 3, 'transform', 'Whether CHAR is transformed before analysis.', ...
    @(g) makeDrop(g, 'transform', {'0: No transformation', '1: Base-10 log', ...
        '2: Natural log'}, [0 1 2], init(10), 'transform'), ...
    'For log transformation, 1 is added to all CHAR values first.', ...
    {'Log transformation compresses large CHAR values relative to small ones.'}, 0);

% ---- Smoothing tab ---------------------------------------------------------
gSmo = newTab('Smoothing', 2);
addRow(gSmo, 1, 'method', 'How low-frequency CHAR (C_background) is estimated.', ...
    @(g) makeDrop(g, 'method', {'1: Lowess', '2: Lowess, robust to outliers', ...
        '3: Moving average', '4: Moving median', '5: Moving mode'}, 1:5, ...
        init(11), 'method'), ...
    ['The smoother used to estimate the slowly varying background component ' ...
     'of the record. Click "i" for what each option does.'], ...
    {'Lowess: locally weighted regression within the window.', ...
     'Lowess, robust to outliers: lowess that iteratively down-weights outliers, reducing the influence of peaks on C_background.', ...
     'Moving average: mean of the values in the window.', ...
     'Moving median: running median of the window, then smoothed with lowess.', ...
     'Moving mode: running mode of the window (from a histogram of the values), then smoothed with lowess.'}, 0);
addRow(gSmo, 2, 'yr', 'Window over which C_background is estimated (yr).', ...
    @(g) makeNum(g, 'smoYr', init(12), 'smoYr'), ...
    ['Years to smooth the record over. C_background sensitivity (Results ' ...
     'tab) shows how results change with this window.'], ...
    {'For local thresholds, this window also sets the length of the moving window in which thresholds are defined.'}, 0);

% ---- Peaks tab -------------------------------------------------------------
gPk = newTab('Peaks', 5);
addRow(gPk, 1, 'cPeak', 'How high-frequency CHAR (C_peak) is calculated.', ...
    @(g) makeDrop(g, 'cPeak', {'1: Residuals', '2: Ratios'}, [1 2], init(13), 'cPeak'), ...
    ['Residuals: C_peak = C_interpolated - C_background. ' ...
     'Ratios: C_peak = C_interpolated / C_background.'], ...
    {'Residuals: C_peak is measured in CHAR units above C_background.', ...
     'Ratios: C_peak is measured relative to C_background.'}, 0);
addRow(gPk, 2, 'threshType', ...
    'Whether one threshold or locally varying thresholds separate peaks from noise.', ...
    @(g) makeDrop(g, 'threshType', {'1: Global', '2: Local'}, [1 2], init(14), 'threshType'), ...
    ['Global: one threshold for the whole record. Local: thresholds are ' ...
     'defined within a moving window (equal to the length selected in ' ...
     'Smoothing) and thus vary through time.'], ...
    {'Global: one noise distribution and threshold for the whole record.', ...
     'Local: a noise distribution and threshold are estimated within each window, so thresholds vary through time.'}, 0);
addRow(gPk, 3, 'threshMethod', 'How threshold values are determined.', ...
    @(g) makeDrop(g, 'threshMethod', {'1: User-defined', '2: Gaussian noise', ...
        '3: Gaussian mixture'}, [1 2 3], init(15), 'threshMethod'), ...
    ['For methods 2 and 3, threshold values are percentiles of the noise ' ...
     'distribution (e.g., 0.95). For method 1, values are in C_peak units.'], ...
    {'User-defined: you supply thresholds in C_peak units (global thresholds only).', ...
     'Gaussian noise: noise is modeled as a Gaussian distribution with mean 0 (residuals) or 1 (ratios).', ...
     'Gaussian mixture: C_peak values are modeled as a mixture of two Gaussian populations, noise and signal; thresholds are percentiles of the noise population.'}, 0);
addRow(gPk, 4, 'threshValues', 'Thresholds evaluated; Peaks Final is used for analysis.', ...
    @makeThreshInputs, ...
    ['Peaks 1-3 are for comparison. Peaks Final is used for peak ' ...
     'identification and all calculations in the Results tab.'], ...
    {'Names match the peaks 1, peaks 2, peaks 3, and peaks Final columns in the _charResults.csv output.'}, 100);
addRow(gPk, 5, 'minCountP', 'Cut-off for screening peaks that could arise from counting noise.', ...
    @(g) makeNum(g, 'minCountP', init(20), 'minCountP'), ...
    ['A peak is flagged and excluded if the minimum count in the 75 yr ' ...
     'before it has more than this probability of coming from the same ' ...
     'Poisson distribution as the peak''s maximum count. Set to 0.99 to ' ...
     'turn off.'], ...
    {'Screens out peaks that could result from sampling variability when charcoal counts are low. Flagged peaks are still plotted.'}, 0);

% ---- Results tab -----------------------------------------------------------
yn = {'0: No', '1: Yes'};
gRes = newTab('Results', 5);
addRow(gRes, 1, 'peakFrequ', ...
    'Window for smoothing fire frequency and fire return intervals (yr).', ...
    @(g) makeNum(g, 'peakFrequ', init(21), 'peakFrequ'), ...
    ['Years over which fire frequency and fire return intervals, from ' ...
     'Peaks Final, are smoothed.'], {}, 0);
addRow(gRes, 2, 'C_background sensitivity', ...
    'Whether the analysis is repeated across smoothing windows.', ...
    @(g) makeDrop(g, 'bkgSens', yn, [0 1], init(22), 'bkgSens'), ...
    ['Shows how sensitive results are to the C_background window ' ...
     '(Figure 10). Recommended for the first analysis of a record.'], ...
    {'Repeats the analysis across a range of smoothing windows and plots how the results change.'}, 0);
addRow(gRes, 3, 'saveData', 'Save the results table.', ...
    @(g) makeDrop(g, 'saveData', yn, [0 1], init(24), 'saveData'), ...
    sprintf(['Writes %s_charResults.csv to the folder holding the parameter ' ...
     'file, replacing any earlier results file. Results can always be ' ...
     'recreated from the saved parameters.'], site), {}, 0);
addRow(gRes, 4, 'saveFigures', 'Save the figures.', ...
    @(g) makeDrop(g, 'saveFigures', yn, [0 1], init(23), 'saveFigures'), ...
    ['Writes figures as PDF and TIFF files to the same folder, replacing ' ...
     'earlier copies.'], {}, 0);
addRow(gRes, 5, 'allFigures', 'Show and save the diagnostic figures.', ...
    @(g) makeDrop(g, 'allFigures', yn, [0 1], init(25), 'allFigures'), ...
    ['Yes: all figures are shown and, if saveFigures is Yes, saved. No: ' ...
     'Figure 9 is not drawn, and Figures 1 and 2 are not brought forward ' ...
     'or saved.'], ...
    {'Figures 1 and 2 show pretreatment and threshold determination; Figure 9 shows threshold details. They are useful for checking an analysis.'}, 0);

% ---- Review tab ------------------------------------------------------------
tabRev = uitab(tg, 'Title', 'Review');
gRev = uigridlayout(tabRev, [1 1]);
txtReview = uitextarea(gRev, 'Editable', 'off', 'Value', {''});

tg.SelectedTab = tabs.Pretreatment;
result = '';
refresh();
uiwait(fig);

if isvalid(fig)
    delete(fig)
end
paramsFile = result;
if isempty(paramsFile)
    disp('Parameter interface cancelled; no file written.')
else
    fprintf('Parameter file: %s\n', paramsFile)
    fprintf('To run (or re-run) this analysis from the command window:\n  CharAnalysis(''%s'')\n', ...
        paramsFile)
end

%% ======================================================================
%% Nested functions (share the variables above)
%% ======================================================================

    function g = newTab(name, nRows)
        t = uitab(tg, 'Title', name);
        tabs.(name) = t;
        g = uigridlayout(t, [nRows + 2, 3]);
        g.ColumnWidth = {250, 250, '1x'};
        g.RowHeight   = [{22}, repmat({60}, 1, nRows), {'1x'}];
        g.RowSpacing  = 6;
        hdr = {'Parameter', 'Setting', 'Description'};
        for c = 1:3
            h = uilabel(g, 'Text', hdr{c}, 'FontWeight', 'bold', ...
                'FontColor', [0.33 0.33 0.33]);
            h.Layout.Row = 1; h.Layout.Column = c;
        end
        notes.(name) = uitextarea(g, 'Editable', 'off', 'Value', {''});
        notes.(name).Layout.Row = nRows + 2;
        notes.(name).Layout.Column = [1 3];
    end

    function addRow(g, i, name, def, makeInput, desc, tip, inputH)
        r = i + 1;
        defLines  = wrapText(def, wrapDef);
        descLines = wrapText(desc, wrapDesc);
        h = max([22 + numel(defLines) * lineH, numel(descLines) * lineH + 4, inputH]);
        g.RowHeight{r} = h + 6;

        % Column 1: name, info button, definition
        g1 = uigridlayout(g, [2 2]);
        g1.Layout.Row = r; g1.Layout.Column = 1;
        g1.RowHeight = {18, '1x'}; g1.ColumnWidth = {'1x', 24};
        g1.Padding = [0 0 0 0]; g1.RowSpacing = 2; g1.ColumnSpacing = 2;
        hn = uilabel(g1, 'Text', name, 'FontWeight', 'bold');
        hn.Layout.Row = 1; hn.Layout.Column = 1;
        if ~isempty(tip)
            hb = uibutton(g1, 'Text', 'i', 'FontWeight', 'bold', ...
                'ButtonPushedFcn', @(~, ~) uialert(fig, strjoin(tip, ...
                    [newline newline]), name, 'Icon', 'info'));
            hb.Layout.Row = 1; hb.Layout.Column = 2;
        end
        hd = uilabel(g1, 'Text', defLines, 'VerticalAlignment', 'top', ...
            'FontColor', [0.27 0.27 0.27]);
        hd.Layout.Row = 2; hd.Layout.Column = [1 2];

        % Column 2: input(s)
        g2 = uigridlayout(g, [1 1]);
        g2.Layout.Row = r; g2.Layout.Column = 2;
        g2.Padding = [0 0 0 0];
        makeInput(g2);

        % Column 3: description
        hs = uilabel(g, 'Text', descLines, 'VerticalAlignment', 'top', ...
            'FontColor', [0.33 0.33 0.33]);
        hs.Layout.Row = r; hs.Layout.Column = 3;
    end

    function h = makeNum(g, id, val, key)
        str = charParamsFmtNum(val);
        g.RowHeight = {24, '1x'};
        h = uieditfield(g, 'text', 'Value', str{1}, ...
            'ValueChangedFcn', @(~, ~) refresh());
        h.Layout.Row = 1; h.Layout.Column = 1;
        ctl.(id) = h;
        h.UserData = h.BackgroundColor;
        if isfield(ctlOf, key), ctlOf.(key){end+1} = h; else, ctlOf.(key) = {h}; end
    end

    function h = makeDrop(g, id, items, data, val, key)
        if isnan(val)
            items{end+1} = '(missing in file)';
            data(end+1)  = -99999;
            val = -99999;
        elseif ~ismember(val, data)
            items{end+1} = sprintf('%g: (value in file; not a valid option)', val);
            data(end+1)  = val;
        end
        g.RowHeight = {24, '1x'};
        h = uidropdown(g, 'Items', items, 'ItemsData', data, 'Value', val, ...
            'ValueChangedFcn', @(~, ~) refresh());
        h.Layout.Row = 1; h.Layout.Column = 1;
        ctl.(id) = h;
        h.UserData = h.BackgroundColor;
        if isfield(ctlOf, key), ctlOf.(key){end+1} = h; else, ctlOf.(key) = {h}; end
    end

    function makeZoneInputs(g)
        g.RowHeight = {16, 24, 8, 32, 24, 8, 16, 24, '1x'};
        g.RowSpacing = 1;
        stackLabel(g, 1, 'Start');
        ctl.zStart = numField(g, 2, zStart);
        stackLabel(g, 4, {'Zone boundaries (optional; separate', ...
                          'multiple values with commas)'});
        ctl.zBounds = uieditfield(g, 'text', 'Value', zIntStr, ...
            'ValueChangedFcn', @(~, ~) refresh());
        ctl.zBounds.Layout.Row = 5; ctl.zBounds.Layout.Column = 1;
        ctl.zBounds.UserData = ctl.zBounds.BackgroundColor;
        stackLabel(g, 7, 'End');
        ctl.zEnd = numField(g, 8, zEnd);
        ctlOf.zoneDiv = {ctl.zStart, ctl.zBounds, ctl.zEnd};
    end

    function makeThreshInputs(g)
        g.RowHeight   = {16, 24, 8, 16, 24, '1x'};
        g.ColumnWidth = {'1x', '1x'};
        g.RowSpacing  = 1; g.ColumnSpacing = 8;
        labs = {'Peaks 1', 'Peaks 2', 'Peaks 3', 'Peaks Final'};
        ids  = {'t1', 't2', 't3', 't4'};
        rr   = [1 1 4 4]; cc = [1 2 1 2];
        ctlOf.threshValues = {};
        for k = 1:4
            stackLabel(g, rr(k), labs{k}, cc(k));
            ctl.(ids{k}) = numField(g, rr(k) + 1, init(15 + k), cc(k));
            ctlOf.threshValues{end+1} = ctl.(ids{k});
        end
    end

    function stackLabel(g, row, txt, col)
        if nargin < 4, col = 1; end
        h = uilabel(g, 'Text', txt, 'FontSize', 11, ...
            'FontColor', [0.33 0.33 0.33], 'VerticalAlignment', 'bottom');
        h.Layout.Row = row; h.Layout.Column = col;
    end

    function h = numField(g, row, val, col)
        if nargin < 4, col = 1; end
        str = charParamsFmtNum(val);
        h = uieditfield(g, 'text', 'Value', str{1}, ...
            'ValueChangedFcn', @(~, ~) refresh());
        h.Layout.Row = row; h.Layout.Column = col;
        h.UserData = h.BackgroundColor;
    end

    function [v, boundsErr] = collectValues()
        [bv, boundsErr] = parseBounds(ctl.zBounds.Value);
        zd = [str2double(ctl.zStart.Value); bv(:); str2double(ctl.zEnd.Value)];
        zd = [zd; -9999 * ones(max(0, 8 - numel(zd)), 1)];
        zd = zd(1:8);
        num = @(id) str2double(ctl.(id).Value);
        dd  = @(id) dropVal(ctl.(id).Value);
        v = [zd; num('yrInterp'); dd('transform'); dd('method'); num('smoYr'); ...
             dd('cPeak'); dd('threshType'); dd('threshMethod'); ...
             num('t1'); num('t2'); num('t3'); num('t4'); num('minCountP'); ...
             num('peakFrequ'); dd('bkgSens'); dd('saveFigures'); ...
             dd('saveData'); dd('allFigures')];
    end

    function chk = checkValues(v, boundsErr)
        E = {}; W = {};       % each row: {key, text}
        if summ.nMissing > 0
            W(end+1, :) = {'data', sprintf(['%d sample(s) have missing values ' ...
                '(volume <= 0); CharAnalysis will interpolate across them.'], ...
                summ.nMissing)};
        end
        if summ.nGaps > 0
            W(end+1, :) = {'data', sprintf(['%d break(s) in sample depths ' ...
                '(cmTop not equal to the previous cmBot); resampling will ' ...
                'occur across them. Consider filling missing rows with -999 ' ...
                'to interpolate instead.'], summ.nGaps)};
        end
        if ~isempty(boundsErr)
            E(end+1, :) = {'zoneDiv', boundsErr};
        end
        if any(isnan(v))
            for k = find(isnan(v))'
                E(end+1, :) = {posKey{k}, sprintf('%s is empty or not a number.', ...
                    posLabel{k})};                                  %#ok<AGROW>
            end
            chk = struct('E', {E}, 'W', {W});
            return
        end
        zd = v(1:8); zd = zd(zd ~= -9999);
        if numel(zd) > 2
            inner = zd(2:end-1);
            if any(inner <= zd(1) | inner >= zd(end))
                E(end+1, :) = {'zoneDiv', 'Zone boundaries must fall between Start and End.'};
            end
        end

        P.zoneDiv = zd; P.yrInterp = v(9); P.transform = v(10);
        S.method = v(11); S.yr = v(12);
        PA.cPeak = v(13); PA.threshType = v(14); PA.threshMethod = v(15);
        PA.threshValues = v(16:19); PA.minCountP = v(20);
        PA.peakFrequ = v(21); PA.bkgSens = v(22);
        R.saveFigures = v(23); R.save = v(24); R.allFigures = v(25);
        [errMsg, warnMsgs] = runValidate(charData, P, S, PA, R);
        if ~isempty(errMsg)
            m = regexprep(errMsg, '^CharValidateParams:\s*', '');
            E(end+1, :) = {keyOfMessage(m), m};
        end
        for k = 1:numel(warnMsgs)
            m = regexprep(warnMsgs{k}, '^CharValidateParams:\s*', '');
            W(end+1, :) = {keyOfMessage(m), m};                    %#ok<AGROW>
        end

        if numel(zd) >= 2
            zEndV = zd(end);
            if zd(1) < summ.ageTop
                W(end+1, :) = {'zoneDiv', sprintf(['Start (%g) is younger than ' ...
                    'the top of the record (%g).'], zd(1), summ.ageTop)};
            end
            if zd(1) > summ.ageTop
                W(end+1, :) = {'zoneDiv', sprintf(['Start (%g) is older than the ' ...
                    'top of the record (%g); samples above it will be excluded.'], ...
                    zd(1), summ.ageTop)};
            end
            if zEndV > summ.ageBot
                W(end+1, :) = {'zoneDiv', sprintf(['End (%g) is older than the ' ...
                    'bottom of the record (%g); CharAnalysis will reset it to %g.'], ...
                    zEndV, summ.ageBot, summ.ageBot)};
            end
            if zEndV < summ.ageBot
                W(end+1, :) = {'zoneDiv', sprintf(['End (%g) is younger than the ' ...
                    'bottom of the record (%g); samples below it will be excluded.'], ...
                    zEndV, summ.ageBot)};
            end
        end
        if v(9) > 0 && v(9) < summ.resMedian
            W(end+1, :) = {'yrInterp', sprintf(['The interpolation interval (%g yr) ' ...
                'is finer than the median sample resolution (%.1f yr), so ' ...
                'interpolated values will not be independent measurements.'], ...
                v(9), summ.resMedian)};
        end
        chk = struct('E', {E}, 'W', {W});
    end

    function ch = changeList(v)
        ch = {};
        if isempty(old), return, end
        for k = 1:25
            a = old(k); b = v(k);
            if xor(isnan(a), isnan(b)) || (~isnan(a) && ~isnan(b) && a ~= b)
                sa = charParamsFmtNum(a); sb = charParamsFmtNum(b);
                if isempty(sb{1}), sb{1} = '(empty)'; end
                if isempty(sa{1}), sa{1} = '(empty)'; end
                ch{end+1} = sprintf('%s: %s -> %s', posLabel{k}, sa{1}, sb{1}); %#ok<AGROW>
            end
        end
    end

    function refresh()
        [v, boundsErr] = collectValues();
        chk = checkValues(v, boundsErr);

        % Reset highlights
        keys = fieldnames(ctlOf);
        for k = 1:numel(keys)
            hs = ctlOf.(keys{k});
            for m = 1:numel(hs)
                hs{m}.BackgroundColor = hs{m}.UserData;
            end
        end
        % Notes by tab
        tabNames = {'Data', 'Pretreatment', 'Smoothing', 'Peaks', 'Results'};
        txt = struct();
        for k = 1:numel(tabNames), txt.(tabNames{k}) = {}; end
        for k = 1:size(chk.W, 1)
            key = chk.W{k, 1}; t = keyTab.(key);
            txt.(t){end+1} = ['Note: ' chk.W{k, 2}];
            highlight(key, colNote);
        end
        for k = 1:size(chk.E, 1)
            key = chk.E{k, 1}; t = keyTab.(key);
            txt.(t) = [{['Must fix: ' chk.E{k, 2}]}, txt.(t)];
            highlight(key, colErr);
        end
        for k = 1:numel(tabNames)
            if isempty(txt.(tabNames{k}))
                notes.(tabNames{k}).Value = {'No notes for this tab.'};
            else
                notes.(tabNames{k}).Value = txt.(tabNames{k});
            end
        end

        % Review tab
        rv = {};
        if size(chk.E, 1) > 0
            rv{end+1} = sprintf('%d problem(s) must be fixed before saving.', size(chk.E, 1));
        else
            rv{end+1} = 'No problems prevent saving.';
        end
        rv{end+1} = '';
        for k = 1:size(chk.E, 1)
            rv{end+1} = sprintf('Must fix (%s): %s', keyTab.(chk.E{k, 1}), chk.E{k, 2}); %#ok<AGROW>
        end
        for k = 1:size(chk.W, 1)
            rv{end+1} = sprintf('Note (%s): %s', keyTab.(chk.W{k, 1}), chk.W{k, 2}); %#ok<AGROW>
        end
        if strcmp(mode, 'revise')
            ch = changeList(v);
            rv{end+1} = '';
            rv{end+1} = 'Changes from the existing file:';
            if isempty(ch)
                rv{end+1} = '  No values changed; the file will not be rewritten.';
            else
                rv = [rv, strcat({'  '}, ch)];
            end
        end
        txtReview.Value = rv;
    end

    function highlight(key, col)
        if isfield(ctlOf, key)
            hs = ctlOf.(key);
            for m = 1:numel(hs)
                if ~isequal(hs{m}.BackgroundColor, colErr)   % errors win
                    hs{m}.BackgroundColor = col;
                end
            end
        end
    end

    function onSave()
        [v, boundsErr] = collectValues();
        chk = checkValues(v, boundsErr);
        if size(chk.E, 1) > 0
            msg = cell(1, size(chk.E, 1));
            for k = 1:size(chk.E, 1)
                msg{k} = sprintf('%s: %s', keyTab.(chk.E{k, 1}), chk.E{k, 2});
            end
            uialert(fig, strjoin(msg, [newline newline]), ...
                'Fix these problems before saving');
            return
        end
        if strcmp(mode, 'revise')
            ch = changeList(v);
            if isempty(ch)
                fprintf('No parameter values changed; %s was not rewritten.\n', ...
                    [site '_charParams.csv'])
                finish(paramsPath);
                return
            end
            uiconfirm(fig, sprintf(['The previous values will be kept in a new, ' ...
                'time-stamped Parameter Log column.\n\n%s'], strjoin(ch, newline)), ...
                'Confirm changes', 'Options', {buttonLabel, 'Go back'}, ...
                'DefaultOption', 1, 'CancelOption', 2, ...
                'CloseFcn', @(~, evt) onConfirm(evt, v));
        else
            doSave(v);
        end
    end

    function onConfirm(evt, v)
        if strcmp(evt.SelectedOption, buttonLabel)
            doSave(v);
        end
    end

    function doSave(v)
        try
            if strcmp(mode, 'new')
                CharParamsNew(v, paramsPath);
            else
                CharParamsRevise(paramsPath, v);
            end
        catch ME
            uialert(fig, ME.message, 'Could not save');
            return
        end
        finish(paramsPath);
    end

    function finish(p)
        result = p;
        uiresume(fig);
    end

    function onCancel()
        result = '';
        uiresume(fig);
    end

end % CharParamsGUI

%% ======================================================================
%% Local functions
%% ======================================================================

function p = charAbsPath(p)
% Make a path absolute relative to the current folder.
if ~(numel(p) >= 2 && p(2) == ':') && ~startsWith(p, '\\') && ~startsWith(p, '/')
    p = fullfile(pwd, p);
end
end

function s = charDataSummary(charData)
% Summary of the charcoal record used for display, suggestions, and notes.
n   = size(charData, 1);
res = diff(charData(:, 3));
s.n          = n;
s.cmTop      = charData(1, 1);  s.cmBot  = charData(n, 2);
s.ageTop     = charData(1, 3);  s.ageBot = charData(n, 4);
s.resMedian  = median(res);
s.resMin     = min(res);        s.resMax = max(res);
s.resSuggest = round(median(res));
s.nMissing   = sum(charData(:, 5) <= 0);
s.nGaps      = sum((charData(2:end, 1) - charData(1:end-1, 2)) ~= 0);
end

function v = dropVal(x)
% Map the "(missing in file)" dropdown sentinel back to NaN.
v = x;
if v == -99999
    v = NaN;
end
end

function [vals, err] = parseBounds(txt)
% Comma-separated zone boundaries -> numeric vector (or an error message).
vals = []; err = '';
txt = strtrim(txt);
if isempty(txt), return, end
parts = strtrim(strsplit(txt, {',', ';'}));
parts = parts(~cellfun(@isempty, parts));
num   = str2double(parts);
if any(isnan(num))
    err = sprintf(['Zone boundaries must be numbers separated by commas; ' ...
        'could not read: %s.'], strjoin(parts(isnan(num)), ', '));
    return
end
if numel(num) > 6
    err = sprintf('At most six zone boundaries can be entered (found %d).', numel(num));
    return
end
vals = num(:);
end

function key = keyOfMessage(msg)
% Which parameter row a CharValidateParams message belongs to.
rules = {'smoothing window', 'smoYr'; 'smoothing method', 'method'; ...
         'zoneDiv', 'zoneDiv'; 'yrInterp', 'yrInterp'; ...
         'locally defined threshold', 'threshMethod'; ...
         'threshValues', 'threshValues'; 'threshType', 'threshType'; ...
         'threshMethod', 'threshMethod'; 'cPeak', 'cPeak'; ...
         'minCountP', 'minCountP'; 'peakFrequ', 'peakFrequ'};
key = 'data';
for k = 1:size(rules, 1)
    if ~isempty(strfind(lower(msg), lower(rules{k, 1})))   %#ok<STREMP>
        key = rules{k, 2};
        return
    end
end
end

function [errMsg, warnMsgs] = runValidate(charData, P, S, PA, R)
% Run CharValidateParams, capturing its error and warnings without printing.
errMsg = ''; warnMsgs = {};
ws = warning;
warning('on', 'all'); warning('off', 'backtrace');
lastwarn('');
try
    out = evalc('CharValidateParams(charData, P, S, PA, R)');
catch ME
    out = '';
    errMsg = ME.message;
end
warning(ws);
tok = regexp(out, 'Warning: ([^\n]*)', 'tokens');
warnMsgs = cellfun(@(c) strtrim(c{1}), tok, 'UniformOutput', false);
if isempty(warnMsgs)
    lw = lastwarn;
    if ~isempty(lw)
        warnMsgs = {lw};
    end
end
end

function lines = wrapText(str, width)
% Wrap a string into a cell array of lines of at most width characters.
words = strsplit(strtrim(str), ' ');
lines = {};
cur = '';
for k = 1:numel(words)
    w = words{k};
    if isempty(cur)
        cur = w;
    elseif numel(cur) + 1 + numel(w) <= width
        cur = [cur ' ' w];                                     %#ok<AGROW>
    else
        lines{end+1} = cur;                                    %#ok<AGROW>
        cur = w;
    end
end
if ~isempty(cur)
    lines{end+1} = cur;
end
if isempty(lines)
    lines = {''};
end
end
