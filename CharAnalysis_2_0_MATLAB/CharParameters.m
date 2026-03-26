function [charData, Pretreatment, Smoothing, PeakAnalysis, ...
          Results, site, plotData] = CharParameters(fileName)
% CharParameters    Retrieve parameters and charcoal data from input file.
%   [charData, Pretreatment, Smoothing, PeakAnalysis,
%    Results, site, plotData] = CharParameters(fileName)
%
%   Reads the .xls/.xlsx or .csv input file and returns the charcoal
%   dataset and all analysis parameters for use in CharAnalysis.
%
%   INPUTS
%     fileName : path to the .xls/.xlsx or .csv parameter file
%
%   OUTPUTS
%     charData     : raw charcoal data matrix (samples x 6 columns)
%     Pretreatment : struct  (zoneDiv, yrInterp, transform)
%     Smoothing    : struct  (method, yr)
%     PeakAnalysis : struct  (cPeak, threshType, threshMethod,
%                             threshValues, minCountP, peakFrequ, bkgSens)
%     Results      : struct  (save, saveFigures, allFigures)
%     site         : char string, site name read from cell G1 of charData
%     plotData     : scalar flag (1 = generate diagnostic plots).
%                    Returned as an explicit output in v2.0; was a global
%                    variable in v1.x.  Pass this value forward to every
%                    function that previously read it from global scope.
%
%   v2.0 changes vs v1.1
%     - xlsread()  replaced by readmatrix() / readcell()  (R2019a+).
%     - global plotData   removed; returned as output argument instead.
%     - global bkgSensIn  removed entirely; initialised as a plain local
%       variable in CharAnalysis.m.
%     - Duplicated parameter-unpack block (one per file type) consolidated
%       into a single shared block executed after both read paths.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── READ DATA AND RAW PARAMETER VECTOR ───────────────────────────────────
%
%   After this block, regardless of file type, two variables exist:
%     charData   : numeric matrix  (samples x 6+)
%     charParams : numeric vector  (25 x 1)
%     site       : char string

if min(fileName(end-2:end) == 'xls') > 0
    %% ---- Excel input (.xls or .xlsx) ------------------------------------

    % Data sheet
    charData = readmatrix(fileName, 'Sheet', 'charData');

    % Parameter column: rows 2-26, column C (col 3)
    % readmatrix returns the numeric values; text rows come back as NaN
    charParams = readmatrix(fileName, 'Sheet', 'charParams', ...
                            'Range', 'C2:C26');

    % Site name lives in a text cell (G1 of charData sheet)
    raw = readcell(fileName, 'Sheet', 'charData', 'Range', 'G1:G1');
    if ~isempty(raw) && ~ismissing(raw{1})
        site = char(raw{1});
    else
        site = 'UnknownSite';
        warning('CharParameters: site name not found in cell G1 of charData sheet. Using ''UnknownSite''.');
    end

else
    %% ---- CSV input ------------------------------------------------------
    %
    %   Convention (unchanged from v1.1):
    %     <siteName>_CharParams.csv   holds the parameter rows
    %     <siteName>_charData.csv     holds the charcoal data
    %   The site name is the filename stem before '_CharParams'.

    charData = csvread([fileName(1:end-15) '_charData.csv'], 1, 0);
    site     = fileName(1:end-15);

    % Read parameter file: skip header row, grab col 2 (the numeric value)
    % and col 1 (the parameter name, discarded).
    fid  = fopen(fileName);
    temp = textscan(fid, '%*q %q %n %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q %*q', ...
                    25, 'delimiter', ',', 'HeaderLines', 1);
    fclose(fid);
    charParams = temp{2};

end

%% ── UNPACK PARAMETER VECTOR INTO STRUCTS ─────────────────────────────────
%
%   Row assignments match the charParams worksheet layout (rows 2-26 in the
%   .xls template, positions 1-25 in the vector here).  This single block
%   replaces the two duplicated unpack blocks from v1.1.

%   -- Pretreatment (rows 1-10) --
Pretreatment.zoneDiv  = charParams(1:8);

% .xls path: xlsread left unused zoneDiv slots as NaN
Pretreatment.zoneDiv(isnan(Pretreatment.zoneDiv)) = [];

% .csv path: unused slots are filled with the sentinel value -9999
Pretreatment.zoneDiv(Pretreatment.zoneDiv == -9999) = [];

Pretreatment.yrInterp  = charParams(9);
Pretreatment.transform = charParams(10);

%   -- Smoothing (rows 11-12) --
Smoothing.method = charParams(11);
Smoothing.yr     = charParams(12);

%   -- PeakAnalysis (rows 13-22) --
PeakAnalysis.cPeak        = charParams(13);
PeakAnalysis.threshType   = charParams(14);
PeakAnalysis.threshMethod = charParams(15);
PeakAnalysis.threshValues = charParams(16:19);
PeakAnalysis.minCountP    = charParams(20);
PeakAnalysis.peakFrequ    = charParams(21);
PeakAnalysis.bkgSens      = charParams(22);

%   -- Results (rows 23-25) --
Results.saveFigures = charParams(23);
Results.save        = charParams(24);

if numel(charParams) >= 25 && ~isnan(charParams(25))
    Results.allFigures = charParams(25);
else
    Results.allFigures = 1;     % default: show all diagnostic figures
end

%% ── PLOT FLAG ─────────────────────────────────────────────────────────────
%
%   In v1.1 this was set as a global variable here and read by CharSmooth,
%   CharPretreatment, CharThreshGlobal, and CharThreshLocal.
%   In v2.0 it is a plain output argument passed forward explicitly so that
%   none of those functions need global scope.
%
%   Value is always 1 on a normal run; bkgCharSensitivity sets its own
%   local copy to 0 to suppress per-iteration diagnostic plots.

plotData = 1;

end
