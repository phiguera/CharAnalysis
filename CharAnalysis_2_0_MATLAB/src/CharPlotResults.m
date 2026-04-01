function CharPlotResults(r)
% CharPlotResults   Generate all output figures for CharAnalysis.
%
%   CharPlotResults(results)
%
%   Produces Figures 3-9 (and saves them if Results.saveFigures == 1).
%   This function is a thin wrapper that calls each individual figure
%   function in sequence. All analytical values come from the Post struct
%   produced by CharPostProcess. No computation is performed here.
%
%   Individual figures can also be called directly, for example:
%       results = CharAnalysis('mysite.csv');
%       CharPlotFig7_ContinuousFireHistory(results)
%
%   For interactive figure selection, call CharAnalysis with the modular
%   option:
%       CharAnalysis('mysite.csv', 'modular')

%% -- Unpack needed fields for save logic ------------------------------------
Results = r.Results;

%% -- Generate figures -------------------------------------------------------
CharPlotFig3_CintCbackCpeak(r)
CharPlotFig4_ThresholdSNI(r)
CharPlotFig5_CumulativePeaks(r)
CharPlotFig6_FRIDistributions(r)
CharPlotFig7_ContinuousFireHistory(r)
CharPlotFig8_ZoneComparisons(r)

if Results.allFigures == 1
    CharPlotFig9_ThresholdDetails(r)
end

%% -- Save figures -----------------------------------------------------------
if Results.saveFigures == 1
    if Results.allFigures == 1
        saveFig(1, '01_pretreatment')
        saveFig(2, '02_threshold_determination')
    end
    saveFig(3, '03_CHAR_analysis')
    saveFig(4, '04_CHAR_peak_sens')
    saveFig(5, '05_cum_peaks_through_time')
    saveFig(6, '06_FRI_dists')
    saveFig(7, '07_continuous_fire_hx')
    saveFig(8, '08_CHAR_dists')
    if Results.allFigures == 1
        saveFig(9, '09_threshold_details')
    end
end

end % end CharPlotResults

%% -- Local helper: save one figure as PDF + TIFF ---------------------------
function saveFig(figNum, baseName)
    figure(figNum)
    set(gcf, 'PaperPositionMode', 'auto', 'PaperType', 'uslegal')
    orient(gcf, 'landscape')
    print('-dpdf', '-r300', [baseName, '.pdf'])
    print('-dtiff', '-r300', [baseName, '.tif'])
end
