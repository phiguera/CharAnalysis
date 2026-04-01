function CharFigureMenu(r, figSelection, saveFigs)
% CharFigureMenu   Interactive figure menu for the CharAnalysis modular path.
%
%   CharFigureMenu(results)
%
%   Presents a command-window menu allowing the user to select which
%   CharAnalysis figures to generate. All selected figures are generated
%   at once. An optional save step saves only the figures generated in
%   the current session.
%
%   Called automatically when CharAnalysis is run with the 'modular' option:
%       CharAnalysis('mysite.csv', 'modular')
%
%   Figures available:
%       3  - C_int, C_back, and C_peak
%       4  - Sensitivity to alternative thresholds and SNI
%       5  - Cumulative peaks through time
%       6  - FRI distributions by zone
%       7  - Continuous fire history
%       8  - Between-zone comparisons
%       9  - Alternative threshold displays

%% -- Menu definition --------------------------------------------------------
menuItems = { ...
    '3',  'Figure 3  -- C_int, C_back, and C_peak'; ...
    '4',  'Figure 4  -- Sensitivity to alternative thresholds and SNI'; ...
    '5',  'Figure 5  -- Cumulative peaks through time'; ...
    '6',  'Figure 6  -- FRI distributions by zone'; ...
    '7',  'Figure 7  -- Continuous fire history'; ...
    '8',  'Figure 8  -- Between-zone comparisons'; ...
    '9',  'Figure 9  -- Alternative threshold displays'; ...
};
nItems = size(menuItems, 1);

%% -- Display menu -----------------------------------------------------------
disp(' ')
disp('=================================================================')
disp('  CharAnalysis -- Modular Figure Menu')
disp('=================================================================')
disp('  Enter figure numbers separated by spaces to generate figures,')
disp('  or enter ''all'' to generate all figures.')
disp(' ')
for i = 1:nItems
    fprintf('    %s  --  %s\n', menuItems{i,1}, menuItems{i,2})
end
disp(' ')
disp('  Examples:')
disp('    3 5 7      (generate Figures 3, 5, and 7)')
disp('    all        (generate all figures)')
disp('=================================================================')

%% -- Get user selection -----------------------------------------------------
if ~isempty(figSelection)
    % Programmatic call -- skip interactive menu
    selectedNums = unique(figSelection);
    fprintf('  Figure selection: %s\n', num2str(selectedNums))
else
    % Interactive call
    selection = strtrim(input('  Your selection: ', 's'));
    if strcmpi(selection, 'all')
        selectedNums = cellfun(@str2double, menuItems(:,1))';
    else
        tokens = strsplit(selection);
        selectedNums = [];
        for i = 1:length(tokens)
            n = str2double(tokens{i});
            if ~isnan(n) && ismember(num2str(round(n)), menuItems(:,1))
                selectedNums(end+1) = round(n); %#ok<AGROW>
            else
                fprintf('  Skipping unrecognised entry: %s\n', tokens{i})
            end
        end
        selectedNums = unique(selectedNums);
    end
end

%% -- Generate selected figures ----------------------------------------------
disp(' ')
disp('  Generating figures...')
for figNum = selectedNums
    switch figNum
        case 3
            fprintf('  Generating Figure 3...\n')
            CharPlotFig3_CintCbackCpeak(r)
        case 4
            fprintf('  Generating Figure 4...\n')
            CharPlotFig4_ThresholdSNI(r)
        case 5
            fprintf('  Generating Figure 5...\n')
            CharPlotFig5_CumulativePeaks(r)
        case 6
            fprintf('  Generating Figure 6...\n')
            CharPlotFig6_FRIDistributions(r)
        case 7
            fprintf('  Generating Figure 7...\n')
            CharPlotFig7_ContinuousFireHistory(r)
        case 8
            fprintf('  Generating Figure 8...\n')
            CharPlotFig8_ZoneComparisons(r)
        case 9
            fprintf('  Generating Figure 9...\n')
            CharPlotFig9_ThresholdDetails(r)
    end
end
disp('  Done.')

%% -- Offer to save ----------------------------------------------------------
disp(' ')
if ~isempty(saveFigs)
    saveChoice = saveFigs;
else
    saveChoice = strcmpi(strtrim(input('  Save generated figures as PDF and TIFF? (y/n): ', 's')), 'y');
end

if saveChoice
    figNames = containers.Map( ...
        {3, 4, 5, 6, 7, 8, 9}, ...
        {'03_CHAR_analysis', '04_CHAR_peak_sens', '05_cum_peaks_through_time', ...
         '06_FRI_dists', '07_continuous_fire_hx', '08_CHAR_dists', ...
         '09_threshold_details'});
    disp('  Saving figures...')
    for figNum = selectedNums
        baseName = figNames(figNum);
        figure(figNum)
        set(gcf, 'PaperPositionMode', 'auto', 'PaperType', 'uslegal')
        orient(gcf, 'landscape')
        print('-dpdf', '-r300', [baseName, '.pdf'])
        print('-dtiff', '-r300', [baseName, '.tif'])
        fprintf('  Saved: %s\n', baseName)
    end
    disp('  All figures saved.')
else
    disp('  Figures not saved.')
end

disp(' ')
disp('  CharFigureMenu complete.')
disp(' ')

end
