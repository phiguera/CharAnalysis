function [charData, Pretreatment, Smoothing, PeakAnalysis,...
    Results, site] = CharParameters(fileName)
% CharParameters    Retrieve parameters and charcoal data from input file.
%   [charData, Pretreatment, Smoothing, PeakAnalysis,...
%   Results, site] = CharParameters(fileName)
%
%   Reads the parameters and data in the input variable fileName, 
%   referencing the input file in .xls format, and returns parameters and
%   charcoal dataset for use in CharAnalysis.

%% CREATE CHARCOAL DATA AND PARAMETER VARIALBES
if min(fileName(end-2:end) == ['xls']) > 0  % If using an .xls file...
    charData = xlsread(fileName,'charData');
      
    charParams = xlsread(fileName,'charParams','c2:c26');
    Pretreatment.zoneDiv = charParams(1:8);
        Pretreatment.zoneDiv(isnan(Pretreatment.zoneDiv)) = [];    
    Pretreatment.yrInterp = charParams(9);
    Pretreatment.transform = charParams(10);
    Smoothing.method = charParams(11);
    Smoothing.yr = charParams(12);
    PeakAnalysis.cPeak = charParams(13);
    PeakAnalysis.threshType = charParams(14);
    PeakAnalysis.threshMethod = charParams(15);
    PeakAnalysis.threshValues = charParams(16:19);
    PeakAnalysis.minCountP = charParams(20);
    PeakAnalysis.peakFrequ = charParams(21);
    PeakAnalysis.bkgSens = charParams(22);
    Results.save = charParams(24);
    Results.saveFigures = charParams(23);
    if length(charParams) == 25
        Results.allFigures = charParams(25);
    else
        Results.allFigures = 1;
    end

    [trash site] = xlsread(fileName,'charData','g1:g2');
else    % Else input file is in .csv format
   
    charData = csvread([fileName(1:end-15) '_charData.csv'],1,0);
    site = fileName(1:end-15);
    
    fid = fopen(fileName);
    temp = textscan(fid,'%*q%q%n%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q%*q', 25, ...
    'delimiter', ',', 'HeaderLines', 1);
	fclose(fid);
    charParams = temp{2};
    Pretreatment.zoneDiv = charParams(1:8);
        Pretreatment.zoneDiv(Pretreatment.zoneDiv == -9999) = [];    
    Pretreatment.yrInterp = charParams(9);
    Pretreatment.transform = charParams(10);
    Smoothing.method = charParams(11);
    Smoothing.yr = charParams(12);
    PeakAnalysis.cPeak = charParams(13);
    PeakAnalysis.threshType = charParams(14);
    PeakAnalysis.threshMethod = charParams(15);
    PeakAnalysis.threshValues = charParams(16:19);
    PeakAnalysis.minCountP = charParams(20);
    PeakAnalysis.peakFrequ = charParams(21);
    PeakAnalysis.bkgSens = charParams(22);
    Results.save = charParams(24);
    Results.saveFigures = charParams(23);
    if length(charParams) == 25
        Results.allFigures = charParams(25);
    else
        Results.allFigures = 1;
    end
end
    clear trash   
    clear global plotData
    clear global bkgSenesIn
    
    bkgSensIn = 0;
    global bkgSensIn
    
    plotData = 1; 
    global plotData




