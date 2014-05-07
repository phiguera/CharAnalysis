function [CharAnalysisResults] = CharAnalysis (fileName)
% CharAnalysis   Analyze charcoal time series based on selected parameters.
%   
%    *****************************************************************  
%                           CharAnalysis 1.1                       
%                      (c) 2004 - 2010, P.E. Higuera                   
%              MATLAB (r). (c) 1984 - 2010 The MathWorks, Inc.           
%                        phiguera@uidaho.edu                   
%                       Please read documentation at:              
%                  http://code.google.com/p/charanalysis/               
%    '*****************************************************************
% CharAnalysis requires an input file in .xls or .csv format. This file 
% includes the selected parameters for peak analysis, and it either 
% includes (.xls) or references (.csv) the input charcoal dataset. See the 
% template file, templateChar.xls, for details. 

disp('                                                                 ') 
disp('*****************************************************************')  
disp('                       CharAnalysis 1.1                          ')
disp('                   (c) 2004 - 2010, P.E. Higuera                 ')
disp('          MATLAB(r). (c) 1984 - 2010 The MathWorks, Inc.         ')
disp('                       phiguera@uidaho.edu                       ')
disp('                    Please read documentation:                   ')
disp('              http://code.google.com/p/charanalysis/             ')
disp('*****************************************************************')  
disp('                                                                 ')  

if nargin == 0
  disp('  CharAnalysis requires an input file in .xls or .csv format.    ')
  disp('   This file includes the selected parameters for peak analysis, ')
  disp('   and it either includes (.xls) or references (.csv) the input  ')
  disp('   charcoal dataset. The input file must be in the working       ')
  disp('   directory for the program to run.                             ')
  disp('                                                                 ')
  disp('  If you choose to save the results (figures and/or data), they  ')
  disp('   will be saved in the directory containing the input file.     ')
  disp('   *NOTE*: you must close the input file before running          ')
  disp('   CharAnalysis for the .xls or .csv file to be updated with new ')
  disp('   results.                                                      ')
  disp('                                                                 ')
  disp('*****************************************************************')
  fileName = inputwd('\n Input the file name OR the full path to the site directory, \n bounded with single quotations and including the file \n extension (if file name): ','');
end
warning off all
clear charData

%% HANDLE DIRECTORY NAME INPUT
if (~min(fileName(end-2:end) == ['xls']) ... % If not using an .xls file...
    && ~min(fileName(end-2:end) == ['csv'])) % If not using an .xls file,
    % then assume it's a directory.
    cd(fileName); % Change to site directory.
    i1 = regexp(fileName,'\'); % Find where site directory starts in path (PC).
    i2 = regexp(fileName,'/'); % Find where site directory starts in path (MAC).
    i = max([i1 i2]);
    site = fileName((i(end) + 1):end);      % Set site name.
    fileName = [site '_CharParams.csv'];    % Set params file name based on 
        % site name.
end

%% READ INPUT FILE AND CREATE PARAMETERS VARIABLE
disp ('(1) Reading input file...')
[charData, Pretreatment, Smoothing, PeakAnalysis,...
    Results, site] = CharParameters(fileName);

%% SCREEN FOR INCOMPATABLE PARAMETERS
if PeakAnalysis.threshMethod == 1 && PeakAnalysis.threshType > 1
    error ('A locally defined threshold cannot be user defined. Change input parameters.')
end
if PeakAnalysis.threshMethod == 2 &&...
       (sum(PeakAnalysis.threshValues(1))>1 ||...
       sum(PeakAnalysis.threshValues(2))>1 ||...
       sum(PeakAnalysis.threshValues(3))>1 ||...
       sum(PeakAnalysis.threshValues(4))>1)
    error ('Threshold values must be < 1.0 when using data-defined thresholds. Change input parameters.')
end
disp('      ...done.')

%% RESAMPLE DATA, CALCULATE CHAR, AND TRANSFORM DATA (IF SELECTED).
disp ('(2) Pretreating charcoal data...') 
[Charcoal Pretreatment gapIn] = CharPretreatment(charData,site,...
    Pretreatment, Results);
disp('      ...done.')

%% SMOOTH CHARCOAL RECORD TO ESTIMATE LOW-FREQUENCY TRENDS
disp('(3) Smoothing resampled CHAR to estimate low-frequency trends')
disp('    and calculating peak CHAR...')
[Charcoal] = CharSmooth (Charcoal,Pretreatment,Smoothing,...
    Results);
if min(Charcoal.accIS) == 0 && PeakAnalysis.cPeak == 2
    error ('Cannot calculate C_peak when C_background values = 0; change parameters.')
end

%% CALCULATE PEAK CHAR COMPONENT BY REMOVING BACKGROUND CHAR
if PeakAnalysis.cPeak == 1
    Charcoal.peak = Charcoal.accI - Charcoal.accIS; % Residual charcoal.
else
    Charcoal.peak = Charcoal.accI ./ Charcoal.accIS;% Standardized charcoal.
end
disp('      ...done.')

%% DEFINE POSSIBLE THRESHOLD FOR PEAK IDENTIFICATION
disp('(4) Defining possible thresholds for peak identification...')
if  PeakAnalysis.threshType == 1   % If threshold is defined globally.
    [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment,...
    PeakAnalysis, site, Results);
end   
if  PeakAnalysis.threshType == 2  % If threshold is defined locally.
    [CharThresh] = CharThreshLocal(Charcoal,...
    Smoothing, PeakAnalysis, site, Results);
end
disp('      ...done.')

%% IDENTIFY CHAROCAL PEAKS BASED ON POSIBLE THRESHOLDS
disp('(5) Identifying peaks based on possible thresholds...')
[Charcoal, CharThresh] = CharPeakID (Charcoal,Pretreatment,PeakAnalysis,...
    CharThresh);
disp('      ...done.')

%% PLOT RESULTS FROM CharPeakID
if Results.save == 1
    disp('(6) Plotting results and saving data...')
else
    disp('(6) Plotting results...')
end
[Charcoal] = CharPeakAnalysisResults (Charcoal, Pretreatment,...
    PeakAnalysis, CharThresh, Smoothing, site, Results, fileName, gapIn);
disp('      ...CharAnalysis finished.')

if Results.save == 1
    inputFilePath = cd;
    disp('                            ')
    if min(fileName(end-2:end) == ['xls']) > 0  % If using an .xls file...
        disp('   Results saved to input file:')
        disp(['     ' cd '\' fileName])
    else    % If using .csv file
        disp('   Results saved to file:')
        disp(['     ' cd '\' site '_CharResults.csv'])
    end
    disp('                            ')
end

%% CREATE VARIABLES TO RETURN TO THE WORKSPACE
CharAnalysisResults.Charcoal = Charcoal;
CharAnalysisResults.CharThresh = CharThresh;
CharAnalysisResults.Pretreatment = Pretreatment;
CharAnalysisResults.Smoothing = Smoothing;
CharAnalysisResults.PeakAnalysis = PeakAnalysis;
CharAnalysisResults.Results = Results;

%% RUN BACKGROUND CHAR SENSITIVITY, IF SELECTED
if PeakAnalysis.bkgSens == 1
%     commandwindow
    disp('(7) Running C_background sensitivity analysis:')
    [z, GOF_i, SNI_i] = bkgCharSensitivity (Charcoal, CharThresh,...
    PeakAnalysis, Pretreatment, Smoothing, Results, site);
    disp('C_background sensitivity analysis finished')
end

%% ORGANIZE FIGURES
    if  PeakAnalysis.bkgSens == 1
        figIn = 1:10;
    else
        figIn = 1:9;
    end
    if  Results.allFigures ~= 1
        figIn([1,2,9]) = [];
    end
    for i = 1:length(figIn)
        figure (figIn(i))
    end