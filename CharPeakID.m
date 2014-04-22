function [Charcoal CharThresh] = CharPeakID (Charcoal,Pretreatment,...
    PeakAnalysis,CharThresh)
% CharPeakID    Identify charcoal peaks based on threshold values.
%   [Charcoal CharThresh] = CharPeakID (Charcoal,Pretreatment,...
%    PeakAnalysis,CharThresh)
%
%   Identifies charcoal samples that exceeds threshold value(s) determined 
%   in CharThrshLocal or CharThreshGlobal, and screens these value
%   according to the minimum-count criterion selected. 

%% CREATE LOCAL VARIABLES
r = Pretreatment.yrInterp;

%% PEAK IDENTIFICATION ALGORITHM
% Create space for peaks, Charcoal.charPeaks
if PeakAnalysis.threshType == 1   % If threshold is global.
    CharThreshPossible = CharThresh.possible(CharThresh.possible>0);
    Charcoal.charPeaks = 0*ones(length(Charcoal.peak),...
        length(CharThreshPossible));
    thresholdValues = NaN*ones(size(Charcoal.charPeaks));
    for i = 1:length(CharThreshPossible)
        thresholdValues(:,i) = CharThreshPossible(i);
    end
else                                % Else, threshold is local.
    Charcoal.charPeaks = 0*ones(size(CharThresh.pos));   % Create 
        % Charcoal.charPeaks matrix
    thresholdValues = CharThresh.pos;
end
[l nThresholds] = size(thresholdValues);  % Number of threhsolds

% Flagg values exceeding thresholds
for i = 1:(length(Charcoal.peak))   % For each value in Charcoal.peak.
    for j = 1:nThresholds  % For each threshold value.
        if Charcoal.peak (i,1) > thresholdValues (i,j)  % If Charcoal.peak. 
                                                 % exceeds threshold...
            Charcoal.charPeaks (i,j) = 2;        % Charcoal.charPeaks == 2,
        else
            Charcoal.charPeaks (i,j) = 0;        % else Charcoal.charPeaks
                                                 % == 0. 
        end
    end
end
% Remove consecutive Charcoal.charPeaks:
for i = 1:(length(Charcoal.peak)-1)        % For each Charcoal.peak value                 
    for j = 1:nThresholds           % For each threshold 
                                        % value
        if Charcoal.charPeaks (i,j) > 0 & Charcoal.charPeaks (i+1,j) % If 
                % two consecutive values are > 0 (i.e. flagged as 
                % Charcoal.charPeaks)...
           Charcoal.charPeaks (i,j) = Charcoal.charPeaks (i,j)/2; %...then 
                % make consecutive Charcoal.charPeaks 1, keep first peak as
                % 2.
        else
        end
    end
end
for i = 1:(length(Charcoal.peak))
    for j = 1:nThresholds
        if Charcoal.charPeaks (i,j) < 2    % if value is 0 or 1...
           Charcoal.charPeaks (i,j) = 0;   % make value 0
        else
            Charcoal.charPeaks (i,j) = 1;  % otherwise make value 1; by 
                                    % here Charcoal.charPeaks = 1, and only 
                                    % the first value of a peak is 
                                    % identified.
        end                      
    end
end

% Make a variable to hold the threshold value for each peak identified.
% i.e. instead of 1 / 0 for peak / no peak, 1 is replaced with the 
% threshold value used to identify the peak. This is purely for graphing 
% purposes.
Charcoal.charPeaksThresh = NaN*ones(size(Charcoal.charPeaks)); 
for i = 1:length(Charcoal.peak) % For each Charcoal.peak value
    [in1 in2] = size(thresholdValues);
    for j = 1:in2  % For each threshold
        Charcoal.charPeaksThresh (i,j) = Charcoal.charPeaks (i,j)*...
            thresholdValues (1,j);
           % Multiplying the Charcoal.charPeaks matrix by threshold values.
    end
end

%% MINIMUM-COUNT ANALYSIS
mcWindow = round(150/r)*r; % [yr] Years before and after a peak to look 
    % for the min. and max. value
d = zeros(length(Charcoal.accI),nThresholds);
CharThresh.minCountP = NaN*ones(length(Charcoal.accI),nThresholds);
alphaPeak = PeakAnalysis.minCountP;
for j = 1:nThresholds
    peakIndex = find(Charcoal.charPeaks(:,j)); % Index for the location of 
                                        % identified peaks.
    if length(peakIndex)>1    % Only proceed if there is more than 1 peak.
    for i = 1:length(peakIndex)     % For each peak identified
        peakYr = Charcoal.ybpI(peakIndex);  % Years of charcoal peaks.
        windowTime = [(max(Charcoal.ybpI(Charcoal.ybpI<=peakYr(i)+...
            mcWindow))) (min(Charcoal.ybpI(Charcoal.ybpI>=peakYr(i)-...
            mcWindow)))]; % Years defining peak window.
        windowTime_in = [find(Charcoal.ybpI == windowTime(1))...
            find(Charcoal.ybpI == windowTime(2))]; % Index for years  
            % defining peak window.
        if i == 1   % find the year of the two adjacent Charcoal.charPeaks,
            % unless first peak, then use windowTime(2) as youngest age OR
            % unless last peak, then use windowTime(1) as oldest age.
            windowPeak_in = [find(Charcoal.ybpI ==...
            Charcoal.ybpI(peakIndex(i+1)))...
            find(Charcoal.ybpI == windowTime(2))];
        else
            if i == length(peakIndex)
            windowPeak_in = [find(Charcoal.ybpI == windowTime(1))...
                find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i-1)))];
            else
            windowPeak_in = [find(Charcoal.ybpI ==...
                Charcoal.ybpI(peakIndex(i+1)))...
                find(Charcoal.ybpI == Charcoal.ybpI(peakIndex(i-1)))];
            end
        end
        if windowTime_in(1) > windowPeak_in(1)   % If a peak falls 
                % within the time window defined by mcWindow.
                windowTime_in(1) = windowPeak_in(1);  % Replace the 
                % windowTime_in with the windowPeak_in.
        end
        if windowTime_in(2) < windowPeak_in(2)    % If a peak falls
                % within the time window defined by mcWindow.
                windowTime_in(2) = windowPeak_in(2);  % Replace the 
                % windowTime_in with the windowPeak_in.
        end
        windowSearch = [windowTime_in(1) windowTime_in(2)]; % Final index
        % value for search window: window (1) defines oldest sample and
        % window(2) defines youngest sample; search for max and min 
        % Charcoal.charPeaks within this window.
        countMax = round(max(Charcoal.countI(windowSearch(2):...
            peakIndex(i)))); 
            %[# cm^-3] Max charcoal concentration after peak.
            countMaxIn = windowSearch(2)-1+...
            find(round(Charcoal.countI(windowSearch(2):...
            peakIndex(i))) == countMax);    % Index for location of max 
            % count.
        countMin = round(min(Charcoal.countI(peakIndex(i):...
            windowSearch(1)))); 
            %[# cm^-3] Min charcoal concentration before peak.
            countMinIn = peakIndex(i)-1+...
            find(round(Charcoal.countI(peakIndex(i):...
            windowSearch(1))) == countMin);    % Index for location of max 
                % count.
%             if countMin == 0   % If the minimum count is a 0, then...
%                 countMin = 1;  % Assume a true minimum count of 1; there's
%                     % a 37% chance of obtaining a sample of 0 from a
%                     % Poisson distrubtion with Lambda = 1. This is
%                     % conservative with respect to identifying peaks.
%             end
        volMax = Charcoal.volI(countMaxIn(1));
%         max(Charcoal.volI(windowSearch(2):peakIndex(i))); %[cm^-3] 
            % Sample volume of maximum count.
        volMin = Charcoal.volI(countMinIn(1));
%         min(Charcoal.volI(peakIndex(i):windowSearch(1))); %[cm^-3] 
            % Sample volume of maximum count.           
        d(peakIndex(i),j) = (abs(countMin-(countMin+countMax)*...
            (volMin./(volMin+volMax)))-0.5)/(sqrt((countMin+countMax)...
            *(volMin/(volMin+volMax))*(volMax/(volMin+volMax))));
            % Test statistic.
        CharThresh.minCountP(peakIndex(i),j) =...
            1-tcdf(d(peakIndex(i),j),10^10);  
                            % Inverse of the Student's T cdf at 
                            % CharThresh.minCountP, with 10^10 degrees of 
                            % freedom.
        % From Charster (Gavin 2005):
        % This is the expansion by Shuie and Bain (1982) of the equation by 
        % Detre and White (1970) for unequal 'frames' (here, sediment 
        % volumes). The significance of d is based on the t distribution 
        % with an infinite degrees of freedom, which is the same as the 
        % cumulative normal distribution.     
    end
    end
end

% REMOVE PEAKS THAT DO NOT PASS THE MINIMUM COUNT CRITERION
for j = 1:nThresholds
    in = find(Charcoal.charPeaks(:,j)>0 &...
        CharThresh.minCountP(:,j)>alphaPeak);  % Index for 
        % Charcoal.charPeaks values that also have p-values > alphaPeak. 
    Charcoal.charPeaks(in,j) = 0; % Replace these values with 0.
    Charcoal.charPeaksThresh (in,j) = 0; % Replace these values with 0. 
        % (i.e. no peak)
end

%% CALCULATE SENSITIVITY INDICIES
for j = 1:(nThresholds) % For each threshold value.
    Charcoal.peaksTotal(j) = sum (Charcoal.charPeaks(:,j)); % peaksTotal =  
        % vector of total Charcoal.charPeaks identified for each threshold 
        % value and for each smoothing window.
   inFRI = diff(Charcoal.ybpI(Charcoal.charPeaks(:,j)>0)); % FRIs for 
        % threshold (j).
   if length(inFRI) > 0 
        Charcoal.threshFRI(1:length(inFRI),j)=inFRI; % FRIs for each 
            % threshold.
   end
end

%% FIND PEAKS AT EACH TREHSOLD, IF THRESHOLD IS DEFINED GLOBALLY
% if PeakAnalysis.threshType == 1
% threshIn(1) = find(CharThreshPossible == CharThresh.pos(1,1));
% threshIn(2) = find(CharThreshPossible == CharThresh.pos(1,2));
% threshIn(3) = find(CharThreshPossible == CharThresh.pos(1,3));
% threshIn(4) = find(CharThreshPossible == CharThresh.pos(1,4));
% Charcoal.charPeaks = Charcoal.charPeaks(:,threshIn);
% Charcoal.charPeaksThresh = Charcoal.charPeaksThresh(:,threshIn);
% end