function [z, GOF_i, SNI_i] = bkgCharSensitivity (Charcoal, CharThresh,...
    PeakAnalysis, Pretreatment, Smoothing, Results, site)
global plotData
global bkgSensIn

% bkgCharSensitivity    Run CharAnalysis with multiple background windows.
%   [z, GOF_i, SNI_i] = bkgCharSensitivity (Charcoal, CharThresh,...
%   PeakAnalysis, Pretreatment, Smoothing, Results, site)
%
%   Analyzes a charcoal record using multiple smoothing window widths (but
%   all with the same smoothing method) and plots a series of graphs
%   illustrating the sensitivity of results to smoothing windows.

%% CREATE GLOBAL VARIABLE TO LET CharThreshGlobal.m KNOW THIS HAS RUN ONCE
if PeakAnalysis.threshType == 1 % If global threshold.                                
bkgSensIn = 1; %global bkgSensIn
end
plotData = 0;

%% FIND MIN AND MAX BACKGROUND SMOOTHING VALUES, ROUNDED TO NEAREST 100 YR
% Maximum value for background window must be less than the total years in 
% the record.
if max(Charcoal.ybpI) >= 1000
    bkgMax = 1000;
else
    bkgBin = [100:100:900];
    [x,n] = hist(max(Charcoal.ybpI),bkgBin);
    bkgMax = n(x);    % Find largest value for background window.
end

if PeakAnalysis.threshType == 1 % If global threshold.
    bkgMin = 100;
else
% Smallest value for background window such that each local threhsold has
% 30 or more samples:
bkgMin = 30*Pretreatment.yrInterp;
bkgBin = [100:100:500];
[x,n] = hist(bkgMin,bkgBin);
    bkgMin = n(x>0);    % Find smallest value for background window.
end

%% CREATE LOCAL VARIABLES 
if PeakAnalysis.threshType == 1 % If global threshold.
bkSmooth = [100:100:bkgMax];    % [yr] Background smoothing windows to 
                                % evaluate.
y = bkSmooth;  % Possible smoothing window lengths.

figure (2); xlimIn = get(gca,'xlim');
x = xlimIn(1):range(xlimIn)/250:xlimIn(2);  % Possible threshold values.
    x = x(x>0);
z = NaN*ones(length(y),length(x));  % Space for number of peaks identifed
           % at threshold x and background smoothing window y.
SNI_i = NaN*ones(length(y));
GOF_i = NaN*ones(length(y));
else                                 
bkSmooth = [bkgMin:100:bkgMax];     % [yr] Background smoothing windows to 
                                    % evaluate.
x = [1:3];      % Possible threhsold criterion.
y = bkSmooth;   % Possible smoothing window lengths.
z = NaN*ones(length(y),4);  % Space for number of peaks identified.
mFRI = NaN*ones(length(y),4,2);   % Space for mean and std FRI and 95% CIs,
    % as a function of smoothing window. 
SNI_i = NaN*ones(length(y),length(CharThresh.SNI));
GOF_i = NaN*ones(length(y),length(CharThresh.GOF));  % Space for GOF restuls
           % for each smoothing window. 
end

SmoothingI = Smoothing;
    
%% CALCULATE z, NUMBER OF PEAKS, OR KS-GOF RESULTS FOR EACH X-Y COMBO.
for i = 1:length(bkSmooth) % For each potential smoothing window.    
SmoothingI.yr = bkSmooth(i); % Input smoothing window for this iteration.
disp (['    C_background sensitivity iteration ' num2str(i) ' of ' num2str(length(bkSmooth)) ': window = ' num2str(bkSmooth(i)) ' yrs.'])
    %% SMOOTH CHARCOAL RECORD TO ESTIMATE LOW-FREQUENCY TRENDS
    [Charcoal] = CharSmooth (Charcoal,Pretreatment,SmoothingI,...
        Results);
    %% CALCULATE PEAK CHAR COMPONENT BY REMOVING BACKGROUND CHAR
    if PeakAnalysis.cPeak == 1
        Charcoal.peak = Charcoal.accI - Charcoal.accIS; % Residual charcoal.
    else
        Charcoal.peak = Charcoal.accI ./ Charcoal.accIS;% Standardized 
            % charcoal.
    end
    %% DEFINE POSSIBLE THRESHOLD FOR PEAK IDENTIFICATION
    if  PeakAnalysis.threshType == 1   % If threshold is defined globally.
        [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment,...
        PeakAnalysis, site, Results);
    end   
    if  PeakAnalysis.threshType == 2  % If threshold is defined locally.
        [CharThresh] = CharThreshLocal(Charcoal,...
        SmoothingI, PeakAnalysis, site, Results);
    end
    %% IDENTIFY CHARCOAL PEAKS BASED ON POSIBLE THRESHOLDS
    [Charcoal, CharThresh] = CharPeakID (Charcoal,Pretreatment,PeakAnalysis,...
        CharThresh); 
   
    if PeakAnalysis.threshType == 1 % If global threshold.
        z(i,:) = sum(Charcoal.charPeaks);% Sum of charcoal peaks at each
            % threshold value
        GOF_i (i) = CharThresh.GOF(1);   % Goodness-of-fit results for 
            % noise distributions. 
        SNI_i (i) = CharThresh.SNI;
    else                            % If local threshold.
        z(i,:) = sum(Charcoal.charPeaks);% Sum of charocal peaks at each
            % thrshold value. 
        mFRI(i,:,1) = mean(diff(Charcoal.ybpI(Charcoal.charPeaks(:,end)>0)));
        mFRI(i,:,2) = std(diff(Charcoal.ybpI(Charcoal.charPeaks(:,end)>0)));
        GOF_i (i,:) = CharThresh.GOF;    % Goodness-of-fit results for 
            % noise distributions, for each smoothing window.
        SNI_i(i,:) = CharThresh.SNI;     % SNI results for for each 
            % smoothing window
    end
    clear CharThresh
end

%% PLOT RESULTS
figure (10); clf; set(gcf,'color','w','name',...
    'Sensitivity to different background windows','units','normalized',...
    'position',[0.0435    0.0571    0.8649    0.6943])

if PeakAnalysis.threshType == 1 % If global threshold.

contours = [1:10:max(z)];
[X,Y] = meshgrid(x,y);
contourf(X,Y,z)%,contours)
colormap (1-gray)
brighten(.4)
colorbar
x_lim = get(gca,'xlim');
x_lim = [x_lim(1) prctile(x,50)];
xlim (x_lim)
set(gca,'TickDir','out')
xlabel ('threshold (C_p_e_a_k units)')
ylabel ('smoothing window width (yr)')
text(1.2*max(x),min(y)+range(y)/2,'peaks identified (#)','Rotation',270,...
    'HorizontalAlignment','Center')
title ('Peaks identified with varying threshold and C_b_a_c_k_g_r_o_u_n_d criteria')
box off

else                            % If local threshold.

x_lim = [min(bkSmooth)-0.5*mean(diff(bkSmooth)) max(bkSmooth)+...
    0.5*mean(diff(bkSmooth))];
subplot(2,2,1)
boxplot(GOF_i')
set(gca,'XTick',[1:length(bkSmooth)],'XTickLabel',bkSmooth,'TickDir','out')
ylim ([0 1]);
% xlabel ('smoothing window width (yr)')
ylabel ('KS-test results (p-value)')
title ('(a) Noise distribution goodness of fit')
box off

subplot(2,2,2)
boxplot(SNI_i')
set(gca,'XTick',[1:length(bkSmooth)],'XTickLabel',bkSmooth,'yscale',...
    'linear','TickDir','out')
y_lim = ([0 min([max(max(SNI_i)) 10])]);
ylim (y_lim); 
hold on
plot(ones(length(bkSmooth),1)*3,'k--')
% xlabel ('smoothing window width (yr)')
ylabel ('signal-to-noise index')
title ('(b) Signal-to-noise index')
box off

subplot(2,2,3)
plot(bkSmooth,[nanmedian(GOF_i,2)+nanmedian(SNI_i,2)],'-ok','linewidth',2)
xlim(x_lim)
xlabel ('smoothing window width (yr)')
ylabel ('sum of median SNI & GOF')
set(gca,'TickDir','out')
box off
title ('(c) Sum of median SNI and GOF value')

subplot(2,2,4)
errorbar(bkSmooth,mFRI(:,end,1),mFRI(:,end,2),'ko-','linewidth',2)
y_lim = get(gca,'ylim');
hold on
[AX,H1,H2] = plotyy(bkSmooth,mFRI(:,end,1),bkSmooth,z(:,end));
set(AX(1),'ycolor','k','xlim',x_lim)
set(AX(2),'ycolor','b','xticklabel','','xlim',x_lim)
set(H1,'color','k')
set(H2,'color','b','marker','o','linestyle','-')
xlabel ('smoothing window width (yr)')
ylabel('fire-return interval (yr)','color','k')
text (1.05*x_lim(2),mean(range(y_lim)),...
    '# of peaks identified','color','b','rotation',270)
set(gca,'TickDir','out','yscale','linear')
box off
title ('(d) Mean fire-return interval +/- standard deviation')

% subplot(3,6,[5 6 11 12])
% contours = [1:10:max(z)];
% [X Y] = meshgrid(x,y);
% Z = z(:,1:3);
% contourf(Y,X,Z);%,contours)
% colormap (1-gray)
% brighten(.4)
% colorbar
% set(gca,'XTick',bkSmooth,'YTick',[1 2 3],'YTickLabel',...
%     PeakAnalysis.threshValues(1:3))
% ylabel ('threshold criterion (percentile of noise distribution)')
% xlabel ('smoothing window width (yr)')
% text(1150,2,'peaks identified (#)','Rotation',270,...
%     'HorizontalAlignment','Center')
% title ('(c) Peaks identified with varying threshold criteria')
% box off
end

%%  SAVE FIGURE
if Results.saveFigures == 1
    set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 10_sensitivity_to_C_background.pdf
    print -dtiff -r300 10_sensitivity_to_C_background.tif
end