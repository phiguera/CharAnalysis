function [FRIyr FRI smFRIyr smFRI smFRIci] = smoothFRI (yr,peaks,winWidth,...
    alpha,nBoot,FRI_param,makeFigure)
% function [FRIyr FRI smFRIyr smFRI smFRIci] = smoothFRI (yr,peaks,winWidth,...
%     alpha,nBoot,FRI_param,makeFigure)
%
% yr = year for each sample.
% peaks = binary seris, same size as yr, with 1 identifying peaks.
% winWidth = yr to smooth over.
% alpha = confidence level for confidence intervals (if selected).
% nBoot = # of bootstrapped samples for confidence intervals.
% FRI_param: 0 = stdev, 1 = 1-alpha% confidence intervals.
% makeFigure: 1 = make figure, 0 = don't make figure.
%
% PEH, July 2010. 
    
%% Input parameters
if nargin < 7
    makeFigure = 0;
end
stepLength = 100;   % [yr] Time step to move window across records. 
maxFRIs = 100;      % Maximum number of FRIs (conservative) possible in one
                    % window.

%% Derrived parameters
steps = -70:stepLength:max(max(yr));        % [cal yr BP] Time steps.
nSteps = length(steps);                     % Number of time steps.
peakYr = yr(peaks == 1);                    % [cal yr BP] Yr with peaks.
spaceIn = round(winWidth/stepLength - 1);
x_lim = [-60 ceil(max(peakYr)/1000)*1000];  % [cal yr BP] X-axis limits.
%% Make spcae for variables
FRImatrix.FRI = NaN*ones(maxFRIs, nSteps);   % Space for matrix  
        % of FRIs(i), for each of n steps (j), for each window (k).
FRImatrix.params = NaN*ones(3,nSteps);   % Space for matrix of
        % mean or median FRI (i = 1) and upper and lower CI (i = 2:3), for
        % each of n steps (j), for each window (k).
FFmatrix = zeros*ones(nSteps,1);   % Space for matrix of
        % mean or median FRI (i = 1) and upper and lower CI (i = 2:3), for
        % each of n steps (j), for each window (k).
p = NaN*ones(nSteps,1);  % Space for p-values for between-step
        % comparisons of FRI distributions for each time step (i)
p2 = NaN*ones(nSteps,1);  % Space for p-values for between-step
        % comparisons of FRI distributions for each time step (i)
p3 = NaN*ones(nSteps,1);  % Space for p-values for between-step
        % comparisons of FRI distributions for each time step (i)      
%% Calculate FRI-distribution parameters and p-values for FRI comparisons.

for i = 1:nSteps%-ceil(winWidth(j)/stepLength)   % For each time step.
        startYr = steps(i)-0.5*winWidth; 
        endYr = steps(i)+0.5*winWidth;
        lengthIn = yr >= startYr & yr < endYr;
        yrIn = find(peakYr >= startYr & peakYr < endYr);
        fris = diff(peakYr(yrIn));  % [yr] FRIs for this i.
        if length(fris) > 1 % If there is more than 1 fri...
%             cumFirefit (i,1:2) = polyfit(peakYr(yrIn),cumPeaks(yrIn),1);
%             cumFirefit (i,3:4) = [startYr endYr];
            FRImatrix.FRI(1:length(fris),i) = fris;
            FFmatrix(i) = length(peakYr(yrIn))/range(yr(lengthIn));
            if FRI_param == 0   % Standard deveation
                FRImatrix.params(1,i) = mean(fris);  
                FRImatrix.params(2,i) = mean(fris)+ [std(fris)/(length(fris)^0.5)];
                FRImatrix.params(3,i) = mean(fris) - [std(fris)/(length(fris)^0.5)];
%                 boot_means = bootstrp(nBoot,'mean',fris);
%                 FRImatrix.params(2:3,i) = prctile(boot_means,...
%                     [alpha/2*100 (1-alpha/2)*100]);
            else % Confidence intervals
                FRImatrix.params(1,i) = mean(fris);  
                boot_means = bootstrp(nBoot,'mean',fris);
                FRImatrix.params(2:3,i) = prctile(boot_means,...
                    [alpha/2*100 (1-alpha/2)*100]); 
            end
            % AD-test between non-overlapping FRI diststributions
            if i > spaceIn %&& i <= nSteps - spaceIn + 1
                x1 = FRImatrix.FRI(:,i-spaceIn); 
                    % More recent FRIs.
                    x1 = x1(x1>-999);
                if ~isempty(x1)
                    x2 = fris;  % Older fris.
                    x12 = NaN*ones(sum([length(x1) length(x2)]),2);
                        x12(1:length(x1),1) = x1;
                        x12(1:length(x1),2) = 1;
                        x12(length(x1)+1:end,1) = x2;
                        x12(length(x1)+1:end,2) = 2;
                    p(i) = AnDarksamtest(x12,alpha);
                    [~, p2(i) K] = kstest2(x1,x2,alpha);
                    p3(i) = ranksum(x1,x2);
                end
            end
        end
end

%% FRI param through time
spaceIn = winWidth/stepLength - 1;

    y_lim = [50 1000];
    in = find (FRImatrix.params(1,:)>-999);
    x = (steps(in));%+winWidth(i)/2);
    smoothIn = (winWidth/stepLength)/length(x);
    y = smooth(FRImatrix.params(1,in),smoothIn,'lowess');
    y4 = FRImatrix.params(1,in);
    y3 = smooth(FRImatrix.params(2,in),smoothIn,'lowess');
    y2 = smooth(FRImatrix.params(3,in),smoothIn,'lowess');
    X = [x'; flipud(x')];
    Y = [y2; flipud(y3)];
    
    % Collect output parameters.
    FRIyr = peakYr(1:end-1);
    FRI = diff(peakYr);
    smFRIyr = x;
    smFRI = y;
    smFRIci = [y3 y2];
    
if makeFigure == 1
figure (10); clf; set(gcf,'color','w')
    plot(peakYr(1:end-1),diff(peakYr),'s','markerfacecolor',[0.6 0.6 0.6],...
        'markeredgecolor',[0.6 0.6 0.6])
    hold on
    plot(x,y,'-k','linewidth',2)
    h = fill(X,Y,[0.8 0.8 0.8]);
    set(h,'edgecolor',[0.8 0.8 0.8])
    set(gca,'xdir','rev','yscale','log','xlim',x_lim,'ylim',y_lim)
    if FRI_param == 0
        ylabel([num2str(winWidth) '-yr mean FRI'])
        legend('FRI','mean FRI +/- 1 stdev.')
    else
        ylabel([num2str(winWidth) '-yr median FRI'])
        legend('FRI','mean FRI',[num2str(100*(1-alpha)) '% CI'])
    end
    plot(peakYr(1:end-1),diff(peakYr),'s','markerfacecolor',[0.6 0.6 0.6],...
        'markeredgecolor',[0.6 0.6 0.6])
    plot(x,y,'-k','linewidth',2)
    plot(x,y4,'-k','linewidth',1)
    grid on
end