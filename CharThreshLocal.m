function [CharThresh] = CharThreshLocal(Charcoal,...
    Smoothing, PeakAnalysis, site, Results)
% CharThreshLocal  Calculate local threshold value for charcoal peak id.
%   [CharThresh] = CharThreshlocal(Charcoal, Smoothing, PeakAnalysis,...
%   site, Results)
%
%   Determines a threshold value for each interpolated sample, based on the
%   distribution of CHAR values within the selected window (yr) and either
%   a Gaussian mixture model or the assumption that the noise component 
%   of the peak charcoal record (C_peak) is normally distributed around 0
%   (if C_peak is defined by residuals) or 1 (if C_peak is defined by
%   ratios). 


%% CREATE LOCAL VARIABLES
threshYr = Smoothing.yr; % [yr] Years over which to define threshold.
figPosition = [0.1178    0.2158     0.8648    0.6943];
global plotData

%% PREALLOCATE SPACE FOR VARIABLES
CharThresh.pos = NaN*ones(length(Charcoal.peak),...
    length(PeakAnalysis.threshValues));  % Space for threshold values.
CharThresh.neg = NaN*ones(length(Charcoal.peak),...
    length(PeakAnalysis.threshValues));    % Space for negative 
    % threshold value.
muHat = NaN*ones(length(Charcoal.peak),2);
    % Space for mean of noise distributions.
sigmaHat = NaN*ones(length(Charcoal.peak),2); % Space for standard 
    % deviation of noise distribution.
propN = NaN*ones(length(Charcoal.peak),2);    % Space for the proportion of
    % each CLUSTER-identified distribution.
CharThresh.SNI = NaN*ones(length(Charcoal.peak),1);  % Space for signal-to-
    % noise index values.
CharThresh.GOF = NaN*ones(length(Charcoal.peak),1); % Space for 
    % goodness-of-fit results.

%% DEFINE INTERNAL VARIABLES FOR PLOTTING
FS = 8;                     % Font size, for graphs
r = mean(diff(Charcoal.ybpI));    % Resolution of record
P = PeakAnalysis.threshValues(4); % Threshold selected

%% SELECT Charcoal.peak VALUES TO EVALUATE, BASED ON Smoothing.yr
    if plotData ~= 0 && Results.allFigures == 1
    plotIn = 1; figure (2); clf; set(gcf,'color','w','name',...
    'Local distributions of C_peak values','units','normalized',...
    'position',[figPosition]); 
    end
for i = 1:length(Charcoal.peak) % For each value in Charcoal.peak, find the
                                % threshold.
    if sum(i/100 == [1:100])
        disp (['    Calculating local threshold ' num2str(i) ' of ' num2str(length(Charcoal.peak))]);
    end
    if i < round(0.5*(threshYr/r))+1  % First 'threshYr' samples.
%         X = Charcoal.peak(1:round(0.5*(threshYr/r))); % Pre June 2009.
%         X = Charcoal.peak(1:round(threshYr/r)); % Modified, June 2009, PEH.
        X = Charcoal.peak(1:round(0.5*(threshYr/r))+i); % Modified, 
                % June 2009, PEH.
    else
        if i > length(Charcoal.peak)-round(0.5*(threshYr/r)) % Last 
                % 'threshYr' samples.
%             X = Charcoal.peak(length(Charcoal.peak)-...
%                 round((threshYr/r)):end);   % Pre June 2009.
            X = Charcoal.peak(i-round(0.5*(threshYr/r)):end);   % Modified,
            % June 2009, PEH. As recommended by RK, this uses samples from 
            % a half-window before i, all the way to end of record. 
        else
            X = Charcoal.peak(i-round(0.5*(threshYr/r)):...
                i+round(0.5*(threshYr/r)));   % All samples 
                % between first and last 'thrshYr' samples.
        end
    end

%% ESTIMATE LOCAL NOISE DISTRIBUTION
    if PeakAnalysis.threshMethod == 2 &&...
            PeakAnalysis.cPeak == 1    % Estimate noise distribution with 
            % 0-mean Guassian.
    sigmaHat(i) = std([X(X<=0);abs(X(X<0))]); % Estimated standard 
            % deviation of noise distribution: based on Charcoal.peak
            % values <= 0.
    muHat(i) = 0;  % Assume that mean = 0.
    else    % Charcoal.peak defiend by ratios, therfore estimate noise 
            % distribution with 1-mean Gaussian.
    sigmaHat(i) = std([[X(X<=1)]-1;abs((X(X<1)-1))]+1);  % Estimated 
            % standard deviation of noise distribution: based on 
            % Charcoal.peak values <= 1.
    muHat(i) = 1;  % Assume that mean = 1.
    end
    if PeakAnalysis.threshMethod == 3 % Estimate noise distribution with
            % Guassian mixture model
    if sum(X) == 0
        disp ('NOTE: All C_peak values = 0; cannot fit noise distribution.')
        disp ('      Mean and standard deviation forced to equal 0.') 
        disp ('      Consider longer smoothing window or alternative for') 
        disp ('      threshMethod parameter.')
        muHat(i,:) = [0 0]; sigmaHat(i,:) = [10^-100 10^-100]; 
            propN(i,:) = [0 0];
    else    
    [muHat(i,:),sigmaHat(i,:),temp,propN(i,:)] =...
        GaussianMixture(X,2,2,false);  % Estimate mean and standard 
        % deviation of noise distribution using the CLUSTER Gaussian 
        % mixture model:
        % http://cobweb.ecn.purdue.edu/~bouman/software/cluster/
        
        % Added, October, 2009, PEH:
        if muHat(1) == muHat(2)
             beep
             disp('WARNING: poor fit of Gaussian mixture model;')
             disp('         re-fit starting with three classes')
             [mu,sig] = GaussianMixture(Charcoal.peak,3,2,false);
        end
   noiseIn = min(find(muHat(i,:) == min(muHat(i,:))));% Index to smaller
        % muHat value, noise.
   signalIn = max(find(muHat(i,:) == max(muHat(i,:))));% Index to larger
        % muHat value, signal.
   muHat(i,:) = [muHat(i,noiseIn) muHat(i,signalIn)];  
   sigmaHat(i,:) = [sigmaHat(i,noiseIn) sigmaHat(i,signalIn)];  
   propN(i,:) = [propN(i,noiseIn) propN(i,signalIn)];  
        % Make noise value first, signal value second.
    end
    end
    
%% DEFINE LOCAL THRESHOLD, SIGNAL-TO-NOISE INDEX, AND GOODNESS-OF-FIT
    % Define range of threshold values, plus threhsold value selected.
    CharThresh.pos(i,:) =  norminv(PeakAnalysis.threshValues',...
        muHat(i,1),sigmaHat(i,1)); 
    
    % Define negative threshold value, based on threhsold value selected.
    CharThresh.neg(i,:) = norminv(1-PeakAnalysis.threshValues,...
        muHat(i,1),sigmaHat(i,1)); 
    
    % Define signal-to-noise index
    sig_i = X(X>norminv(PeakAnalysis.threshValues(4),muHat(i,1),...
        sigmaHat(i,1)));
    noise_i = X(X<=norminv(PeakAnalysis.threshValues(4),...
        muHat(i,1),sigmaHat(i,1)));
    
    % NEW May 2010: SNI calculation based on Kelly et al., in review.
    %   - PEH
    if ~isempty(sig_i)
    CharThresh.SNI(i) = (1/length(sig_i)) .* sum((sig_i - mean(noise_i))./ ...
        std(noise_i)) .* [(length(noise_i)-2)/length(noise_i)];
    else
    CharThresh.SNI(i) = 0;
    end
    
%     Pre-May 2010 method for SNI calculation, based on Higuera et al. 2009
%     (Ecological Monogrphs, 79: 201-213). 
%
%     CharThresh.SNI (i) = var(sig_i) / [var(sig_i) + var(noise_i)];
%      CharThresh.SNI (i) = [nanmean(sig_i) - nanmean(noise_i)]./...
%          std(noise_i);

    % Evaluate goodness-of-fit between modeled noise distribution and 
    % Charcoal.peak data for this time window
        if length(X(X<=norminv(P,muHat(i,1),sigmaHat(i,1)))) > 3
          ksX = X(X<=norminv(P,muHat(i,1),sigmaHat(i,1)));
            ksBin = min(ksX):range(ksX)/100:max(ksX);
            ksCdf = normcdf(ksBin,muHat(i,1),sigmaHat(i,1));
            [ksH,ksP,ksk] = kstest(ksX,[ksBin' ksCdf']);
            CharThresh.GOF(i) = ksP;
        end

%%  PLOT SELECTED Charcoal.peak DISTRIBUTIONS
    if plotData ~= 0 & Results.allFigures == 1
    if sum(i == [round(threshYr/r):...
            round(length(Charcoal.peak)/25):10000]) >= 1 & plotIn < 25
    subplot(5,5,plotIn)
    [n,x] = hist(X,50); 
    H1 = bar(x,n/sum(n),1);set(H1,'facecolor',[.75 .75 .75],...
    'edgecolor',[.75 .75 .75]);
    hold on
    if PeakAnalysis.threshMethod == 3 % If using Gaussian mixture 
       pdf1 = normpdf(x,muHat(i,1),sigmaHat(i,1))*mean(diff(x))*propN(i,1);
       pdf2 = normpdf(x,muHat(i,2),sigmaHat(i,2))*mean(diff(x))*propN(i,2);
       plot(x,pdf1,'k','linewidth',1)   % Gaussian 1
       plot(x,pdf2,'k','linewidth',1)   % Gaussian 2
       plot(x,pdf1+pdf2,'b','linewidth',1)  % Mixed Gaussian
    else
       plot(x,normpdf(x,muHat(i,1),sigmaHat(i,1))*mean(diff(x)), 'k',...
        'linewidth',2)
    end
    plot([norminv(P,muHat(i,1),sigmaHat(i,1)) norminv(P,muHat(i,1),...
        sigmaHat(i,1))],[0 0.5],'r')
    set(gca,'tickdir','out','ylim',[0 0.25],'fontsize',10)
    box off
    x_lim = get(gca,'xlim'); y_lim = get(gca,'ylim');
    
    title([num2str(round(Charcoal.ybpI(mean(i-round(0.5*(threshYr/r)):...
        i+round(0.5*(threshYr/r)))))) ' yr BP'],'fontsize',FS)
    
    text(x_lim(2),y_lim(2),[{['SNI_i = ' num2str(round(...
        CharThresh.SNI(i)*100)/100)]},{['KS p-val ='...
        num2str(round(ksP*100)/100)]},{['t_i = '...
        num2str(round(norminv(P,muHat(i,1),sigmaHat(i,1))*1000)/1000)]}],...
        'horizontalalignment','right','verticalalignment','top',...
        'FontSize',FS);%,'BackgroundColor','w')
    if plotIn == 1
        title([char(site),': ',...
        num2str(round(Charcoal.ybpI(mean(i-round(0.5*(r/Smoothing.yr)):...
        i+round(0.5*(r/Smoothing.yr)))))) ' yr BP'],'fontsize',FS)
    end
    if plotIn == 11
        ylabel ('proportion OR density (scaled)','FontSize',FS)
    end
    if plotIn == 23
        xlabel ('CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    plotIn = plotIn + 1;
    end 
    end% End if plotData > 0
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%% END PLOTTING Charcoal.peak DISTRIBUTIONS
end     % End loop for each Charcoal.peak value

%% SMOOTH THREHSOLDS WITH LOWESS SMOOTHER
CharThresh.SNI (:,1) = smooth(CharThresh.SNI(:,1),...
    Smoothing.yr/r,'lowess');
CharThresh.SNI (CharThresh.SNI < 0,1) = 0;
for i = 1:length(PeakAnalysis.threshValues)
CharThresh.pos(:,i) = smooth(CharThresh.pos(:,i),...
    Smoothing.yr/r,'lowess');
CharThresh.neg (:,i)= smooth(CharThresh.neg(:,i),...
    Smoothing.yr/r,'lowess');
end