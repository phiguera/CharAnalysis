function [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment,...
    PeakAnalysis, site, Results)
% CharThreshGlobal  Calculate global threshold value for charcoal peak id.
%   [CharThresh] = CharThreshGlobal(Charcoal, Pretreatment,...
%   PeakAnalysis, site, Results)
%
%   Determines a threshold value for each interpolated sample, based on the
%   distribution of CHAR values within the entire record and either
%   a Gaussian mixture model or the assumption that the noise component 
%   of the peak charcoal record (C_peak) is normally distributed around 0
%   (if C_peak is defined by residuals) or 1 (if C_peak is defined by
%   ratios). 

%% CREATE LOCAL VARIABLES
zoneDiv = Pretreatment.zoneDiv;
figPosition = [0.1178    0.2158    0.8648    0.6943];
global plotData
bkgSensIn = 0;
% global bkgSensIn

%% DEFINE POSSIBLE THRESHOLDS
if bkgSensIn == 0   % Normally, do the following...
    if PeakAnalysis.cPeak == 1    % If cPeak defined by residuals
      posThreshBins = min(Charcoal.peak):range(Charcoal.peak)/250:...
          max(Charcoal.peak);   
        % Possible threshold values.
    else                        % If cPeak defined by ratio
      posThreshBins = min(Charcoal.peak):range(Charcoal.peak)/250:...
          max(Charcoal.peak);
            % Possible threshold values.
    end
else    % ...but if background sensitivity test is underway, do this:
    if Results.allFigures == 1
    figure (2); xlimIn = get(gca,'xlim');
    posThreshBins = xlimIn(1):range(xlimIn)/250:xlimIn(2);   
            % Possible threshold values.
    end
end
CharThresh.possible = posThreshBins;

%% IF THRESHOLD IS USER-DEFINED, FORCE THRESHOLDS
if PeakAnalysis.threshMethod == 1 
    CharThresh.pos = ones(length(Charcoal.peak),...
        length(PeakAnalysis.threshValues));
    for i = 1:length(PeakAnalysis.threshValues)
    in1 = find(CharThresh.possible >= PeakAnalysis.threshValues(i), 1 );
    in2 = find(CharThresh.possible <= PeakAnalysis.threshValues(i), 1, 'last' );
        if abs(CharThresh.possible(in1)-PeakAnalysis.threshValues(i)) <...
                abs(CharThresh.possible(in2)-PeakAnalysis.threshValues(i))
            inFinal = in1;
        else
            inFinal = in2;
        end
        % Index to find the closest match between selected threshold, and 
        % threhsolds possible, given bins dividing up Charcaol.peak.
        CharThresh.pos(:,i) = CharThresh.possible(inFinal);
    end
    CharThresh.neg = -99*ones(length(Charcoal.peak),...
        length(PeakAnalysis.threshValues));
    CharThresh.noisePDF = -99;
end
%% IF THRESHOLD IS DATA-DEFINED, DEFINE NOISE DISTRIBUTION AND THRESHOLDS 
if PeakAnalysis.threshMethod > 1 % Use 0- or 1-mean normal distribution
    if PeakAnalysis.threshMethod == 2
        if PeakAnalysis.cPeak == 1    % If cPeak defined by residuals, then
                % use zero-mean normal distribution.
            sigmaHat = std([Charcoal.peak(Charcoal.peak<=0);...
                abs(Charcoal.peak(Charcoal.peak<0))]); 
                % Estimated standard deviation of noise distribution: based on 
                % Charcoal.peak values <= 0.
            muHat = 0;  % Assume mean = 0.
            CharThresh.noisePDF = normpdf(posThreshBins,muHat,sigmaHat);  
            % Noise PDF 
        end
        if PeakAnalysis.cPeak == 2;   % If cPeak defeind by ratio, then 
                % use one-mean normal distribution.
            sigmaHat = std([[Charcoal.peak(Charcoal.peak<=1)]-1;...
                abs((Charcoal.peak(Charcoal.peak<1)-1))]+1); 
                % Estimated standard deviation of noise distribution: based on 
                % Charcoal.peak values <= 1.
            muHat = 1;  % Assume mean = 1.
            CharThresh.noisePDF = normpdf(posThreshBins,muHat,sigmaHat);  
            % Noise PDF
        end
    end
    
    if PeakAnalysis.threshMethod == 3 % Use Gaussian Mixture Model to
                                        % estimate noise distribution
        [mu,sig] = GaussianMixture(Charcoal.peak,2,2,false); % Estimate mean  
            % and standard deviation of noise distribution using the CLUSTER 
            % mixture model:
            % http://cobweb.ecn.purdue.edu/~bouman/software/cluster/
        
        % Added, October, 2009, PEH: 
        if mu(1) == mu(2)
             beep
             disp('WARNING: poor fit of Gaussian mixture model;')
             disp('         re-fit starting with three classes')
             [mu,sig] = GaussianMixture(Charcoal.peak,3,2,false);
         end
        in = find(mu == min(mu), 1 ); % Index to find smallest mu
        sigmaHat = sig(in); % Standard deviation of noise distribution
        muHat = mu(in);    % Mean of noise distribution
        CharThresh.noisePDF = normpdf(posThreshBins,muHat,sigmaHat); % Noise 
            % PDF 
    end

    % Given mean and standard deviations, select threhsold 
    thresh = norminv([PeakAnalysis.threshValues],muHat,sigmaHat); 
            % upper threhsold values
        CharThresh.pos = ones(length(Charcoal.peak),...
            length(PeakAnalysis.threshValues));
        for i = 1:length(PeakAnalysis.threshValues)
            in1 = find(CharThresh.possible >= thresh(i), 1 );
            in2 = find(CharThresh.possible <= thresh(i), 1, 'last' );
            if abs(CharThresh.possible(in1)-thresh(i)) <...
                    abs(CharThresh.possible(in2)-PeakAnalysis.threshValues(i))
                inFinal = in1;
            else
                inFinal = in2;
            end
            % Index to find the closest match between selected threshold, and 
            % threhsolds possible, given bins dividing up Charcaol.peak.
            CharThresh.pos(:,i) = CharThresh.possible(inFinal);
        end
    threshNeg = norminv([1-PeakAnalysis.threshValues(3)],...
            muHat,sigmaHat);    % lower threshold value
        CharThresh.neg = threshNeg*ones(length(Charcoal.peak),1);
end    
%% DEFINE SIGNAL-TO-NOISE INDEX
signal = Charcoal.peak(Charcoal.peak > CharThresh.pos(4));    
    % Charcoal.peak values above threshold.
noise = Charcoal.peak(Charcoal.peak <= CharThresh.pos(4));    
    % Charcoal.peak values equal to or below threshold.
% NEW May 2010: SNI calculation based on Kelly et al., in review.
    %   - PEH
    if ~isempty(signal)
    CharThresh.SNI = (1/length(signal)) .* sum((signal - mean(noise))./ ...
        std(noise)) .* [(length(noise)-2)/length(noise)];
    else
    CharThresh.SNI = 0;
    end
    
%     Pre-May 2010 method for SNI calculation, based on Higuera et al. 2009
%     (Ecological Monogrphs, 79: 201-213). 
% CharThresh.SNI = var(signal) / [var(noise) + var(signal)];

%% CREATE DUMMY VARIABLE, CharThresh.GOF
CharThresh.GOF = -999*ones(size(Charcoal.peak));

%% PLOTTING, IF SELECTED
if plotData == 1 && Results.allFigures == 1
    figure (2); clf; set(gcf,'color','white','name',...
    'Peak CHAR distribution, threshold values, and noise distribution (if selected)',...
    'units','normalized','position',[figPosition])
    [n,x] = hist(Charcoal.peak,CharThresh.possible);
    h = bar(x,n/sum(n),1);
    set(h,'faceColor',[.5 .5 .5],'edgeColor',[.5 .5 .5])
    hold on
    if PeakAnalysis.threshMethod > 1
    plot(CharThresh.possible,CharThresh.noisePDF*...
        mean(diff(CharThresh.possible)),'k','lineWidth',2)
    end
    plot([CharThresh.pos(1,1) CharThresh.pos(1,1)],[0 max(n/sum(n))],'--k')
    plot([CharThresh.pos(1,4) CharThresh.pos(1,4)],[0 max(n/sum(n))],'r')
    plot([CharThresh.pos(1,2) CharThresh.pos(1,2)],[0 max(n/sum(n))],'--k')
    plot([CharThresh.pos(1,3) CharThresh.pos(1,3)],[0 max(n/sum(n))],'--k')
    plot([CharThresh.pos(1,4) CharThresh.pos(1,4)],[0 max(n/sum(n))],'r')
    set(gca,'tickdir','out','box','off')
    ylabel ('proportion or scaled density')
    xlabel ('peak CHAR (# cm^-^2 yr^-^1)')
    title ([char(site), ': ', num2str(zoneDiv(1)),' to ',...
        num2str(zoneDiv(end)),' cal. yr BP'])
    if PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak == 1
        legend ('peak CHAR distribution',...
        'estimated noise PDF, via 0-mean Gaussian',...
        'possible threshold values','selected threshold value')
    end
    if PeakAnalysis.threshMethod == 2 && PeakAnalysis.cPeak == 2
        legend ('peak CHAR distribution',...
        'estimated noise PDF, via 1-mean Gaussian',...
        'possible threshold values','selected threshold value')
    end
    if PeakAnalysis.threshMethod == 3 
        legend ('peak CHAR distribution',...
        'estimated noise PDF, via Gaussian mixture model',...
        'possible threshold values','selected threshold value')
    end
    if PeakAnalysis.threshMethod == 1
         legend ('peak CHAR distribution',...
        'possible threshold values','selected threshold value')
        text(1.5*CharThresh.pos(1,4),0.75*max(n/sum(n)),['Threshold = ',...
        num2str(CharThresh.pos(1,4))],'backgroundColor','w')
        text(1.75*CharThresh.pos(1,4),0.6*max(n/sum(n)),...
            ['signal-to-noise index = ',num2str(CharThresh.SNI)],...
            'backgroundColor','w')
    else
        text(1.5*CharThresh.pos(1,4),0.75*max(n/sum(n)),['Threshold = ',...
        num2str(CharThresh.pos(1,4)),'; ',...
        num2str(round(PeakAnalysis.threshValues(4)*100)),...
        '^t^h percentile'],'backgroundColor','w')
        text(1.75*CharThresh.pos(1,4),0.6*max(n/sum(n)),...
            ['signal-to-noise index = ',num2str(CharThresh.SNI)],...
            'backgroundColor','w')
    end
end

    