function [Charcoal Pretreatment gapIn] =...
    CharPretreatment(charData,site,Pretreatment,Results)
% CharPretreatment  Interpolate record, derrive CHAR, transform (if selec.).
%     [Charcoal Pretreatment gapIn] =...
%         CharPretreatment(charData,site,Pretreatment,Results)
%
% Interpolates charcoal data to resolution defined by yrInterp, derives raw
% and resampled charcoal accumulation rates, and log transforms charcoal
% accumulations rates, if desired. 

% Input variables:
% charData -- charcoal data in Charster format, without column header
% zoneDiv -- years to start and stop the record, and any years identifying
%       zone divisions within the record. e.g. zoneDiv = [-53 5000 10000];
% yrInterp -- years to interpolate record to.
% transform -- 0 for not transformations, 1 for base 10 log transform,
%                 and 2 for natural log transform.
% Output variables:
% Charcoal -- data stucture with 12 fields:
% Charcoal.cm: sample depths [cm]
% Charcoal.count: charcoal counts [#]
% Charcoal.vol: sample volume [cm^3]
% Charcoal.con: charcoal concentration [# cm^-3] 
% Charcoal.ybp: year of each charcoal samples [cal. ybp]
% Charcoal.ybpI: resampled year of each sample [cal. ybp]
% Charcoal.cmI: resampled sample depths [cm]
% Charcoal.countI: resampled charcoal counts [#]
% Charcoal.volI: resampled sample volume [cm^3]
% Charcoal.conI: resampled charcoal concentration [# cm^-3]
% Charcoal.acc: charcoal accumulation rate [# cm^-2 yr^-1]
% Charcoal.accI: resampled charcoal accumulation rate [# cm^-2 yr^-1]
                          
%% CREATE LOCAL VARIABLES
zoneDiv = Pretreatment.zoneDiv;
yrInterp = Pretreatment.yrInterp;
transform = Pretreatment.transform;
figPosition = [0.1343    0.2383    0.8648    0.6943];
global plotData

%% TRIM RECORD AT ybpStart AND ybpStop:
if length(zoneDiv) > 1  % Record will be trimmed if > 2 zoneDiv values 
ybpStop = zoneDiv(1);                % [cal ybp] year to stop record at. 
ybpStart = zoneDiv(length(zoneDiv)); % [cal ybp] year to start record at.
charData = charData(((charData(:,3)>= ybpStop)...
    &(charData(:,4)<= ybpStart)),:);% trim charData to reflect time 
                                    % interval from ybpStart to ybpStop.
end

%% SCREEEN RECORD FOR MISSING VALUES:
missingValuesIndex = find(charData(:,5)<=0); % Index for rows with sample 
    % volume <= 0. 
nMissingValues = length(missingValuesIndex); % Number of missing values.
if nMissingValues > 0  % if some levels were not sampled...
    startIn = missingValuesIndex...
        (find([99; diff(missingValuesIndex)] > 1)); % Index start of gaps
        % created by missing values.
    endIn = missingValuesIndex...
        (find([diff(missingValuesIndex)] > 1)); % Index end of gaps   
        % created by missing values. 
    endIn(end+1,1) = missingValuesIndex(end); % Add on last sample as last
        % end point for gaps.
    gapIn = [startIn endIn]; % Index values for start (j = 1) and end 
        % (j = 2) of each gap.
    
    nGaps = length(startIn);  % Number of gaps in the record  
    cmGaps = sum(charData(gapIn(:,2)+1,1)-charData(gapIn(:,1)-1,1)); % Sum 
        % of all cm of gap(s).
    yrGaps = sum(charData(gapIn(:,2)+1,3)-charData(gapIn(:,1)-1,3)); % Sum 
        % of all cm of gap(s).
    disp(['NOTE: ' num2str(nMissingValues) ' missing value(s); '...
            num2str(length(gapIn)) ' gap(s) totaling '...
        num2str(cmGaps) ' cm and ' num2str(round(yrGaps)) ' years.'])
    disp('      Values created via interpolation.')
    
    for i = 1:nGaps % Fill in gaps via interpolation.
        x = [charData(gapIn(i,1)-1,1) charData(gapIn(i,2)+1,1)]; % [cm] End 
            % depths for interpolation.
        cmInterp = (diff(charData(gapIn(i,1),1:2)));
        xi = [charData(gapIn(i,1)-1):cmInterp:charData(gapIn(i,2)+1,1)]; 
        % [cm] Desired interpolated depths. 
        y = [charData(gapIn(i,1)-1,5) charData(gapIn(i,2)+1,5)]; % [cm^3] 
            % End volume for interpolation.
        y2 = [charData(gapIn(i,1)-1,6) charData(gapIn(i,2)+1,6)]; % [pieces] 
            % End charcoal count for interpolation.    
        yi = interp1(x,y,xi);
        y2i = interp1(x,y2,xi);
        if length(yi(2:end-1)) ~= length(gapIn(i,1):gapIn(i,2))
            yi(2) = [];     % Trim inerpolated values if there are more 
            y2i(2) = [];    % interpolated values than needed, given 
                            % variable sampling intervals around gap. 
        end
        charData(gapIn(i,1):gapIn(i,2),5) = yi(2:end-1);    % Fill in 
            % interpolated sample volumes.
        charData(gapIn(i,1):gapIn(i,2),6) = y2i(2:end-1);   % Fill in 
            % interpolated sample counts.
    end
else
    gapIn = NaN(0,1);  % If no missing values, create empty 
        % variable to pass to CharPeakAnalysisResults.
end
remainingGap = find(charData(2:end,1)-charData(1:end-1,2));
remainingGapCm = charData(remainingGap+1,1)-charData(remainingGap,1);% [cm]
    % Total cm of remaining gap(s).
remainingGapYr = charData(remainingGap+1,3)-charData(remainingGap,3); % [yr]
    % Total yr of remaining gap(s).
if length(remainingGap) > 0
    disp(['NOTE: ' num2str(length(remainingGap))...
        ' gap(s) in record totaling '...
        num2str(sum(remainingGapCm)) ' cm and ' num2str(sum(remainingGapYr)) ' yr.'])
    disp('      Resampling will occur across this (these) gap(s); consider filling in values with')
    disp('      -999 to interpolate over gap.')
end

%% RETRIEVE VARIABLES FROM INPUT FILES:
Charcoal.cm = charData(:,1);    % [cm] sample depths.
Charcoal.count = charData(:,6); % [#] charcoal counts
Charcoal.vol = charData(:,5);   % [cm^3] sample volumes.
Charcoal.con = Charcoal.count ./ Charcoal.vol; % [# cm-3] charcoal con.
Charcoal.ybp = charData(:,3);   % [cal ybp] age at top of sample

%% DERIVE SEDIMENT ACCUMULATOIN RATE:
sedAcc = zeros(length(Charcoal.cm),1); % space for sediment acc. rate          
for i = 1:length(Charcoal.ybp)-1
    sedAcc(i) = (Charcoal.cm(i+1)-Charcoal.cm(i))./...
        (Charcoal.ybp(i+1)-Charcoal.ybp(i)); 
    % [cm/yr] sediment accumulation rate
end

%% CALCULATE yrInterp IF NOT USER-DEFINED
if Pretreatment.yrInterp == 0
    yrInterp = round(median(diff(charData(:,3))));
    Pretreatment.yrInterp = yrInterp;
    disp('                                             ')
    disp(['Record resampled to median resolution of '...
        num2str(Pretreatment.yrInterp) ' years.'])
end

%% INTERPOLATE CHAROCAL DATA TO yrInterp INTERVALS:
cmTop = charData(:,1);  % [cm] Depth at top of sample.
cmBot = charData(:,2);  % [cm] Depth at bottom of sample.
ageTop = charData(:,3); % [yr BP] Age at top of sample.
ageBot = charData(:,4); % [yr BP] Age at bottom of sample.

Charcoal.ybpI = [Pretreatment.zoneDiv(1):Pretreatment.yrInterp:...
    Pretreatment.zoneDiv(end)]'; % [yr BP] Years to 
    % resample record to.

propMatrix = NaN*ones(length(Charcoal.ybpI),length(Charcoal.ybp));

%  Determine the proportion of each sample contained in resampled sample
    % figure (100); clf; set(gcf,'color','w')
    %         plot([1 2],[ageTop(j) ageTop(j)],'k','linewidth',2); hold on
    %         plot([3 4],[rsAgeTop rsAgeTop],'b','linewidth',2); 
    %         plot([1 2],[ageBot(j) ageBot(j)],'k','linewidth',2); 
    %         plot([3 4],[rsAgeBot rsAgeBot],'b','linewidth',2); 
    %         set(gca,'box','off','ydir','rev'); grid on
    %         ylim ([min([ageTop(j) ageBot(j) rsAgeTop rsAgeBot]-10),...
    %             10+max([ageTop(j) ageBot(j) rsAgeTop rsAgeBot])])
    %         legend ('record','re-sample')
    %         text(1.5,ageTop(j),'ageTop','backgroundcolor','w')
    %         text(1.5,ageBot(j),'ageBot','backgroundcolor','w')
    %         text(3.5,rsAgeTop,'rsAgeTop','backgroundcolor','w')
    %         text(3.5,rsAgeBot,'rsAgeBot','backgroundcolor','w')

for i = 1:length(Charcoal.ybpI)  % For each resampled sample.
    rsAgeTop = Charcoal.ybpI(i);
    rsAgeBot = rsAgeTop+Pretreatment.yrInterp;
    
    for j = 1:length(Charcoal.ybp)  % For each raw sample.
       % If raw sample straddles rsAgeBot
       if ageTop(j) >= rsAgeTop && ageTop(j) < rsAgeBot &&...
               ageBot(j) > rsAgeBot
           propMatrix(i,j) = rsAgeBot - ageTop(j);
       end
       % If raw sample straddles rsAgeTop
       if ageTop(j) < rsAgeTop && ageBot(j) <= rsAgeBot &&...
               ageBot(j) > rsAgeTop
           propMatrix(i,j) = ageBot(j) - rsAgeTop;
       end
       % If raw sample is entirely within resampled sample.
       if ageTop(j) >= rsAgeTop && ageBot(j) <= rsAgeBot
           propMatrix(i,j) = ageBot(j) - ageTop(j);
       end
       % If raw sample is entirely outside of resamples sample (i.e.
       % resampling finer than actual sample).
       if ageTop(j) < rsAgeTop && ageBot(j) > rsAgeBot
           propMatrix(i,j) = rsAgeBot - rsAgeTop;
       end
    end
end % End making porportion matrix
propMatrix = propMatrix ./ Pretreatment.yrInterp;

% Determine values for each resampled interval.
Charcoal.cmI = NaN*ones(length(Charcoal.ybpI),1);
Charcoal.countI = NaN*ones(length(Charcoal.ybpI),1);
Charcoal.volI = NaN*ones(length(Charcoal.ybpI),1);
Charcoal.conI = NaN*ones(length(Charcoal.ybpI),1);
sedAccI = NaN*ones(length(Charcoal.ybpI),1);

for i = 1:length(Charcoal.ybpI)  % For each resampled sample
    in = find (propMatrix(i,:)>0); % Index for raw samples contributing to 
        % resampled sample.
    
%     Charcoal.cmI(i) = sum(Charcoal.cm(in) .* propMatrix(i,in)');
    Charcoal.countI(i) = sum(Charcoal.count(in) .* propMatrix(i,in)');
    Charcoal.volI(i) = sum(Charcoal.vol(in) .* propMatrix(i,in)');
    Charcoal.conI(i) = sum(Charcoal.con(in) .* propMatrix(i,in)');
    sedAccI(i) = sum(sedAcc(in) .* propMatrix(i,in)'); 
end
    Charcoal.cmI = interp1(Charcoal.ybp,Charcoal.cm,Charcoal.ybpI); % [cm] 
    % interpolated depths.
    
%% DERIVE CHARCOAL ACCUMULATOIN RATE:
Charcoal.acc = (Charcoal.con.*sedAcc)'; % [#/cm2/yr] take sedAcc and 
    % multiply by Charcoal.con to get  Charcoal.acc
Charcoal.accI = (Charcoal.conI .* sedAccI); % [#/cm2/yr] take sedAccI and 
    % multiply by Charcoal.conI to get Charcoal.accI   
    
%% TRANSFORM DATA, IF DESIRED:
if Pretreatment.transform == 1
    Charcoal.accI = log10(Charcoal.accI+1);
end
if Pretreatment.transform == 2
    Charcoal.accI = log(Charcoal.accI+1);
end
    
%%  PLOT RAW CHAR AND RESAMPLED SERIES, IF DESIRED
if plotData == 1 & Results.allFigures == 1
    figure (1); clf; set(gcf,'name','(a) C_raw and C_resampled; (b) C_raw and different options for C_background',...
        'units','normalized','position',[figPosition],'color','w')
    subplot(2,1,1)
    h = bar(Charcoal.ybp,Charcoal.acc,1);
    set(h,'facecolor',[.5 .5 .5],'edgecolor',[.5 .5 .5])
    ylim ([0, prctile(Charcoal.acc,99)]);
    hold on
    stairs(Charcoal.ybpI-0.5*Pretreatment.yrInterp,Charcoal.accI,...
        '.-k','linewidth',1)
   if length(gapIn) > 0    % If there are gaps in the record, then
        % fill in the graph with white
        buffer = 0;%PeakAnalysis.peakFrequ/2; % [yr] Years to pad the years 
            % with missing values. 
        for i = 1:nGaps
            y = get(gca,'ylim');
            x = [Charcoal.ybp(gapIn(i,1))-buffer ...
                Charcoal.ybp(gapIn(i,1))-buffer ...
                Charcoal.ybp(gapIn(i,2))+buffer ...
                Charcoal.ybp(gapIn(i,2))+buffer];
            patch(x,[y(1) y(2) y(2) y(1)],'w','edgecolor',[.75 .75 .75],...
                'facecolor',[.75 .75 .75])
        end
         legend('C_r_a_w',['C_i_n_t_e_r_p_o_l_a_t_e_d: '...
        num2str(yrInterp) ' yr'],'missing values') 
   else
        legend('C_r_a_w',['C_i_n_t_e_r_p_o_l_a_t_e_d: '...
        num2str(yrInterp) ' yr']) 
   end
    xlabel ('time (cal. yr BP)')
    ylabel ('CHAR (# cm^-^2 yr^-^1)')
    set(gca,'xlim',[min(Charcoal.ybp)-100 max(Charcoal.ybp)],'box','off',...
        'tickdir','out','xdir','rev')
    title ([{char(site)}, {'(a) C_r_a_w and C_i_n_t_e_r_p_o_l_a_t_e_d'}])
end