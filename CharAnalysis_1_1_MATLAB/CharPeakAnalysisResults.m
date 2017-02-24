function [Charcoal] = CharPeakAnalysisResults (Charcoal, Pretreatment,...
    PeakAnalysis, CharThresh, Smoothing, site, Results, fileName, gapIn)
% CharPeakAnalysisResults   Summarize and plot results from CharAnalysis.m
%   [Charcoal] = CharPeakAnalysisResults (Charcoal, Pretreatment,...
%    PeakAnalysis, CharThresh, Smoothing, site, Results, fileName, gapIn)
%
%    Plots 9 or 10 figures, depending on parameter choices, illustrating
%    results of CharAnalysis. Data and figures are saved, if desired.

%% CREATE LOCAL VARIABLES
zoneDiv = Pretreatment.zoneDiv;
r = Pretreatment.yrInterp;  % Resolution of resampled record.
yr = Smoothing.yr;          % Smoothing window for low-frequency CHAR.
transform = Pretreatment.transform;
[temp1 temp2] = size(gapIn);% Temporary variable to hold size of gapIn.
nGaps = temp1;              % Number of gaps in record.

if max (Charcoal.ybp)>5500
wm = 1.25*2.8e-005; % time factor, for scaling width of graphs (make 
                    % bigger to make graphs wider)
else
    if max (Charcoal.ybp)>2500
        wm = 1.5*2.8e-005;      
    else
        wm = 5*2.8e-005;
    end
end
if max (Charcoal.ybp)>15000
wm = 0.5*2.8e-005; % time factor, for scaling width of graphs (make 
                    % bigger to make graphs wider)
end

FS = 8;             % font size for tick labels
FW = 'bold';        % font weight 
zoneText = {'Zone 1' 'Zone 2' 'Zone 3' 'Zone 4' 'Zone 5' 'Zone 6' 'Zone 7'};
figPosition = [0.1013    0.1933   0.8648    0.6943];

%% FIGURE 3
figure (3); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Resampled charocal (C_interpolated), low-frequency trends (C_background), and detrended series (C_peak)'); 
subplot(3,5,1:4)    % C_raw and C_background
    x = Charcoal.ybpI;
    y = Charcoal.accI;
    y2 = Charcoal.accIS;
    H = bar(x,y); hold on;
    set(H,'facecolor',[0 0 0],'edgecolor',[0 0 0])
    plot (x,y2,'color',[.5 .5 .5],'linewidth',2); 
    xlim ([zoneDiv(1), zoneDiv(length(zoneDiv))]);
    ylim ([0, 1.1*max(y)]);
%     ylim ([0, prctile(y,99)]);
    if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv],[max(y)*1.01 max(y)*1.1],'-k','color',...
        [.5 .5 .5],'linewidth',2)
    end
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out','XTick',[0:1000:max(zoneDiv)],...
        'XTickLabel',[],'Position',[0.1 .55 wm*range(zoneDiv) .175])
    box off
    if transform == 0
            ylabel ('CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 1
            ylabel ('log CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 2
            ylabel ('ln CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(x)*1.01 max(x)*1.1],'-k',...
                'color',[.5 .5 .5],'linewidth',2); 
        end
        if length(zoneDiv) > 2
            text(mean(zoneDiv(z-1:z)),[max(y)*1.05],zoneText(z-1),...
            'horizontalalignment','center','fontweight','normal',...
            'fontsize',FS)
        end
    end
    title([{['(a) ',char(site),': C_i_n_t_e_r_p_o_l_a_t_e_d (',num2str(r),...
        ' yr)',' and C_b_a_c_k_g_r_o_u_n_d defined by ',...
        num2str(yr(length(yr))),'-yr trends']}],'FontSize',FS,...
        'FontWeight',FW)
    
subplot (3,5,10);   % pCHAR and threshold values.
    y = Charcoal.peak;
    y2 = CharThresh.pos(:,end);
    y3 = CharThresh.neg(:,end);
    if PeakAnalysis.cPeak == 1  % If cPeak defined by residuals.
        [H1] = bar(x,y,1);
    else                        % If cPeak defined by ratios.
        [H1] = bar(x,y,1,'BaseValue',1);
    end
    set(H1,'facecolor',[0 0 0],'edgecolor',[0 0 0])
    hold on 
    plot (x,y2,'r')
    plot(x,y3,'r')
    if PeakAnalysis.cPeak == 1 % If detrending by residuals
        plot([min(x) max(x)],[0 0],'k');
    else
        if PeakAnalysis.cPeak == 2 % If detrending by ratio
        plot([min(x) max(x)],[1 1],'k');
        end
        if PeakAnalysis.cPeak == -99 % If not using pCHAR
            plot(x,thresholds(:,1),'color',[.5 .5 .5],'linewidth',2);
            % plot median 
        end
    end
    if PeakAnalysis.threshType == 1    % If global threshold
        threshIn1 = sum(Charcoal.charPeaksThresh) ./...
            sum(Charcoal.charPeaks);    % Index for positive threshold
            % values.
        threshIn(1) = min(find(threshIn1 >= CharThresh.pos(1,1)));
        threshIn(2) = min(find(threshIn1 >= CharThresh.pos(1,2)));
        threshIn(3) = min(find(threshIn1 >= CharThresh.pos(1,3)));
        threshIn(4) = min(find(threshIn1 >= CharThresh.pos(1,4)));
        peakIn = find(Charcoal.charPeaks(:,threshIn(end)));
        CharcoalCharPeaks = Charcoal.charPeaks(:,threshIn);
        
        peakScreenIn = find(CharThresh.minCountP(:,threshIn(end)) >...
        PeakAnalysis.minCountP);
        
    else
        peakScreenIn = find(CharThresh.minCountP(:,end) >...
        PeakAnalysis.minCountP);
        peakIn = find(Charcoal.charPeaks(:,end));
        CharcoalCharPeaks = Charcoal.charPeaks;
    end
    
%     plot (x(peakScreenIn),logical(peakScreenIn)*0.8*prctile(y,99),'.k',...
    plot (x(peakScreenIn),logical(peakScreenIn)*0.8*max(y),'.k',...
        'color',[.75 .75 .75]);
    plot (x(peakIn),0.8*max(y),'+k')
%     plot (x(peakIn),[0.8*prctile(y,99)],'+k')
    xlim ([zoneDiv(1), zoneDiv(length(zoneDiv))]);
    ylim ([min(y3), 1.1*max(y)]);
%     ylim ([min(y3), prctile(y,99)]);
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick',...
        'on','Position',[0.1 .32 wm*range(zoneDiv) .175],'TickDir','out',...
        'XTick',[0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000])
    if PeakAnalysis.cPeak == 1
        title([{'(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d - C_b_a_c_k_g_r_o_u_n_d), thresholds defining C_n_o_i_s_e, and peaks',...
        'identified (gray dots fail to pass peak-magnitude test)'}],...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
    else
        if PeakAnalysis.cPeak == 2
            title(['(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d / C_b_a_c_k_g_r_o_u_n_d), ',...
        'thresholds defining C_n_o_i_s_e, and peaks identified',...
        ' (gray peaks fail to pass peak-magnitude test)'],...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
        end
        if PeakAnalysis.cPeak == -99
        title(['(b) C_i_n_t_e_r_p_o_l_a_t_e_d, ',...
        'thresholds defining C_n_o_i_s_e, and peaks identified',...
        ' (gray peaks fail to pass peak-magnitude test)'],...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
        end
    end
    xlabel ('time (cal. yr BP x 1000)','FontSize',FS)
    if transform == 0
        ylabel ('CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 1
        ylabel ('log CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 2
        ylabel ('ln CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    box off

%% FIGURE 4
figPosition = figPosition - [0.0165 0.0225 0 0];
figure (4); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Sensitvity to alternative threhsolds and quality of record'); 
subplot(3,5,1:4)    % C_raw, C_background, and threshold value(s)
    x = Charcoal.ybpI;    
    y = Charcoal.accI;
    y2 = Charcoal.accIS;
    y3 = CharThresh.pos(:,end);
    y4 = CharThresh.neg(:,end);
    y5 = CharcoalCharPeaks;
    y5(y5 == 0) = -99;
    H = bar(x,y,1); hold on;
    set(H,'facecolor',[0 0 0],'edgecolor',[0 0 0])
    plot (x,y2,'color',[.5 .5 .5],'linewidth',2); 
    xlim ([zoneDiv(1), zoneDiv(length(zoneDiv))]);
    ylim ([0, 1.15*max(y)]);
%     ylim ([0,prctile(y,99)]);
    if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv],[max(y)*1.01 max(y)*1.1],'-k','color',...
        [.5 .5 .5]);
    end
    hold on
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick',...
        'on','TickDir','out','XTick',[0:1000:max(zoneDiv)],'XTickLabel',...
        [0:max(zoneDiv)/1000],'Position',[0.1 .75 wm*range(zoneDiv) .175])
        if transform == 0
            ylabel ('CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
        end
        if transform == 1
            ylabel ('log CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
        end
        if transform == 2
            ylabel ('ln CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
        end
    if PeakAnalysis.cPeak == 1  % If using residuals.
        plot (x,y2+y3,'r')
        plot(x,y2+y4,'r')
    else                        % Otherwise, using ratios.
        plot(x,y2.*y3,'r')
        plot(x,y2.*y4,'r')
    end
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(y)*1.01 max(y)*1.1],'-k',...
                'color',[.5 .5 .5],'linewidth',2); 
        end
        if length(zoneDiv) > 2
        text(mean(zoneDiv(z-1:z)),[max(y)*1.05],char(zoneText(z-1)),...
        'horizontalalignment','center','fontweight','normal','fontsize',FS)
        end
    end
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick',...
    'on','TickDir','out','XTick',[0:1000:max(zoneDiv)],'XTickLabel',[],...
    'Position',[0.1 .55 wm*range(zoneDiv) .175])
    plot(x,max(y)*0.78*y5(:,1),'.k','color',[.5 .5 .5])
    plot(x,max(y)*0.85*y5(:,2),'.k','color',[.5 .5 .5])
    plot(x,max(y)*0.92*y5(:,3),'.k','color',[.5 .5 .5])

    if sum(y5(:,1) - y5(:,end)) == 0; mIndex = 0.78; end
    if sum(y5(:,3) - y5(:,end)) == 0; mIndex = 0.92; end
    if sum(y5(:,2) - y5(:,end)) == 0; mIndex = 0.85; end
    plot(x,max(y)*mIndex*y5(:,end),'.k','color',[1 1 1])
    plot(x,max(y)*mIndex*y5(:,end),'+k')
    if transform == 0
        ylabel ('CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 1
        ylabel ('log CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if transform == 2
        ylabel ('ln CHAR (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    title([{char(site)},{''},{'(a) C_i_n_t_e_r_p_o_l_a_t_e_d, C_b_a_c_k_g_r_o_u_n_d, and peak ID (+)'}],...
        'FontSize',FS,'FontWeight',FW)
    box off

subplot(3,5,5)  % Sensitivity to alternative thresholds
if PeakAnalysis.threshType == 1   % If threshold is defined globally.
    if PeakAnalysis.cPeak == 1    % If cPeak defined by residuals.
        x = CharThresh.possible;  % Potential threhsolds.
        x2 = x(x>0);
        y = Charcoal.peak;              % Peak CHAR.
        y3 = sum(Charcoal.charPeaks);   % Total peaks for each threshold.
    else                          % else, cPeak defined by ratios.
        x = CharThresh.possible;  % Potential threhsolds.
        x2 = x(x>1);
        y = Charcoal.peak;              % Peak CHAR.
        y3 = sum(Charcoal.charPeaks);   % Total peaks for each threshold.
    end
    [n,x] = hist(y,x);
    [H1] = bar(x,n/sum(n));
    set(H1,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75],'BarWidth',1.0)
    hold on
    % Plot Noise PDF and sensitivity of peaks identified.
    if PeakAnalysis.threshMethod > 1  % If noise PDF used.
        if PeakAnalysis.cPeak == 1    % If cPeak defined by residuals.
        y2 = CharThresh.noisePDF;       % Noise PDF.
        [AX,H1,H2] = plotyy(x,y2*mean(diff(x)),...
        x(x>0),y3,'plot','plot'); 
        else                            % Else, cPeak defined by ratios.
            xIn = length(x)-length(x(x>=1));
            xIn2 = xIn - [length(x)-length(y3)];
        y2 = CharThresh.noisePDF;       % Noise PDF.
        [AX,H1,H2] = plotyy(x,y2*mean(diff(x)),...
        x(xIn:end),y3(xIn2:end),'plot','plot'); 
        end
    else                               % Else, thresholds are user defined.
        [AX,H1,H2] = plotyy(-99,-99,...
        x(x>0),y3,'plot','plot');  
    end
    
    % Plot peaks for all threshold cirteria.
    for t = 1:4              
        plot([CharThresh.pos(1,t) CharThresh.pos(1,t)],...
        [0 max(n/sum(n))],'--','color',[0.75 0.75 0.75]); 
        if t == 4
            plot([CharThresh.pos(1,t) CharThresh.pos(1,t)],...
            [0 max(n/sum(n))],'-k');
            text(CharThresh.pos(1,t),0,['<'],...
            'FontSize',FS+2,'Rotation',90,'FontWeight','bold')
        end
    end
    % Label threshold value, and percentile cut-off if relevant.
    if PeakAnalysis.threshMethod > 1  % If noise PDF used
        text(CharThresh.pos(1,end)*1.25,0.9*max(n/sum(n)),...
        [{[num2str(PeakAnalysis.threshValues(end)*100),...
        '^t^h percentile']},...
        {['= ',num2str(round(max(x(x<=CharThresh.pos(1,end)))*...
        10000)/10000)]}],'FontSize',FS,'HorizontalAlignment','left',...
        'backgroundcolor','w')
    else
    text(CharThresh.pos(1,end)*1.25,0.9*max(n/sum(n)),...
        [{['threshold value']},...
        {['= ',num2str(round(max(x(x<=CharThresh.pos(1,end)))*...
        10000)/10000)]}],...
        'FontSize',FS,'HorizontalAlignment','left','backgroundcolor','w')
    end
    % Set figure properties
    set(AX(1),'TickDir','out','FontSize',FS,'Ycolor','k','YMinorTick','off',...
        'XLim',[min(Charcoal.peak) 0.75*max(Charcoal.peak)],'YLim',...
        [0 1.01*max([max(n/sum(n)), max(n/sum(n))])],...
        'YTick',[0:round(10*max(n/sum(n)))/10/4:max(n/sum(n))],'position',...
        [0.17+wm*range(zoneDiv) .55 .16 .175],'XMinorTick','on');
    set(AX(2),'TickDir','out','FontSize',FS,'Ycolor','k','YTick',...
        [0:round(y3/40)*10:0.9*...
        max(y3)],'YMinorTick','on','XLim',...
        [min(Charcoal.peak) 0.75*max(Charcoal.peak)],'YLim',...
        [0 0.99*max(y3)],...
        'position',[0.17+wm*range(zoneDiv) .55 .16 .175]);
    set(H1,'color',[0 0 0],'linestyle','--','linewidth',1.5)
    set(H2,'color','k','linestyle','-','linewidth',1.5)
    set(get(AX(1),'Ylabel'),'String','relative frequency','fontsize',FS)
    set(get(AX(2),'Ylabel'),'String','# of peaks identified','fontsize',FS,...
        'rotation',270,'verticalalignment','bottom','color','k')
    ylabel ('relative frequency','FontSize',FS)
    if PeakAnalysis.cPeak == 1 
        xlabel ('residual CHAR value (pieces cm^-^2 yr^-^1)','FontSize',FS)
    else
        xlabel ('CHAR ratio (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    box off
    if PeakAnalysis.threshMethod > 1
        title ([{'(b) High-frequency CHAR distribution,'},...
        {'est. noise dist., and threshold values'}],'FontWeight',FW)
    else
        title ([{'(b) High-frequency CHAR distribution'},...
        {'and threshold values'}],'FontWeight',FW)
    end
end   % End if threshold defined globally.
if PeakAnalysis.threshType == 2   % If threshold is defined locally.
    plot_in = fliplr(1:length(zoneDiv)-1);  
    for i = 1:length(plot_in)
    xPlot1 = Charcoal.ybpI(find(CharcoalCharPeaks(:,1)>0)); % yr 
            % with peaks, using lower threshold.
    xPlot2 = Charcoal.ybpI(find(CharcoalCharPeaks(:,2)>0)); % middle
    xPlot3 = Charcoal.ybpI(find(CharcoalCharPeaks(:,3)>0)); % upper 

    xPlot1 = xPlot1(xPlot1 >= zoneDiv(i) & xPlot1 < zoneDiv(i+1)); % select
        % only years between zoneDiv(i) and zoneDiv(i+1) -- veg. zones.
    xPlot2 = xPlot2(xPlot2 >= zoneDiv(i) & xPlot2 < zoneDiv(i+1)); 
    xPlot3 = xPlot3(xPlot3 >= zoneDiv(i) & xPlot3 < zoneDiv(i+1)); 
    
    FRI_thresh = NaN*ones(200,3);                
    FRI_thresh(1:length(diff(xPlot1)),1) = diff(xPlot1)';
    FRI_thresh(1:length(diff(xPlot2)),2) = diff(xPlot2)';
    FRI_thresh(1:length(diff(xPlot3)),3) = diff(xPlot3)';

        for j = 1:3 % for each threshold
            mFRI_thresh(j,1) = mean(FRI_thresh(FRI_thresh(:,j)>0,j));
            if sum(FRI_thresh(:,j)>0) > 1
                mFRI_boot = bootstrp(1000,'mean',FRI_thresh(FRI_thresh(:,j)>0,j));
                mFRI_ci(j,1:2) = prctile(mFRI_boot,[2.5 97.5]);% 95% CI for mFRI
            else
                mFRI_ci(j,1:2) = [NaN NaN];
            end
        end
        x = [plot_in(i)-0.25 plot_in(i) plot_in(i)+0.25];
        errorbar(x,mFRI_thresh,mFRI_thresh-mFRI_ci(:,1),...
            mFRI_ci(:,2)-mFRI_thresh,'.k','color',[.5 .5 .5]);
        hold on
        in = min(find(PeakAnalysis.threshValues ==...
            PeakAnalysis.threshValues(end)));
        errorbar(x(in),mFRI_thresh(in),mFRI_thresh(in)-mFRI_ci(in,1),...
            mFRI_ci(in,2)-mFRI_thresh(in),'+k');
    end 
    set(gca,'TickDir','out','FontSize',FS,'xtick',[1:length(zoneDiv)-1],...
        'xticklabel',[fliplr(1:length(zoneDiv)-1)],'position',...
        [0.17+wm*range(zoneDiv) .55 .10 .175]);
        box off
        xlabel('zone')
        ylabel([{'zone-specific mean FRI'};{'+/- 95% ci (years fire^-^1)'}])
        title ([{'(b) Sensitivity to'}; {'alternative thresholds'}],...
                'fontweight','bold')
end     % End for each zone

subplot(3,5,11:13)  % Local signal-to-noise index
    if PeakAnalysis.threshType == 1
        CharThresh.SNI = CharThresh.SNI.*ones(size(Charcoal.peak));
    end
    semilogy(Charcoal.ybpI,CharThresh.SNI,'k')
    plot(Charcoal.ybpI,CharThresh.SNI,'k')
    y_lim = [0 10];
    xlim ([zoneDiv(1), zoneDiv(length(zoneDiv))]);
    ylim (y_lim);
    if y_lim(2) < 20
         y_tick = [0:2:y_lim(2)];
    end
    if y_lim(2) >= 20
         y_tick = [0:5:y_lim(2)];
    end
    if y_lim(2) >= 50
         y_tick = [0:10:y_lim(2)];
    end
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
    'YMinorTick','on','Position',[0.1 .35 wm*range(zoneDiv) .12],'TickDir',...
    'out','XTick',[0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000],...
    'YTick',y_tick)
    hold on
    plot([zoneDiv(1), zoneDiv(length(zoneDiv))],[3 3],'k--')
    box off
    xlabel ('time (cal. yr BP x 1000)','FontSize',FS)
    ylabel ('signal-to-noise index','FontSize',FS)
    title('(c) Local signal-to-noise index',...
        'FontSize',FS,'FontWeight',FW)
   
subplot(3,5,14)     % Global signal-to-noise index
    boxplot(CharThresh.SNI,'colors','kkk',...
        'symbol','.k')
    set(gca,'TickDir','out','FontSize',FS,...
    'YMinorTick','off','position',[0.17+wm*range(zoneDiv) .35 .10 .12],...
    'xtick',[1],'xticklabel','','ylim',y_lim,'YTick',...
    y_tick,'yscale','linear')
    box off; grid off
    hold on
    plot([zoneDiv(1), zoneDiv(length(zoneDiv))],[3 3],'k--')
    title('(d) Global signal-to-noise index','FontSize',FS,'FontWeight',FW)
    xlabel([{'global signal-to-'}; {'noise distribution'}],'FontSize',FS)
    ylabel ('')
        sig_text = [num2str((nanmedian(CharThresh.SNI)))];
    text(1.2,median(CharThresh.SNI),sig_text,...
        'FontSize',FS,'BackgroundColor','w')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% FIGURE 5
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(5); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Cumulative peaks through time')    % Cumulative fires through time.
    xPlot = Charcoal.ybpI(find(CharcoalCharPeaks(:,end)>0)); % Years with 
        % fires.
    yPlot = flipud(cumsum(logical(find(CharcoalCharPeaks(:,end)>0)))); 
        % Cumulative fires through time.
    if length (xPlot) > 0
    plot(xPlot,yPlot,'.k')
    set (gca,'XDir','reverse','TickDir','out','XMinorTick','on',...
        'YMinorTick','on','FontSize',FS,'XTick',...
        [0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000],...
        'Position',[.1 .1 0.8*wm*range(zoneDiv) max(yPlot)*0.005]);
    if nGaps > 0    % If there are gaps in the record, then
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
        legend('Cumulative peaks','missing values') 
    end
    xlabel ('time (cal. yr BP x 1000)','FontSize',FS);
    ylabel ('cumulative number of peaks','FontSize',FS);
    xlim ([min(Charcoal.ybp),max(zoneDiv)])
    title ([char(site),': Cumulative fires as a function of time'],...
        'FontSize',FS,'FontWeight','bold')
    box off
    grid on   
    end
    
%% FIGURE 6
figPosition = figPosition - [0.0165 0.0225 0 0];
figure (6);clf;set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Fire return intervals by zone, with Weibull models if GOF test is passed')
    % FRI distributions by zone
FRI_params_zone = -999*ones(length(zoneDiv)-1,10);  % Space for FRI results: 
    % Nfires, mFRI, mFRI_uCI, mFRI_lCI, WblB, WblB_uCI, WblB_lCI, WblC,
    % WblC_uCI, WblC_lCI.
for i = 1:length(zoneDiv)-1
    inPlot = fliplr([1:length(zoneDiv)-1]);
    subplot(3,5,inPlot(i))
    binWidth = 20;    % [yr] Width for FRI historgram bins.
    xPlot = Charcoal.ybpI(find(CharcoalCharPeaks(:,end)>0)); % Years with 
        % fires.
    xPlot = xPlot(xPlot >= zoneDiv(i) & xPlot < zoneDiv(i+1)); % Select 
        % only years between zoneDiv(i) and zoneDiv(i+1) -- veg. zones
    FRI = diff(xPlot);    % [yr] Fire return intervals.
    if max(FRI) > 5000
        disp(['WARNING, Figure 6: FRIs in Zone ' num2str(i) ' > 5000 years.'])
        disp('     FRI distribution not characterized.')
        text(0.25, 0.5, 'FRIs > 5000 yr')
    else
    y_lim = ([0 0.35]);
    x_lim = ([0 800]);
    if length(FRI) > 1
    [FRI_freq FRI_bin] = hist(FRI,[binWidth:binWidth:1000]); 
            %[yr] fire return interval bins and corresponding frequencies
    param = wblfit(FRI_bin,[],[],FRI_freq); % Weibull a (yr) and b 
                                            % (unitless) parameters.
    % One-sample K-S Goodness-of-fit test:
    FRIBinKS = [0:20:5000]; % Bins for FRIs for KS test.
        % Add a warning if FRIs are longer than FRIBinKS
    wbl_cdf = wblcdf(FRIBinKS,param(1),param(2));    % Estimated wbl cdf.
    [h,pKS,t] = kstest(FRI,[FRIBinKS', wbl_cdf']); % results of test
    % Boot-strapped confidenceintervals on wbl parameters and mFRI:
    mean_mFRI = NaN*ones(1000,1);   % Space for mean_mFRI
    for t = 1:1000
       FRI_t = randsample(FRI,length(FRI),true); % random sample, with 
        % replacement of FRI data
       [FRI_freq_t FRI_bin] = hist(FRI_t,FRI_bin);% FRI_t distribution
       param_t(t,:) = wblfit(FRI_bin,[],[],FRI_freq_t);% wbl params for FRI_t   
       mean_mFRI(t) = mean(FRI_t); % median FRI for this trial
    end
    wbl_a_ci = prctile(param_t(:,1),[2.5 97.5]);% 95% CI for wbl a parameter
    wbl_b_ci = prctile(param_t(:,2),[2.5 97.5]);% 95% CI for wbl b parameter
    mean_mFRI_ci = prctile(mean_mFRI,[2.5 97.5]);   % 95% CI for mFRI
    FRI_params_zone(i,:) = [length(FRI) mean(FRI) mean_mFRI_ci param(1)...
        wbl_a_ci param(2) wbl_b_ci];   
        % store FRI parameters for each zone: wbl a, ci, ci, wbl b, ci, ci,
        % mFRI, ci, ci
    
    wbl_est = binWidth*wblpdf(1:1000,param(1),param(2));   % Multiplied by 
            % binWidth to make proportional         
    [H] = bar(FRI_bin,FRI_freq/sum(FRI_freq));
    set(H,'FaceColor',[.75 .75 .75],'BarWidth',1.0)
    hold on
    if length(FRI)>4    % Only proceed if there are more than 4 FRIs
    if  [length(FRI) < 30 & pKS > 0.10] | [length(FRI) >= 30 &...
                pKS > 0.05]

    plot(1:1000,wbl_est,'k','linewidth',2)
    text (max(x_lim),max(y_lim),[{['Wbl\it b \rm = ',...
        num2str(round(param(1))),' (',num2str(round(wbl_a_ci(1))),...
        '-',num2str(round(wbl_a_ci(2))),')']},...
        {['Wbl\it c \rm = ',num2str(round(param(2)*100)/100),' (',...
        num2str(round(wbl_b_ci(1)*100)/100),'-',...
        num2str(round(wbl_b_ci(2)*100)/100),')']},...
        {['mFRI = ',num2str(round(mean(FRI))),' (',...
        num2str(round(mean_mFRI_ci(1))),'-',...
        num2str(round(mean_mFRI_ci(2))),')']},...
        {['N_F_R_I = ',num2str(length(xPlot)-1)]}],...
        'HorizontalAlignment','right','VerticalAlignment','top',...
        'FontSize',FS)
    end
    else
        text(max(x_lim),0.8*max(y_lim),['N_F_R_I = ',...
            num2str(length(xPlot)-1)],'HorizontalAlignment','right',...
            'FontSize',FS)
    end
    else
        if length(xPlot) < 1
            text(max(x_lim),0.8*max(y_lim),['N_F_R_I = 0'],...
                'HorizontalAlignment','right','FontSize',FS)
        else
        text(max(x_lim),0.8*max(y_lim),['N_F_R_I = ',...
            num2str(length(xPlot)-1)],'HorizontalAlignment','right',...
            'FontSize',FS) 
        end
    end
    xlabel ('FRI (yr)','FontSize',FS)
    if i == length(zoneDiv)-1
    ylabel (['proportion OR density (x',num2str(binWidth),')'],...
        'FontSize',FS)
    text(-300,mean(y_lim),[{[char(site)]}],'fontsize',FS,'FontWeight',...
        'Bold','Rotation',90,'HorizontalAlignment','Center')
    end
    xlim (x_lim); ylim (y_lim)
    set (gca, 'TickDir','out','XTick',[0:200:x_lim(2)],...
        'XMinorTick','on','FontSize',FS);
    box off
    title([char(zoneText(i))],'FontSize',FS,'FontWeight',FW)
    end
end

%% FIGURE 7
figPosition = figPosition - [0.0165 0.0225 0 0];
figure (7); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Continuous fire history: peak magnitude, FRIs through time, and smoothed fire frequency');  
    % Continuous fire history characteristics.
subplot(3,5,1:4); hold off  % Peak Magnitude.
    x = Charcoal.ybpI;
    c_peaks = Charcoal.peak-CharThresh.pos(:,end); % peak CHAR - thresholds 
        % = CHAR exceeding threshold value.
    c_peaks (c_peaks < 0) = 0;    % Make all positives.
    c_peaks(end) = 0; % Make first sample not a peak.
    peak_in = zeros(length(c_peaks),2);
    peak_mag = zeros(length(c_peaks),1);
    for in = 1:length(c_peaks)
        if in == 1 & c_peaks(in) > 0
            peak_in(in,1) = in; % p-stop, if first sample is a peak.
            step = 1;
            while c_peaks(in+step) > 0
            peak_in(in,2) = in+step; % p-start
            step = step + 1;
            end
        else
            if c_peaks(in) > 0 & c_peaks(in-1) == 0
               peak_in(in,1) = in;  % p-stop
               step = 1;
            while c_peaks(in+step) > 0
            peak_in(in,2) = in+step; % p-start
            step = step + 1;
            end
            end
        end
        if peak_in(in,1) > 0 & peak_in(in,2) == 0
            peak_in(in,2) = peak_in(in,1);
        end
        if peak_in(in,2) > 0 
          peak_mag(peak_in(in,2),1) = sum(c_peaks(peak_in(in,1):...
              peak_in(in,2))) * [diff(peak_in(in,:))+1 * r];
          % [# cm^-2 yr^-1] * [yr peak^-1] = [# cm^-2 peak^-1] = Peak 
          % magnitude.
        end
    end
    H = bar(x,peak_mag);
    set(H,'FaceColor','k','BarWidth',1.0)
    hold on
    plot (x(peakScreenIn),logical(peakScreenIn)*0.8*max(peak_mag),'.',...
        'color',[.75 .75 .75]);
    plot (x(peakIn),[0.8*max(peak_mag)],'+r')
    y_lim = [0 1.1*max(peak_mag)];
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out','XTick',[0:1000:max(zoneDiv)],...
        'XTickLabel',[],'XLim',[zoneDiv(1), zoneDiv(length(zoneDiv))],...
        'Ylim',y_lim,'box','off','position',[0.1 .75 wm*range(zoneDiv) .12])
    if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(y_lim)*0.90 max(y_lim)],...
                '-k','color',[.5 .5 .5],'linewidth',2); 
        end
        text(mean(zoneDiv(z-1:z)),[max(y_lim)*0.95],char(zoneText(z-1)),...
            'horizontalalignment','center','fontweight','normal',...
            'fontsize',FS)
    end
    end
    if nGaps > 0    % If there are gaps in the record, then
        % fill in the graph with white
        buffer = 0; % [yr] Years to pad the years 
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
    end
    ylabel ([{'peak magnitude'}, {'(pieces cm^-^2 peak^-^1)'}])
    title ([{'Peak magnitude, FRIs, and fire frequ.'},...
        {''},{char(site)},{''}],...
        'FontSize',FS,'FontWeight','Bold','VerticalAlignment','Bottom');

subplot(3,5,9); % Smoothed fire frequency. 12/2008: added confidence-
    % interval estimates on smoothed fire frequency. 
    % Currently (July 2010), turned off. Confidence intervals on FRIs are
    % more robust. 
    hold off
    nBoot = 100;
    alpha = 0.05;
    nFiresBoot = NaN*ones(nBoot,1);
    Charcoal.peaksFrequCI = zeros(length(Charcoal.ybpI),2);
    Charcoal.peaksFrequ = zeros(size(Charcoal.ybpI));
    ff_smCi = zeros(length(Charcoal.ybpI),2);
    ff_sum_yr = PeakAnalysis.peakFrequ; % [yr] Window to sum peaks over.
    ff_sm_yr = PeakAnalysis.peakFrequ;  % [yr] Window to smooth peak 
        % frequencies over. 
    for i = 1:length(Charcoal.ybpI)
        if i < (ff_sum_yr/r)/2  % If start of record.     
            Charcoal.peaksFrequ (i) = sum(CharcoalCharPeaks...
                (1:floor(ff_sum_yr/r)/2+i,end)) *...
                [round(ff_sum_yr/r) / floor((ff_sum_yr/r)/2+i)];
            for j = 1:nBoot
                randomYears = randsample(CharcoalCharPeaks...
                (1:floor(ff_sum_yr/r)/2+i,end),length(CharcoalCharPeaks...
                (1:floor(ff_sum_yr/r)/2+i,end)),'true');
                nFiresBoot(j) = sum(randomYears);
            end
                Charcoal.peaksFrequCi(i,:) = prctile(nFiresBoot,...
                    [alpha/2*100 100-alpha/2*100])*...
                [round(ff_sum_yr/r) / floor((ff_sum_yr/r)/2+i)];
        else         
        if i > length(Charcoal.ybpI)-(ff_sum_yr/r)/2 % If end of record.
            Charcoal.peaksFrequ (i) = sum(CharcoalCharPeaks(i-...
                (ff_sum_yr/r)/2:end,end)) * [(ff_sum_yr/r) /...
              (length(CharcoalCharPeaks(i-(ff_sum_yr/r)/2:end,end)))];
            for j = 1:nBoot
                randomYears = randsample(CharcoalCharPeaks(i-...
                (ff_sum_yr/r)/2:end,end),length(CharcoalCharPeaks(i-...
                (ff_sum_yr/r)/2:end,end)),'true');
                nFiresBoot(j) = sum(randomYears);
            end
                Charcoal.peaksFrequCi(i,:) = prctile(nFiresBoot,...
                    [alpha/2*100 100-alpha/2*100])*...
                [(ff_sum_yr/r) /...
              (length(CharcoalCharPeaks(i-(ff_sum_yr/r)/2:end,end)))];
        else    % Else, it's the middle of the record.
            Charcoal.peaksFrequ (i) = sum(...
                CharcoalCharPeaks(ceil(i-0.5*(ff_sum_yr/r))+1:...
                ceil(i+0.5*round(ff_sum_yr/r)),end));
%             for j = 1:nBoot
%                 randomYears =...
%                 randsample(CharcoalCharPeaks(ceil(i-0.5*(ff_sum_yr/r))+1:...
%                 ceil(i+0.5*round(ff_sum_yr/r)),end),...
%                 length(CharcoalCharPeaks(ceil(i-0.5*(ff_sum_yr/r))+1:...
%                 ceil(i+0.5*round(ff_sum_yr/r)),end)),'true');
%                 nFiresBoot(j) = sum(randomYears);
%             end
%                 Charcoal.peaksFrequCi(i,:) = prctile(nFiresBoot,...
%                     [alpha/2*100 100-alpha/2*100]);            
        end
        end
    end
    ff_sm = smooth(Charcoal.peaksFrequ,ff_sm_yr/r,'lowess');
%     ff_smCi(:,1) = smooth(Charcoal.peaksFrequCi(:,1),ff_sm_yr/r,'lowess');
%     ff_smCi(:,2) = smooth(Charcoal.peaksFrequCi(:,2),ff_sm_yr/r,'lowess');
    hold on
%     X = [Charcoal.ybpI; flipud(Charcoal.ybpI)];
%     Y = [ff_smCi(:,1); flipud(ff_smCi(:,2))];
%     h = fill(X,Y,[0.8 0.8 0.8]);
%     set(h,'edgecolor',[0.8 0.8 0.8])
    plot(Charcoal.ybpI,ff_sm,'k','linewidth',1)
%     plot(Charcoal.ybpI,ff_smCi,'--k');
    y_lim = [0 1.1*max(ff_sm)];
    if y_lim(1) < 0
        y_lim = [0 y_lim(2)];
    end
    if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(y_lim)*0.90 max(y_lim)],'-k',...
                'color',[.5 .5 .5],'linewidth',2); 
        end
    end
    end
 if nGaps > 0    % If there are gaps in the record, then
        % fill in the graph with white
        buffer = 0;%ff_sum_yr/2; % [yr] Years to pad the years 
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
 end
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out','XTick',[0:1000:max(zoneDiv)],...
        'XTickLabel',[0:max(zoneDiv)/1000],'XLim',[zoneDiv(1), zoneDiv(length(zoneDiv))],...
        'box','off','ylim',y_lim,'position',[0.1 .45 wm*range(zoneDiv) .12])
    ylabel ([{'fire frequency'},...
            {['(fires ' num2str(PeakAnalysis.peakFrequ) ' yr^-^1)']}])
    xlabel ('time (cal. yr BP x 1000)','FontSize',FS)  

subplot(3,5,11:14); hold off    % FRIs through time.
    peak_yrs = Charcoal.ybpI(CharcoalCharPeaks(:,end) > 0); % Years
        % with charocal peaks.
    if length(peak_yrs) > 2
    [FRIyr FRI smFRIyr smFRI smFRIci] = smoothFRI (Charcoal.ybpI,...
        CharcoalCharPeaks(:,end),...
        PeakAnalysis.peakFrequ,alpha,nBoot,1,0);  
    x_lim = [0 1.1*max(FRI)];
    X = [smFRIyr'; flipud(smFRIyr')];
    Y = [smFRIci(:,1); flipud(smFRIci(:,2))];
    plot(FRIyr,FRI,'sk','MarkerFaceColor',[.75 .75 .75],'MarkerSize',4)
    hold on
    plot(smFRIyr,smFRI,'-k')
    H = fill(X,Y,[0.8 0.8 0.8]);
%     legend('FRI','mean FRI',[num2str(100*(1-alpha)) '% CI'],1)
    set(H,'edgecolor',[0.8 0.8 0.8])
    plot(FRIyr,FRI,'sk','MarkerFaceColor',[.75 .75 .75],'MarkerSize',4)
    plot(smFRIyr,smFRI,'-k')
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out','XTick',[0:1000:max(zoneDiv)],...
        'Ylim',x_lim,'XTickLabel',[],'XLim',[zoneDiv(1),...
        zoneDiv(length(zoneDiv))],'box','off','position',...
        [0.1 .6 wm*range(zoneDiv) .12],'yScale','linear')
       % Resampled smoothed FRIs to same intervals as interpolated samples.
       if length(smFRI) > 2
          yis = interp1(smFRIyr,smFRI,Charcoal.ybpI(Charcoal.ybpI<max(x)));   
       else
           yis = -999+zeros(size(smFRI));
           warning('Less than 3 FRIs - smoothing cannot be done; values set to -999')
       end
%     plot(smFRIyr,smFRI,'k','linewidth',1)    % Plot smoothed FRIs.
        if nGaps > 0   % If there are gaps in the record, then fill in the 
                    % graph with white
        buffer = 0; % [yr] Years to pad the years with missing values. 
        for i = 1:nGaps
            y = get(gca,'ylim');
            x = [Charcoal.ybp(gapIn(i,1))-buffer ...
                Charcoal.ybp(gapIn(i,1))-buffer ...
                Charcoal.ybp(gapIn(i,2))+buffer ...
                Charcoal.ybp(gapIn(i,2))+buffer];
            patch(x,[y(1) y(2) y(2) y(1)],'w','edgecolor',[.75 .75 .75],...
                'facecolor',[.75 .75 .75])
        end
        end
        if length(zoneDiv) > 2
        for z = 2:length(zoneDiv)
            if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(x_lim)*0.90 max(x_lim)],...
                '-k','color',[.5 .5 .5],'linewidth',2); 
            end
        end
        end
        ylabel ({['FRI (yr fire^-^1)'],...
        [num2str(PeakAnalysis.peakFrequ) '-yr mean'],...
        [num2str((1-alpha)*100) '% CI']},'FontSize',FS)
        grid on
    else
        yis = NaN;
        smFRI = NaN;
    end
    
%% FIGURE 8
figPosition = figPosition - [0.0165 0.0225 0 0];
figure (8); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Between-zone comparisons of raw charcoal distributions')
    col = [0 0 0; .5 .5 .5; 0 1 0; 0 0 1; 1 0 1; 1 0 0; 1 0.4 0; 0 0.5 0];
    Charcoal.accZone = NaN(1000,length(zoneDiv)-1);    % Space for 
    for i = 1:length(zoneDiv)-1
    CHAR_z = Charcoal.acc(Charcoal.ybp >= zoneDiv(i) &...
        Charcoal.ybp< zoneDiv(i+1)); % Select only years between zoneDiv(i)
        % and zoneDiv(i+1) -- veg. zones.

    Charcoal.accZone(1:length(CHAR_z),i) = CHAR_z; % Store CHAR values for 
                % this zone.
    subplot (3,5,[6 7 8 11 12 13])
    [f,x] = ecdf(CHAR_z); 
    plot(x,f,'color', col(i,:),'linewidth',2);
    hold on
    end
    set(gca,'tickdir','out')
    xlabel ('CHAR (pieces cm^-^2 yr^-^1)')
    ylabel('cum. prop.');
    title ('CDFs of zone-specific raw CHAR, with KS test results')
    box off
    legend (char(zoneText))

    % KS test of CHAR distributions between veg. zones
    if length(zoneDiv) > 2   % If there is more than one zone, compare the 
        % CHAR distributions between zones.
    for i = 1:length(zoneDiv)-2
        for j = 2:length(zoneDiv)-1
        [h(i,j-1) pKS(i,j-1) k(i,j-1)] = kstest2(Charcoal.accZone(:,i),...
            Charcoal.accZone(:,j));
        end
    end
    pKSResults(2:length(pKS)+1,2:length(pKS)+1) = pKS;
    pKSResults(1,2:end) = 2:length(pKS)+1;
    pKSResults(2:end,1) = 1:length(pKS);
    text(0.5*max(Charcoal.acc),0.63,'KS p-value matrix:','FontWeight',...
        'bold');
    text(0.45*max(Charcoal.acc),0.48,...
        [num2str(round(pKSResults*1000)/1000)]);
    text(0.45*max(Charcoal.acc),0.58,'Zone','FontWeight','Bold',...
        'BackgroundColor','w','HorizontalAlignment','Center')
    else
        h = -999; pKS = -999; k = -999;   % dummy variables
    end
    text (0,1.5,[char(site) ': Between-zone comparisons of raw CHAR distributions'],...
        'FontSize',FS+2,'FontWeight','Bold')

    % Box plots of CHAR for each zone
    subplot (3,5,[9 10 14 15])
    boxplot(fliplr(Charcoal.accZone),'colors','kkk',...
        'symbol','.k')
    set(gca,'xticklabel',fliplr(1:length(zoneDiv)-1),'tickdir','out','ylim',...
        [min(min(Charcoal.accZone)) max(max(Charcoal.accZone))])
    ylabel ('CHAR (pieces cm^-^2 yr^-^1)')
    xlabel ('zone number')
    title ('Box plots of raw CHAR per zone')
    box off
    
%% FIGURE 9
if Results.allFigures == 1
figPosition = figPosition - [0.0165 0.0225 0 0];
figure (9); clf; set(gcf,'color','w','units','normalized',...
    'position',[figPosition],'name',...
    'Alternative displays of threshold value(s)')
    CHAR_ratio = Charcoal.accI./Charcoal.accIS;
    CHAR_residuals = Charcoal.accI - Charcoal.accIS;
    y = CharcoalCharPeaks(:,end);  % Charocal peaks
        y(y==0) = -99;
    if PeakAnalysis.cPeak ~= -99
        if PeakAnalysis.cPeak == 1
            y1 = Charcoal.accIS + CharThresh.pos(:,end);
            y2 = y1./Charcoal.accIS;
            y3 = CharThresh.pos(:,end);
        end
        if PeakAnalysis.cPeak == 2
            y1 = Charcoal.accIS .* CharThresh.pos(:,end);
            y2 = CharThresh.pos(:,end);
            y3 = Charcoal.accIS-...
                (Charcoal.accIS ./ CharThresh.pos(:,end));
        end
    else
        y1 = CharThresh.pos(:,end);
        y2 = CharThresh.pos(:,end)./Charcoal.accIS;
        y3 = CharThresh.pos(:,end)-Charcoal.accIS;
    end

    subplot(3,1,1)
    stairs(Charcoal.ybpI,Charcoal.accI,'color',[.75 .75 .75])
    hold on
    plot(Charcoal.ybpI,Charcoal.accIS,'linewidth',1.5,'color','k')
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out',...
        'XTick',[0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000])
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
    ylim([0 max(Charcoal.accI)])
    box off
    plot(Charcoal.ybpI,y1,'r','linewidth',1.5)
    hold on
    plot (Charcoal.ybpI,...%(CharThresh.pos(:,end)>0),...
        max(Charcoal.accI)*0.80*y,'+k','color','r')
    plot (Charcoal.ybpI(peakScreenIn),max(Charcoal.accI)*0.80*logical(peakScreenIn),'.',...
        'color',[.75 .75 .75]);
    ylabel('CHAR (pieces cm^-^2 yr^-^1)')
    title ('(a) Relationship between C_b_a_c_k_g_r_o_u_n_d and peak threshold')
    set(gca,'position',[0.1 .75 wm*range(zoneDiv) .175])
    legend ('C_i_n_t_e_r_p_o_l_a_t_e_d', 'C_b_a_c_k_g_r_o_u_n_d',...
        'C_t_h_r_e_s_h_o_l_d','peaks',...
        'Location','Eastoutside');
    grid on

    subplot(3,1,2)
    stairs(Charcoal.ybpI,CHAR_ratio,'color',[.75 .75 .75])
    hold on
    plot(Charcoal.ybpI,y2,'r','linewidth',1.5)
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out',...
        'XTick',[0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000])
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
%     if max(y2) == Inf
%         ylim([1 NaNmedian(y2)*5])
%     else
%         ylim([1 2*max(y2)])
%     end
    box off; grid on
    ylabel ([{'threshold'} {'(threshold / low-frequency CHAR)'}])
    title ('(b) Peak threshold as ratios of C_b_a_c_k_g_r_o_u_n_d')
    set(gca,'position',[0.1 .50 wm*range(zoneDiv) .175])
    legend ('C_p_e_a_k','C_t_h_r_e_s_h_o_l_d','peaks identified',...
        'Location','Eastoutside')

    subplot(3,1,3)
    stairs(Charcoal.ybpI,CHAR_residuals,'color',[.75 .75 .75])
    hold on
    plot(Charcoal.ybpI,y3,'r','linewidth',1.5)
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on',...
        'YMinorTick','on','TickDir','out',...
        'XTick',[0:1000:max(zoneDiv)],'XTickLabel',[0:max(zoneDiv)/1000])
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
    if min(y3) <= 0
        ylim([0 abs(max(y3)*2)])
    else
        ylim([0 max(y3)*2])
    end
    box off; grid on
    xlabel ('time (cal. yr BP x 1000)')
    ylabel ([{'threshold'} {'(threshold - low-frequency CHAR)'}])
    title ('(c) Peak threshold as residuals of C_b_a_c_k_g_r_o_u_n_d')
    set(gca,'position',[0.1 .20 wm*range(zoneDiv) .175])
     legend ('C_p_e_a_k','C_t_h_r_e_s_h_o_l_d','peaks identified',...
         'Location','Eastoutside')
end

%% ADD DATA TO INPUT DATA STRUCTURES
Charcoal.peakInsig = zeros(size(Charcoal.ybpI));
Charcoal.peakInsig(peakScreenIn) = 1;
Charcoal.peakMagnitude = peak_mag;
Charcoal.smoothedFRI = smFRI;
Charcoal.smoothedFireFrequ = ff_sm; 

%% CREATE OUTPUT VARIABLES AND SAVE DATA
if Results.saveFigures == 1
set(gcf,'PaperPositionMode','auto')
if Results.allFigures == 1
figure (1); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')
    orient(gcf,'landscape')
    print -dpdf -r300 01_pretreatment.pdf
    print -dtiff -r300 01_pretreatment.tif
figure (2); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 02_threshold_determination.pdf
    print -dtiff -r300 02_threshold_determination.tif
end
figure (3); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 03_CHAR_analysis.pdf
    print -dtiff -r300 03_CHAR_analysis.tif    
figure (4); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')
    orient(gcf,'landscape')
    print -dpdf -r300 04_CHAR_peak_sens.pdf
    print -dtiff -r300 04_CHAR_peak_sens.tif
figure (5); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 05_cum_peaks_through_time.pdf
    print -dtiff -r300 05_cum_peaks_through_time.tif
figure (6); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 06_FRI_dists.pdf
    print -dtiff -r300 06_FRI_dists.tif    
figure (7); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 07_continuous_fire_hx.pdf
    print -dtiff -r300 07_continuous_fire_hx.tif 
figure (8); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 08_CHAR_dists.pdf
    print -dtiff -r300 08_CHAR_dists.tif     
figure (9); set(gcf,'PaperPositionMode','auto','PaperType','uslegal')  
    orient(gcf,'landscape')
    print -dpdf -r300 09_threshold_details.pdf
    print -dtiff -r300 09_threshold_details.tif 
end
if Results.save == 1
charResults = NaN*ones(length(Charcoal.cmI),33);
charResults(:,1:22) = [Charcoal.cmI Charcoal.ybpI Charcoal.countI...
Charcoal.volI Charcoal.conI Charcoal.accI Charcoal.accIS...
Charcoal.peak CharThresh.pos CharThresh.neg(:,end) CharThresh.SNI...
CharThresh.GOF(:,1) CharcoalCharPeaks Charcoal.peakInsig...
Charcoal.peakMagnitude Charcoal.smoothedFireFrequ];
charResults(1:length(yis),23) = yis;
charResults(1:length(zoneDiv)-1,24:33) = FRI_params_zone;
    if min(fileName(end-2:end) == ['xls']) > 0  % If using an .xls file...
        xlswrite(fileName,charResults,'charResults','a2:x2500')
    else    % If using a .csv file
        outputResults(charResults,fileName,site);
    end
end