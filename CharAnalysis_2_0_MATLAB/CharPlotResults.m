function CharPlotResults(Charcoal, Pretreatment, PeakAnalysis, ...
    CharThresh, Smoothing, site, Results, gapIn, Post)
% CharPlotResults   Generate all output figures for CharAnalysis.
%   CharPlotResults(Charcoal, Pretreatment, PeakAnalysis, CharThresh, ...
%       Smoothing, site, Results, gapIn, Post)
%
%   Produces Figures 3–9 (and saves them if Results.saveFigures == 1).
%   All analytical values come from the Post struct produced by
%   CharPostProcess. No computation is performed here.
%
%   plotyy() has been replaced throughout with yyaxis left / yyaxis right.
%   smooth()  has been replaced with charLowess().

%% ── Unpack frequently used values ───────────────────────────────────────
zoneDiv  = Pretreatment.zoneDiv;
r        = Pretreatment.yrInterp;
yr       = Smoothing.yr;
transform = Pretreatment.transform;
nGaps    = size(gapIn, 1);

peakIn            = Post.peakIn;
peakScreenIn      = Post.peakScreenIn;
CharcoalCharPeaks = Post.CharcoalCharPeaks;
threshIn          = Post.threshIn;
peak_mag          = Post.peak_mag;
ff_sm             = Post.ff_sm;
FRIyr             = Post.FRIyr;
FRI               = Post.FRI;
smFRIyr           = Post.smFRIyr;
smFRI             = Post.smFRI;
smFRIci           = Post.smFRIci;
alpha             = Post.alpha;
FRI_params_zone   = Post.FRI_params_zone;

%% ── Layout helpers ───────────────────────────────────────────────────────
if max(Charcoal.ybp) > 15000
    wm = 0.5 * 2.8e-005;
elseif max(Charcoal.ybp) > 5500
    wm = 1.25 * 2.8e-005;
elseif max(Charcoal.ybp) > 2500
    wm = 1.5 * 2.8e-005;
else
    wm = 5 * 2.8e-005;
end

FS  = 8;
FW  = 'bold';
zoneText    = {'Zone 1','Zone 2','Zone 3','Zone 4','Zone 5','Zone 6','Zone 7'};
figPosition = [0.1013  0.1933  0.8648  0.6943];

%% ── Shared helper: draw gap patches ─────────────────────────────────────
    function drawGaps(ax, gapIn, Charcoal)
        if isempty(ax), ax = gca; end
        yl = get(ax, 'ylim');
        for gi = 1:size(gapIn,1)
            xg = [Charcoal.ybp(gapIn(gi,1))  Charcoal.ybp(gapIn(gi,1)) ...
                  Charcoal.ybp(gapIn(gi,2))  Charcoal.ybp(gapIn(gi,2))];
            patch(xg, [yl(1) yl(2) yl(2) yl(1)], 'w', ...
                'edgecolor',[.75 .75 .75], 'facecolor',[.75 .75 .75])
        end
    end

%% ── Shared helper: CHAR y-label string ──────────────────────────────────
    function s = charLabel(transform)
        switch transform
            case 1,  s = 'log CHAR (pieces cm^-^2 yr^-^1)';
            case 2,  s = 'ln CHAR (pieces cm^-^2 yr^-^1)';
            otherwise, s = 'CHAR (pieces cm^-^2 yr^-^1)';
        end
    end

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 3 — C_interpolated, C_background, C_peak
%%═══════════════════════════════════════════════════════════════════════════
figure(3); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Resampled charcoal (C_interpolated), low-frequency trends (C_background), and detrended series (C_peak)')

% ── Panel (a): C_int and C_background ────────────────────────────────────
subplot(3,5,1:4)
x  = Charcoal.ybpI;
y  = Charcoal.accI;
y2 = Charcoal.accIS;
H  = bar(x, y); hold on
set(H, 'facecolor',[0 0 0], 'edgecolor',[0 0 0])
plot(x, y2, 'color',[.5 .5 .5], 'linewidth', 2)
xlim([zoneDiv(1), zoneDiv(end)])
ylim([0, 1.1*max(y)])
if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv], [max(y)*1.01 max(y)*1.1], '-k', ...
        'color',[.5 .5 .5], 'linewidth', 2)
end
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on', ...
    'YMinorTick','on','TickDir','out', ...
    'XTick',0:1000:max(zoneDiv),'XTickLabel',[], ...
    'Position',[0.1 .55 wm*range(zoneDiv) .175])
box off
ylabel(charLabel(transform), 'FontSize', FS)
for z = 2:length(zoneDiv)
    if z < length(zoneDiv)
        plot([zoneDiv(z) zoneDiv(z)],[max(y)*1.01 max(y)*1.1], ...
            '-k','color',[.5 .5 .5],'linewidth',2)
    end
    if length(zoneDiv) > 2
        text(mean(zoneDiv(z-1:z)), max(y)*1.05, zoneText(z-1), ...
            'horizontalalignment','center','fontweight','normal','fontsize',FS)
    end
end
title({['(a) ',char(site),': C_i_n_t_e_r_p_o_l_a_t_e_d (',num2str(r),' yr)', ...
    ' and C_b_a_c_k_g_r_o_u_n_d defined by ', ...
    num2str(yr(end)),'-yr trends']}, 'FontSize',FS,'FontWeight',FW)

% ── Panel (b): C_peak and thresholds ─────────────────────────────────────
subplot(3,5,10)
y  = Charcoal.peak;
y2 = CharThresh.pos(:,end);
y3 = CharThresh.neg(:,end);
if PeakAnalysis.cPeak == 1
    H1 = bar(x, y, 1);
else
    H1 = bar(x, y, 1, 'BaseValue', 1);
end
set(H1,'facecolor',[0 0 0],'edgecolor',[0 0 0])
hold on
plot(x, y2, 'r'); plot(x, y3, 'r')
if PeakAnalysis.cPeak == 1
    plot([min(x) max(x)], [0 0], 'k')
elseif PeakAnalysis.cPeak == 2
    plot([min(x) max(x)], [1 1], 'k')
end
plot(x(peakScreenIn), ones(size(peakScreenIn))*0.8*max(y), '.', ...
    'color',[.75 .75 .75])
plot(x(peakIn), 0.8*max(y), '+k')
xlim([zoneDiv(1), zoneDiv(end)])
ylim([min(y3), 1.1*max(y)])
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
    'Position',[0.1 .32 wm*range(zoneDiv) .175],'TickDir','out', ...
    'XTick',0:1000:max(zoneDiv),'XTickLabel',0:max(zoneDiv)/1000)
if PeakAnalysis.cPeak == 1
    title({['(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d - C_b_a_c_k_g_r_o_u_n_d), ' ...
        'thresholds defining C_n_o_i_s_e, and peaks'], ...
        'identified (gray dots fail to pass peak-magnitude test)'}, ...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
elseif PeakAnalysis.cPeak == 2
    title(['(b) C_p_e_a_k (C_i_n_t_e_r_p_o_l_a_t_e_d / C_b_a_c_k_g_r_o_u_n_d), ' ...
        'thresholds defining C_n_o_i_s_e, and peaks identified ' ...
        '(gray peaks fail to pass peak-magnitude test)'], ...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
else
    title(['(b) C_i_n_t_e_r_p_o_l_a_t_e_d, thresholds defining C_n_o_i_s_e, ' ...
        'and peaks identified (gray peaks fail to pass peak-magnitude test)'], ...
        'VerticalAlignment','middle','FontSize',FS,'FontWeight','bold')
end
xlabel('time (cal. yr BP x 1000)','FontSize',FS)
ylabel(charLabel(transform),'FontSize',FS)
box off

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 4 — Sensitivity to alternative thresholds and quality of record
%%═══════════════════════════════════════════════════════════════════════════
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(4); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Sensitivity to alternative thresholds and quality of record')

% ── Panel (a): C_int with C_back and threshold ────────────────────────────
subplot(3,5,1:4)
x  = Charcoal.ybpI;
y  = Charcoal.accI;
y2 = Charcoal.accIS;
y3 = CharThresh.pos(:,end);
y4 = CharThresh.neg(:,end);
y5 = CharcoalCharPeaks;  y5(y5==0) = -99;
H  = bar(x, y, 1); hold on
set(H,'facecolor',[0 0 0],'edgecolor',[0 0 0])
plot(x, y2,'color',[.5 .5 .5],'linewidth',2)
if PeakAnalysis.cPeak == 1
    plot(x, y2+y3,'r');  plot(x, y2+y4,'r')
else
    plot(x, y2.*y3,'r'); plot(x, y2.*y4,'r')
end
xlim([zoneDiv(1), zoneDiv(end)]); ylim([0, 1.15*max(y)])
if length(zoneDiv) > 2
    plot([zoneDiv zoneDiv],[max(y)*1.01 max(y)*1.1],'-k','color',[.5 .5 .5])
end
for z = 2:length(zoneDiv)
    if z < length(zoneDiv)
        plot([zoneDiv(z) zoneDiv(z)],[max(y)*1.01 max(y)*1.1], ...
            '-k','color',[.5 .5 .5],'linewidth',2)
    end
    if length(zoneDiv) > 2
        text(mean(zoneDiv(z-1:z)),max(y)*1.05,char(zoneText(z-1)), ...
            'horizontalalignment','center','fontweight','normal','fontsize',FS)
    end
end
plot(x, max(y)*0.78*y5(:,1), '.k','color',[.5 .5 .5])
plot(x, max(y)*0.85*y5(:,2), '.k','color',[.5 .5 .5])
plot(x, max(y)*0.92*y5(:,3), '.k','color',[.5 .5 .5])
if     sum(y5(:,1) - y5(:,end)) == 0;  mIndex = 0.78;
elseif sum(y5(:,3) - y5(:,end)) == 0;  mIndex = 0.92;
else;                                   mIndex = 0.85;
end
plot(x, max(y)*mIndex*y5(:,end), '.k','color',[1 1 1])
plot(x, max(y)*mIndex*y5(:,end), '+k')
ylabel(charLabel(transform),'FontSize',FS)
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
    'TickDir','out','XTick',0:1000:max(zoneDiv),'XTickLabel',[], ...
    'Position',[0.1 .55 wm*range(zoneDiv) .175])
title([{char(site)},{''},{'(a) C_i_n_t_e_r_p_o_l_a_t_e_d, C_b_a_c_k_g_r_o_u_n_d, and peak ID (+)'}], ...
    'FontSize',FS,'FontWeight',FW)
box off

% ── Panel (b): Threshold sensitivity ─────────────────────────────────────
subplot(3,5,5)
if PeakAnalysis.threshType == 1        % Global threshold
    cpk   = Charcoal.peak;
    y3tot = sum(CharcoalCharPeaks);    % total peaks per threshold column

    if PeakAnalysis.cPeak == 1
        xPoss  = CharThresh.possible;
        xPos2  = xPoss(xPoss>0);
    else
        xPoss  = CharThresh.possible;
        xPos2  = xPoss(xPoss>1);
    end

    [n, xh] = histcounts(cpk, xPoss, 'Normalization','probability');
    xh_c = xh(1:end-1) + diff(xh)/2;   % bin centres
    H1 = bar(xh_c, n, 'BarWidth', 1.0);
    set(H1,'FaceColor',[.75 .75 .75],'EdgeColor',[.75 .75 .75])
    hold on

    % ── yyaxis replaces plotyy ────────────────────────────────────────────
    yyaxis left
    if PeakAnalysis.threshMethod > 1
        y2pdf = CharThresh.noisePDF;
        if PeakAnalysis.cPeak == 1
            plot(xh_c, y2pdf * mean(diff(xh_c)), 'k--', 'linewidth', 1.5)
        else
            xIn  = length(xh_c) - length(xh_c(xh_c>=1));
            xIn2 = xIn - (length(xh_c) - length(y3tot));
            plot(xh_c, y2pdf * mean(diff(xh_c)), 'k--', 'linewidth', 1.5)
        end
    else
        plot(-99, -99, 'k--')           % invisible placeholder
    end
    ylabel('relative frequency','fontsize',FS)
    ylim([0, 1.01*max(n)])
    set(gca,'Ycolor','k')

    yyaxis right
    if PeakAnalysis.cPeak == 1
        xPlotThresh = xh_c(xh_c>0);
        yPlotThresh = y3tot(1:length(xPlotThresh));
    else
        xIn  = length(xh_c) - length(xh_c(xh_c>=1));
        xIn2 = xIn - (length(xh_c) - length(y3tot));
        xPlotThresh = xh_c(xIn:end);
        yPlotThresh = y3tot(xIn2:end);
    end
    plot(xPlotThresh, yPlotThresh, 'k-', 'linewidth', 1.5)
    ylabel('# of peaks identified','fontsize',FS,'rotation',270, ...
        'verticalalignment','bottom','color','k')
    ylim([0, 0.99*max(y3tot)])
    set(gca,'Ycolor','k')

    % Restore left axis for threshold markers
    yyaxis left
    for t = 1:4
        plot([CharThresh.pos(1,t) CharThresh.pos(1,t)],[0 max(n)],'--', ...
            'color',[0.75 0.75 0.75])
        if t == 4
            plot([CharThresh.pos(1,t) CharThresh.pos(1,t)],[0 max(n)],'-k')
            text(CharThresh.pos(1,t), 0, '<', ...
                'FontSize',FS+2,'Rotation',90,'FontWeight','bold')
        end
    end
    if PeakAnalysis.threshMethod > 1
        text(CharThresh.pos(1,end)*1.25, 0.9*max(n), ...
            {[num2str(PeakAnalysis.threshValues(end)*100),'^t^h percentile'], ...
             ['= ',num2str(round(max(xh_c(xh_c<=CharThresh.pos(1,end)))*1e4)/1e4)]}, ...
            'FontSize',FS,'HorizontalAlignment','left','backgroundcolor','w')
    else
        text(CharThresh.pos(1,end)*1.25, 0.9*max(n), ...
            {'threshold value', ...
             ['= ',num2str(round(max(xh_c(xh_c<=CharThresh.pos(1,end)))*1e4)/1e4)]}, ...
            'FontSize',FS,'HorizontalAlignment','left','backgroundcolor','w')
    end
    set(gca,'TickDir','out','FontSize',FS,'XMinorTick','on', ...
        'XLim',[min(Charcoal.peak) 0.75*max(Charcoal.peak)], ...
        'position',[0.17+wm*range(zoneDiv) .55 .16 .175])
    box off
    if PeakAnalysis.cPeak == 1
        xlabel('residual CHAR value (pieces cm^-^2 yr^-^1)','FontSize',FS)
    else
        xlabel('CHAR ratio (pieces cm^-^2 yr^-^1)','FontSize',FS)
    end
    if PeakAnalysis.threshMethod > 1
        title({'(b) High-frequency CHAR distribution,','est. noise dist., and threshold values'},'FontWeight',FW)
    else
        title({'(b) High-frequency CHAR distribution','and threshold values'},'FontWeight',FW)
    end

else    % ── Local threshold: mFRI sensitivity by zone ────────────────────
    plot_in = fliplr(1:length(zoneDiv)-1);
    for i = 1:length(plot_in)
        xPlot1 = Charcoal.ybpI(CharcoalCharPeaks(:,1)>0);
        xPlot2 = Charcoal.ybpI(CharcoalCharPeaks(:,2)>0);
        xPlot3 = Charcoal.ybpI(CharcoalCharPeaks(:,3)>0);
        xPlot1 = xPlot1(xPlot1>=zoneDiv(i) & xPlot1<zoneDiv(i+1));
        xPlot2 = xPlot2(xPlot2>=zoneDiv(i) & xPlot2<zoneDiv(i+1));
        xPlot3 = xPlot3(xPlot3>=zoneDiv(i) & xPlot3<zoneDiv(i+1));
        FRI_thresh = NaN*ones(200,3);
        FRI_thresh(1:length(diff(xPlot1)),1) = diff(xPlot1)';
        FRI_thresh(1:length(diff(xPlot2)),2) = diff(xPlot2)';
        FRI_thresh(1:length(diff(xPlot3)),3) = diff(xPlot3)';
        mFRI_thresh = NaN(3,1);  mFRI_ci = NaN(3,2);
        for j = 1:3
            vals = FRI_thresh(FRI_thresh(:,j)>0, j);
            mFRI_thresh(j) = mean(vals);
            if length(vals) > 1
                bm = NaN(1000,1);
                for b = 1:1000
                    bm(b) = mean(vals(randi(length(vals),length(vals),1)));
                end
                mFRI_ci(j,:) = prctile(bm,[2.5 97.5]);
            end
        end
        xp = [plot_in(i)-0.25, plot_in(i), plot_in(i)+0.25];
        errorbar(xp, mFRI_thresh, mFRI_thresh-mFRI_ci(:,1), ...
            mFRI_ci(:,2)-mFRI_thresh, '.k','color',[.5 .5 .5])
        hold on
        in2 = find(PeakAnalysis.threshValues == PeakAnalysis.threshValues(end),1);
        errorbar(xp(in2), mFRI_thresh(in2), ...
            mFRI_thresh(in2)-mFRI_ci(in2,1), mFRI_ci(in2,2)-mFRI_thresh(in2), '+k')
    end
    set(gca,'TickDir','out','FontSize',FS, ...
        'xtick',1:length(zoneDiv)-1,'xticklabel',fliplr(1:length(zoneDiv)-1), ...
        'position',[0.17+wm*range(zoneDiv) .55 .10 .175])
    box off
    xlabel('zone')
    ylabel({'zone-specific mean FRI','+/- 95% ci (years fire^-^1)'})
    title({'(b) Sensitivity to','alternative thresholds'},'fontweight','bold')
end

% ── Panel (c): Local SNI time series ─────────────────────────────────────
subplot(3,5,11:13)
if PeakAnalysis.threshType == 1
    SNI_plot = CharThresh.SNI .* ones(size(Charcoal.peak));
else
    SNI_plot = CharThresh.SNI;
end
plot(Charcoal.ybpI, SNI_plot, 'k')
y_lim = [0 10];
xlim([zoneDiv(1), zoneDiv(end)]); ylim(y_lim)
if y_lim(2) < 20;   y_tick = 0:2:y_lim(2);
elseif y_lim(2) < 50; y_tick = 0:5:y_lim(2);
else;               y_tick = 0:10:y_lim(2);
end
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
    'Position',[0.1 .35 wm*range(zoneDiv) .12],'TickDir','out', ...
    'XTick',0:1000:max(zoneDiv),'XTickLabel',0:max(zoneDiv)/1000,'YTick',y_tick)
hold on; plot([zoneDiv(1) zoneDiv(end)],[3 3],'k--'); box off
xlabel('time (cal. yr BP x 1000)','FontSize',FS)
ylabel('signal-to-noise index','FontSize',FS)
title('(c) Local signal-to-noise index','FontSize',FS,'FontWeight',FW)

% ── Panel (d): Global SNI boxplot ────────────────────────────────────────
subplot(3,5,14)
boxplot(CharThresh.SNI,'colors','kkk','symbol','.k')
set(gca,'TickDir','out','FontSize',FS,'YMinorTick','off', ...
    'position',[0.17+wm*range(zoneDiv) .35 .10 .12], ...
    'xtick',1,'xticklabel','','ylim',y_lim,'YTick',y_tick,'yscale','linear')
box off; grid off; hold on
plot([zoneDiv(1) zoneDiv(end)],[3 3],'k--')
title('(d) Global signal-to-noise index','FontSize',FS,'FontWeight',FW)
xlabel({'global signal-to-','noise distribution'},'FontSize',FS)
ylabel('')
text(1.2, median(CharThresh.SNI), num2str(nanmedian(CharThresh.SNI)), ...
    'FontSize',FS,'BackgroundColor','w')

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 5 — Cumulative peaks through time
%%═══════════════════════════════════════════════════════════════════════════
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(5); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Cumulative peaks through time')

xPlot = Charcoal.ybpI(CharcoalCharPeaks(:,end) > 0);
yPlot = flipud(cumsum(logical(find(CharcoalCharPeaks(:,end)>0))));
if ~isempty(xPlot)
    plot(xPlot, yPlot, '.k')
    set(gca,'XDir','reverse','TickDir','out','XMinorTick','on', ...
        'YMinorTick','on','FontSize',FS, ...
        'XTick',0:1000:max(zoneDiv),'XTickLabel',0:max(zoneDiv)/1000, ...
        'Position',[.1 .1 0.8*wm*range(zoneDiv) max(yPlot)*0.005])
    if nGaps > 0
        drawGaps(gca, gapIn, Charcoal)
        legend('Cumulative peaks','missing values')
    end
    xlabel('time (cal. yr BP x 1000)','FontSize',FS)
    ylabel('cumulative number of peaks','FontSize',FS)
    xlim([min(Charcoal.ybp), max(zoneDiv)])
    title([char(site),': Cumulative fires as a function of time'], ...
        'FontSize',FS,'FontWeight','bold')
    box off; grid on
end

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 6 — FRI distributions by zone with Weibull models
%%═══════════════════════════════════════════════════════════════════════════
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(6); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Fire return intervals by zone, with Weibull models if GOF test is passed')

binWidth = 20;
nZones   = length(zoneDiv) - 1;
x_lim    = [0 800];
y_lim    = [0 0.35];

for i = 1:nZones
    inPlot = fliplr(1:nZones);
    subplot(3,5,inPlot(i))

    % Retrieve per-zone data from Post (computed by CharPostProcess)
    if ~isfield(Post,'zone') || ~isfield(Post.zone(i),'FRI') || ...
            isempty(Post.zone(i).FRI)
        text(0.25, 0.5, 'insufficient data'); box off; continue
    end

    FRIz        = Post.zone(i).FRI;
    FRI_bin_c   = Post.zone(i).FRI_bin_centers;
    FRI_freq    = Post.zone(i).FRI_freq;
    param       = Post.zone(i).param;
    pKS         = Post.zone(i).pKS;
    wbl_est     = Post.zone(i).wbl_est;
    wbl_a_ci    = Post.zone(i).wbl_a_ci;
    wbl_b_ci    = Post.zone(i).wbl_b_ci;
    mFRIz       = Post.zone(i).mean_mFRI;
    mFRI_ci_z   = Post.zone(i).mean_mFRI_ci;
    xPlotZ      = Post.zone(i).xPlot;

    if max(FRIz) > 5000
        disp(['WARNING, Figure 6: FRIs in Zone ',num2str(i),' > 5000 years.'])
        disp('     FRI distribution not characterised.')
        text(0.25, 0.5, 'FRIs > 5000 yr'); continue
    end

    H = bar(FRI_bin_c, FRI_freq/sum(FRI_freq));
    set(H,'FaceColor',[.75 .75 .75],'BarWidth',1.0); hold on

    if length(FRIz) > 4
        passKS = (length(FRIz) <  30 && pKS > 0.10) || ...
                 (length(FRIz) >= 30 && pKS > 0.05);
        if passKS
            plot(1:1000, wbl_est,'k','linewidth',2)
            text(max(x_lim), max(y_lim), ...
                {['Wbl\it b \rm = ',num2str(round(param(1))),' (', ...
                   num2str(round(wbl_a_ci(1))),'-',num2str(round(wbl_a_ci(2))),')'], ...
                 ['Wbl\it c \rm = ',num2str(round(param(2)*100)/100),' (', ...
                   num2str(round(wbl_b_ci(1)*100)/100),'-', ...
                   num2str(round(wbl_b_ci(2)*100)/100),')'], ...
                 ['mFRI = ',num2str(round(mFRIz)),' (', ...
                   num2str(round(mFRI_ci_z(1))),'-',num2str(round(mFRI_ci_z(2))),')'], ...
                 ['N_F_R_I = ',num2str(length(xPlotZ)-1)]}, ...
                'HorizontalAlignment','right','VerticalAlignment','top','FontSize',FS)
        end
    else
        nStr = num2str(max(0, length(xPlotZ)-1));
        text(max(x_lim), 0.8*max(y_lim), ['N_F_R_I = ',nStr], ...
            'HorizontalAlignment','right','FontSize',FS)
    end

    xlabel('FRI (yr)','FontSize',FS)
    if i == nZones
        ylabel(['proportion OR density (x',num2str(binWidth),')'],'FontSize',FS)
        text(-300, mean(y_lim), char(site),'fontsize',FS,'FontWeight','Bold', ...
            'Rotation',90,'HorizontalAlignment','Center')
    end
    xlim(x_lim); ylim(y_lim)
    set(gca,'TickDir','out','XTick',0:200:x_lim(2),'XMinorTick','on','FontSize',FS)
    box off
    title(char(zoneText(i)),'FontSize',FS,'FontWeight',FW)
end

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 7 — Continuous fire history
%%═══════════════════════════════════════════════════════════════════════════
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(7); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Continuous fire history: peak magnitude, FRIs through time, and smoothed fire frequency')

% ── Panel: Peak magnitude ─────────────────────────────────────────────────
subplot(3,5,1:4); hold off
x = Charcoal.ybpI;
H = bar(x, peak_mag);
set(H,'FaceColor','k','BarWidth',1.0); hold on
plot(x(peakScreenIn), ones(size(peakScreenIn))*0.8*max(peak_mag), ...
    '.','color',[.75 .75 .75])
plot(x(peakIn), 0.8*max(peak_mag), '+r')
y_lim7 = [0, 1.1*max(peak_mag)];
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
    'TickDir','out','XTick',0:1000:max(zoneDiv),'XTickLabel',[], ...
    'XLim',[zoneDiv(1) zoneDiv(end)],'Ylim',y_lim7,'box','off', ...
    'position',[0.1 .75 wm*range(zoneDiv) .12])
if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(y_lim7)*0.90 max(y_lim7)], ...
                '-k','color',[.5 .5 .5],'linewidth',2)
        end
        text(mean(zoneDiv(z-1:z)), max(y_lim7)*0.95, char(zoneText(z-1)), ...
            'horizontalalignment','center','fontweight','normal','fontsize',FS)
    end
end
if nGaps > 0; drawGaps(gca, gapIn, Charcoal); end
ylabel({'peak magnitude','(pieces cm^-^2 peak^-^1)'},'FontSize',FS)
title({'Peak magnitude, FRIs, and fire frequ.','',char(site),''}, ...
    'FontSize',FS,'FontWeight','Bold','VerticalAlignment','Bottom')

% ── Panel: Smoothed fire frequency ───────────────────────────────────────
subplot(3,5,9); hold off
plot(Charcoal.ybpI, ff_sm,'k','linewidth',1); hold on
y_lim9 = [0, max(1.1*max(ff_sm), eps)];
if length(zoneDiv) > 2
    for z = 2:length(zoneDiv)
        if z < length(zoneDiv)
            plot([zoneDiv(z) zoneDiv(z)],[max(y_lim9)*0.90 max(y_lim9)], ...
                '-k','color',[.5 .5 .5],'linewidth',2)
        end
    end
end
if nGaps > 0; drawGaps(gca, gapIn, Charcoal); end
set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
    'TickDir','out','XTick',0:1000:max(zoneDiv), ...
    'XTickLabel',0:max(zoneDiv)/1000,'XLim',[zoneDiv(1) zoneDiv(end)], ...
    'box','off','ylim',y_lim9,'position',[0.1 .45 wm*range(zoneDiv) .12])
ylabel({'fire frequency',['(fires ',num2str(PeakAnalysis.peakFrequ),' yr^-^1)']})
xlabel('time (cal. yr BP x 1000)','FontSize',FS)

% ── Panel: FRIs through time ──────────────────────────────────────────────
subplot(3,5,11:14); hold off
if isnumeric(FRIyr) && ~isnan(FRIyr(1)) && length(smFRI) > 2
    x_lim7 = [0, 1.1*max(FRI)];
    X = [smFRIyr'; flipud(smFRIyr')];
    Y = [smFRIci(:,1); flipud(smFRIci(:,2))];
    plot(FRIyr, FRI,'sk','MarkerFaceColor',[.75 .75 .75],'MarkerSize',4); hold on
    plot(smFRIyr, smFRI,'-k')
    Hf = fill(X, Y,[0.8 0.8 0.8]);
    set(Hf,'edgecolor',[0.8 0.8 0.8])
    plot(FRIyr, FRI,'sk','MarkerFaceColor',[.75 .75 .75],'MarkerSize',4)
    plot(smFRIyr, smFRI,'-k')
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
        'TickDir','out','XTick',0:1000:max(zoneDiv),'Ylim',x_lim7, ...
        'XTickLabel',[],'XLim',[zoneDiv(1) zoneDiv(end)],'box','off', ...
        'position',[0.1 .6 wm*range(zoneDiv) .12],'yScale','linear')
    if nGaps > 0; drawGaps(gca, gapIn, Charcoal); end
    if length(zoneDiv) > 2
        for z = 2:length(zoneDiv)
            if z < length(zoneDiv)
                plot([zoneDiv(z) zoneDiv(z)],[max(x_lim7)*0.90 max(x_lim7)], ...
                    '-k','color',[.5 .5 .5],'linewidth',2)
            end
        end
    end
    ylabel({['FRI (yr fire^-^1)'], ...
        [num2str(PeakAnalysis.peakFrequ),'-yr mean'], ...
        [num2str((1-alpha)*100),'% CI']},'FontSize',FS)
    grid on
end

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 8 — Between-zone comparisons of raw CHAR distributions
%%═══════════════════════════════════════════════════════════════════════════
figPosition = figPosition - [0.0165 0.0225 0 0];
figure(8); clf
set(gcf,'color','w','units','normalized','position',figPosition, ...
    'name','Between-zone comparisons of raw charcoal distributions')

col = [0 0 0; .5 .5 .5; 0 1 0; 0 0 1; 1 0 1; 1 0 0; 1 0.4 0; 0 0.5 0];
nZones = length(zoneDiv)-1;
Charcoal.accZone = NaN(1000, nZones);
for i = 1:nZones
    CHAR_z = Charcoal.acc(Charcoal.ybp >= zoneDiv(i) & ...
                           Charcoal.ybp <  zoneDiv(i+1));
    Charcoal.accZone(1:length(CHAR_z), i) = CHAR_z;
    subplot(3,5,[6 7 8 11 12 13])
    [f, xc] = ecdf(CHAR_z);
    plot(xc, f,'color',col(i,:),'linewidth',2); hold on
end
set(gca,'tickdir','out')
xlabel('CHAR (pieces cm^-^2 yr^-^1)')
ylabel('cum. prop.')
title('CDFs of zone-specific raw CHAR, with KS test results')
box off
legend(char(zoneText))

if nZones > 1
    h_ks = NaN(nZones-1, nZones-1);
    pKS  = NaN(nZones-1, nZones-1);
    for i = 1:nZones-1
        for j = 2:nZones
            [h_ks(i,j-1), pKS(i,j-1)] = kstest2( ...
                Charcoal.accZone(:,i), Charcoal.accZone(:,j));
        end
    end
    pKSResults = zeros(size(pKS)+1);
    pKSResults(2:end,2:end) = pKS;
    pKSResults(1,2:end) = 2:length(pKS)+1;
    pKSResults(2:end,1) = (1:length(pKS))';
    text(0.5*max(Charcoal.acc),0.63,'KS p-value matrix:','FontWeight','bold')
    text(0.45*max(Charcoal.acc),0.48,num2str(round(pKSResults*1000)/1000))
    text(0.45*max(Charcoal.acc),0.58,'Zone','FontWeight','Bold', ...
        'BackgroundColor','w','HorizontalAlignment','Center')
end
text(0,1.5,[char(site),': Between-zone comparisons of raw CHAR distributions'], ...
    'FontSize',FS+2,'FontWeight','Bold')

subplot(3,5,[9 10 14 15])
boxplot(fliplr(Charcoal.accZone),'colors','kkk','symbol','.k')
set(gca,'xticklabel',fliplr(1:nZones),'tickdir','out', ...
    'ylim',[min(min(Charcoal.accZone)) max(max(Charcoal.accZone))])
ylabel('CHAR (pieces cm^-^2 yr^-^1)')
xlabel('zone number')
title('Box plots of raw CHAR per zone')
box off

%%═══════════════════════════════════════════════════════════════════════════
%% FIGURE 9 — Alternative displays of threshold value(s) (allFigures only)
%%═══════════════════════════════════════════════════════════════════════════
if Results.allFigures == 1
    figPosition = figPosition - [0.0165 0.0225 0 0];
    figure(9); clf
    set(gcf,'color','w','units','normalized','position',figPosition, ...
        'name','Alternative displays of threshold value(s)')

    CHAR_ratio     = Charcoal.accI ./ Charcoal.accIS;
    CHAR_residuals = Charcoal.accI  - Charcoal.accIS;
    y_peaks = CharcoalCharPeaks(:,end);  y_peaks(y_peaks==0) = -99;

    if PeakAnalysis.cPeak ~= -99
        if PeakAnalysis.cPeak == 1
            y1 = Charcoal.accIS + CharThresh.pos(:,end);
            y2 = y1 ./ Charcoal.accIS;
            y3 = CharThresh.pos(:,end);
        else  % ratio
            y1 = Charcoal.accIS .* CharThresh.pos(:,end);
            y2 = CharThresh.pos(:,end);
            y3 = Charcoal.accIS - (Charcoal.accIS ./ CharThresh.pos(:,end));
        end
    else
        y1 = CharThresh.pos(:,end);
        y2 = CharThresh.pos(:,end) ./ Charcoal.accIS;
        y3 = CharThresh.pos(:,end) - Charcoal.accIS;
    end

    subplot(3,1,1)
    stairs(Charcoal.ybpI, Charcoal.accI,'color',[.75 .75 .75]); hold on
    plot(Charcoal.ybpI, Charcoal.accIS,'linewidth',1.5,'color','k')
    plot(Charcoal.ybpI, y1,'r','linewidth',1.5)
    plot(Charcoal.ybpI, max(Charcoal.accI)*0.80*y_peaks, '+k','color','r')
    plot(Charcoal.ybpI(peakScreenIn), ...
        max(Charcoal.accI)*0.80*logical(peakScreenIn), '.','color',[.75 .75 .75])
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)]); ylim([0 max(Charcoal.accI)])
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
        'TickDir','out','XTick',0:1000:max(zoneDiv), ...
        'XTickLabel',0:max(zoneDiv)/1000,'position',[0.1 .75 wm*range(zoneDiv) .175])
    box off; grid on
    ylabel('CHAR (pieces cm^-^2 yr^-^1)')
    title('(a) Relationship between C_b_a_c_k_g_r_o_u_n_d and peak threshold')
    legend('C_i_n_t_e_r_p_o_l_a_t_e_d','C_b_a_c_k_g_r_o_u_n_d', ...
        'C_t_h_r_e_s_h_o_l_d','peaks','Location','Eastoutside')

    subplot(3,1,2)
    stairs(Charcoal.ybpI, CHAR_ratio,'color',[.75 .75 .75]); hold on
    plot(Charcoal.ybpI, y2,'r','linewidth',1.5)
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
        'TickDir','out','XTick',0:1000:max(zoneDiv), ...
        'XTickLabel',0:max(zoneDiv)/1000,'position',[0.1 .50 wm*range(zoneDiv) .175])
    box off; grid on
    ylabel({'threshold','(threshold / low-frequency CHAR)'})
    title('(b) Peak threshold as ratios of C_b_a_c_k_g_r_o_u_n_d')
    legend('C_p_e_a_k','C_t_h_r_e_s_h_o_l_d','peaks identified','Location','Eastoutside')

    subplot(3,1,3)
    stairs(Charcoal.ybpI, CHAR_residuals,'color',[.75 .75 .75]); hold on
    plot(Charcoal.ybpI, y3,'r','linewidth',1.5)
    xlim([min(Charcoal.ybpI) max(Charcoal.ybpI)])
    if min(y3) <= 0; ylim([0 abs(max(y3)*2)]); else; ylim([0 max(y3)*2]); end
    set(gca,'FontSize',FS,'XDir','reverse','XMinorTick','on','YMinorTick','on', ...
        'TickDir','out','XTick',0:1000:max(zoneDiv), ...
        'XTickLabel',0:max(zoneDiv)/1000,'position',[0.1 .20 wm*range(zoneDiv) .175])
    box off; grid on
    xlabel('time (cal. yr BP x 1000)')
    ylabel({'threshold','(threshold - low-frequency CHAR)'})
    title('(c) Peak threshold as residuals of C_b_a_c_k_g_r_o_u_n_d')
    legend('C_p_e_a_k','C_t_h_r_e_s_h_o_l_d','peaks identified','Location','Eastoutside')
end

%%═══════════════════════════════════════════════════════════════════════════
%% SAVE FIGURES
%%═══════════════════════════════════════════════════════════════════════════
if Results.saveFigures == 1
    if Results.allFigures == 1
        saveFig(1,'01_pretreatment')
        saveFig(2,'02_threshold_determination')
    end
    saveFig(3,'03_CHAR_analysis')
    saveFig(4,'04_CHAR_peak_sens')
    saveFig(5,'05_cum_peaks_through_time')
    saveFig(6,'06_FRI_dists')
    saveFig(7,'07_continuous_fire_hx')
    saveFig(8,'08_CHAR_dists')
    if Results.allFigures == 1
        saveFig(9,'09_threshold_details')
    end
end

end % ── end CharPlotResults ───────────────────────────────────────────────

%% ── Local helper: save one figure as PDF + TIFF ─────────────────────────
function saveFig(figNum, baseName)
    figure(figNum)
    set(gcf,'PaperPositionMode','auto','PaperType','uslegal')
    orient(gcf,'landscape')
    print('-dpdf','-r300',[baseName,'.pdf'])
    print('-dtiff','-r300',[baseName,'.tif'])
end