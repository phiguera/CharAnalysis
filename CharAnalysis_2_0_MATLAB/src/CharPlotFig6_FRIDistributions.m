function CharPlotFig6_FRIDistributions(r)
% CharPlotFig6_FRIDistributions   Figure 6: FRI distributions by zone with Weibull models.
%
%   CharPlotFig6_FRIDistributions(results)
%
%   Histogram of FRIs within each zone (20-yr bins). If the fitted Weibull
%   model passes the goodness-of-fit test (p > 0.10 if n < 30; p > 0.05
%   if n >= 30), model parameters and 95% confidence estimates are listed
%   along with mFRI, confidence estimates, and n.
%
%   All analytical values are taken from the Post struct produced by
%   CharPostProcess. No computation is performed here.

%% -- Unpack results struct ---------------------------------------------------
Charcoal     = r.Charcoal;
Pretreatment = r.Pretreatment;
site         = r.site;
Post         = r.Post;

zoneDiv  = Pretreatment.zoneDiv;
nZones   = length(zoneDiv) - 1;

S = charPlotSetup(Charcoal);
FS          = S.FS;
FW          = S.FW;
zoneText    = S.zoneText;
figPosition = S.figPosition - [0.0165*3 0.0225*3 0 0];

%% ------------------------------------------------------------------------
%% FIGURE 6 - FRI distributions by zone with Weibull models
%% ------------------------------------------------------------------------
figure(6); clf
set(gcf, 'color', 'w', 'units', 'normalized', 'position', figPosition, ...
    'name', 'Fire return intervals by zone, with Weibull models if GOF test is passed')

binWidth = 20;
x_lim    = [0 800];
y_lim    = [0 0.35];

for i = 1:nZones
    inPlot = fliplr(1:nZones);
    subplot(3,5,inPlot(i))

    % Retrieve per-zone data from Post (computed by CharPostProcess)
    if ~isfield(Post, 'zone') || ~isfield(Post.zone(i), 'FRI') || ...
            isempty(Post.zone(i).FRI)
        text(0.25, 0.5, 'insufficient data'); box off; continue
    end

    FRIz      = Post.zone(i).FRI;
    FRI_bin_c = Post.zone(i).FRI_bin_centers;
    FRI_freq  = Post.zone(i).FRI_freq;
    param     = Post.zone(i).param;
    pKS       = Post.zone(i).pKS;
    wbl_est   = Post.zone(i).wbl_est;
    wbl_a_ci  = Post.zone(i).wbl_a_ci;
    wbl_b_ci  = Post.zone(i).wbl_b_ci;
    mFRIz     = Post.zone(i).mean_mFRI;
    mFRI_ci_z = Post.zone(i).mean_mFRI_ci;
    xPlotZ    = Post.zone(i).xPlot;

    if max(FRIz) > 5000
        disp(['WARNING, Figure 6: FRIs in Zone ', num2str(i), ' > 5000 years.'])
        disp('     FRI distribution not characterised.')
        text(0.25, 0.5, 'FRIs > 5000 yr'); continue
    end

    H = bar(FRI_bin_c, FRI_freq/sum(FRI_freq));
    set(H, 'FaceColor', [.75 .75 .75], 'BarWidth', 1.0); hold on

    if length(FRIz) > 4
        passKS = (length(FRIz) <  30 && pKS > 0.10) || ...
                 (length(FRIz) >= 30 && pKS > 0.05);
        if passKS
            plot(1:1000, wbl_est, 'k', 'linewidth', 2)
            text(max(x_lim), max(y_lim), ...
                {['Wbl\it b \rm = ', num2str(round(param(1))), ' (', ...
                   num2str(round(wbl_a_ci(1))), '-', num2str(round(wbl_a_ci(2))), ')'], ...
                 ['Wbl\it c \rm = ', num2str(round(param(2)*100)/100), ' (', ...
                   num2str(round(wbl_b_ci(1)*100)/100), '-', ...
                   num2str(round(wbl_b_ci(2)*100)/100), ')'], ...
                 ['mFRI = ', num2str(round(mFRIz)), ' (', ...
                   num2str(round(mFRI_ci_z(1))), '-', num2str(round(mFRI_ci_z(2))), ')'], ...
                 ['N_F_R_I = ', num2str(length(xPlotZ)-1)]}, ...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'top', 'FontSize', FS)
        end
    else
        nStr = num2str(max(0, length(xPlotZ)-1));
        text(max(x_lim), 0.8*max(y_lim), ['N_F_R_I = ', nStr], ...
            'HorizontalAlignment', 'right', 'FontSize', FS)
    end

    xlabel('FRI (yr)', 'FontSize', FS)
    if i == nZones
        ylabel(['proportion OR density (x', num2str(binWidth), ')'], 'FontSize', FS)
        text(-300, mean(y_lim), char(site), 'fontsize', FS, 'FontWeight', 'Bold', ...
            'Rotation', 90, 'HorizontalAlignment', 'Center')
    end
    xlim(x_lim); ylim(y_lim)
    set(gca, 'TickDir', 'out', 'XTick', 0:200:x_lim(2), 'XMinorTick', 'on', 'FontSize', FS)
    box off
    title(char(zoneText(i)), 'FontSize', FS, 'FontWeight', FW)
end

end
