function [Charcoal, Post] = CharPostProcess(Charcoal, Pretreatment, ...
    PeakAnalysis, CharThresh, Smoothing, Results, fileName, gapIn, site)
% CharPostProcess   Pure-computation post-processing for CharAnalysis.
%   [Charcoal, Post] = CharPostProcess(Charcoal, Pretreatment, ...
%       PeakAnalysis, CharThresh, Smoothing, Results, fileName, gapIn)
%
%   Performs all analytical calculations that follow peak identification:
%     - Resolves peak/threshold index vectors (peakIn, peakScreenIn, etc.)
%     - Computes peak magnitudes
%     - Computes fire frequency time series
%     - Calls smoothFRI for smoothed FRI curve and CIs
%     - Fits Weibull models and bootstrapped CIs per zone (Figure 6 data)
%     - Assembles the charResults output matrix
%     - Writes output to file if Results.save == 1
%
%   No figure calls are made here. All computed values are returned in the
%   Charcoal struct (fields added) and in the Post struct (plotting inputs).
%
%   Post fields
%   -----------
%   Post.peakIn          - indices of final peaks in Charcoal.ybpI
%   Post.peakScreenIn    - indices of peaks failing minCountP screen
%   Post.CharcoalCharPeaks - [N x 4] peak flag matrix (resolved threshold)
%   Post.threshIn        - [1 x 4] column indices (global threshold only)
%   Post.peak_mag        - [N x 1] peak magnitude time series
%   Post.ff_sm           - [N x 1] smoothed fire frequency
%   Post.FRIyr           - raw FRI years (from smoothFRI)
%   Post.FRI             - raw FRI values (from smoothFRI)
%   Post.smFRIyr         - smoothed FRI years
%   Post.smFRI           - smoothed FRI values
%   Post.smFRIci         - [M x 2] smoothed FRI confidence intervals
%   Post.yis             - smoothed FRI resampled to Charcoal.ybpI grid
%   Post.FRI_params_zone - [nZones x 10] Weibull/FRI params per zone
%   Post.alpha           - alpha used for CIs
%   Post.nBoot           - nBoot used for bootstrapping
%   Post.charResults     - assembled output matrix for file save

%% ?? Parameters ??????????????????????????????????????????????????????????
zoneDiv  = Pretreatment.zoneDiv;
r        = Pretreatment.yrInterp;
alpha    = 0.05;
nBoot    = 100;

%% ?? 1. Resolve peak/threshold index vectors ?????????????????????????????
% Determine which column of Charcoal.charPeaks corresponds to each of the
% four threshold values, and identify final peaks and screened-out peaks.

if PeakAnalysis.threshType == 1        % Global threshold
    threshIn1 = sum(Charcoal.charPeaksThresh) ./ sum(Charcoal.charPeaks);
    threshIn(1) = min(find(threshIn1 >= CharThresh.pos(1,1)));
    threshIn(2) = min(find(threshIn1 >= CharThresh.pos(1,2)));
    threshIn(3) = min(find(threshIn1 >= CharThresh.pos(1,3)));
    threshIn(4) = min(find(threshIn1 >= CharThresh.pos(1,4)));
    peakIn           = find(Charcoal.charPeaks(:, threshIn(end)));
    CharcoalCharPeaks = Charcoal.charPeaks(:, threshIn);
    peakScreenIn     = find(CharThresh.minCountP(:, threshIn(end)) > ...
                            PeakAnalysis.minCountP);
else                                   % Local threshold
    threshIn         = [];             % Not used for local threshold
    peakIn           = find(Charcoal.charPeaks(:, end));
    CharcoalCharPeaks = Charcoal.charPeaks;
    peakScreenIn     = find(CharThresh.minCountP(:, end) > ...
                            PeakAnalysis.minCountP);
end

%% ?? 2. Peak magnitude ???????????????????????????????????????????????????
% For each contiguous run of samples where peak CHAR exceeds the final
% threshold, accumulate (sum * resolution) to get pieces cm^-2 peak^-1.

c_peaks        = Charcoal.peak - CharThresh.pos(:, end);
c_peaks(c_peaks < 0) = 0;
c_peaks(end)   = 0;                    % First sample never a peak

peak_in  = zeros(length(c_peaks), 2);
peak_mag = zeros(length(c_peaks), 1);

for in = 1:length(c_peaks)
    if in == 1 && c_peaks(in) > 0
        peak_in(in,1) = in;
        step = 1;
        while c_peaks(in+step) > 0
            peak_in(in,2) = in + step;
            step = step + 1;
        end
    else
        if c_peaks(in) > 0 && c_peaks(in-1) == 0
            peak_in(in,1) = in;
            step = 1;
            while c_peaks(in+step) > 0
                peak_in(in,2) = in + step;
                step = step + 1;
            end
        end
    end
    if peak_in(in,1) > 0 && peak_in(in,2) == 0
        peak_in(in,2) = peak_in(in,1);
    end
    if peak_in(in,2) > 0
        peak_mag(peak_in(in,2), 1) = ...
            sum(c_peaks(peak_in(in,1):peak_in(in,2))) * ...
            (diff(peak_in(in,:)) + r);
    end
end

%% ?? 3. Fire frequency time series ???????????????????????????????????????
% Counts fires in a sliding window of peakFrequ years, scaled to the
% window actually used at record edges, then smoothed with lowess.

ff_sum_yr = PeakAnalysis.peakFrequ;
ff_sm_yr  = PeakAnalysis.peakFrequ;

peaksFrequ = zeros(size(Charcoal.ybpI));

for i = 1:length(Charcoal.ybpI)
    halfWin = floor(ff_sum_yr/r) / 2;
    if i < halfWin                          % Start of record
        winLen = halfWin + i;
        peaksFrequ(i) = sum(CharcoalCharPeaks(1:floor(winLen), end)) * ...
                        (round(ff_sum_yr/r) / floor(winLen));
    elseif i > length(Charcoal.ybpI) - halfWin   % End of record
        winSlice = CharcoalCharPeaks(i - floor(halfWin):end, end);
        peaksFrequ(i) = sum(winSlice) * ...
                        ((ff_sum_yr/r) / length(winSlice));
    else                                    % Middle of record
        idxLo = max(1,   ceil(i - 0.5*(ff_sum_yr/r)) + 1);
        idxHi = min(length(Charcoal.ybpI), ceil(i + 0.5*round(ff_sum_yr/r)));
        peaksFrequ(i) = sum(CharcoalCharPeaks(idxLo:idxHi, end));
    end
end

ff_sm = charLowess(peaksFrequ, ff_sm_yr/r);

%% ?? 4. Smoothed FRI curve ????????????????????????????????????????????????
peak_yrs = Charcoal.ybpI(CharcoalCharPeaks(:,end) > 0);

if length(peak_yrs) > 2
    [FRIyr, FRI, smFRIyr, smFRI, smFRIci] = smoothFRI( ...
        Charcoal.ybpI, CharcoalCharPeaks(:,end), ...
        PeakAnalysis.peakFrequ, alpha, nBoot, 1, 0);
    if length(smFRI) > 2
        yis = interp1(smFRIyr, smFRI, ...
                      Charcoal.ybpI(Charcoal.ybpI < max(smFRIyr)));
    else
        yis = -999 + zeros(size(smFRI));
        warning('CharPostProcess: fewer than 3 FRIs — smoothing set to -999.')
    end
else
    FRIyr  = NaN;  FRI    = NaN;
    smFRIyr = NaN; smFRI  = NaN;  smFRIci = NaN;
    yis    = NaN;
end

%% ?? 5. Per-zone FRI statistics (Figure 6 data) ??????????????????????????
% FRI_params_zone columns:
%   1  nFRI
%   2  mFRI
%   3  mFRI_uCI
%   4  mFRI_lCI
%   5  WblB (scale, a)
%   6  WblB_uCI
%   7  WblB_lCI
%   8  WblC (shape, b)
%   9  WblC_uCI
%   10 WblC_lCI

nZones = length(zoneDiv) - 1;
FRI_params_zone = -999 * ones(nZones, 10);
binWidth = 20;

for i = 1:nZones
    xPlot = Charcoal.ybpI(CharcoalCharPeaks(:,end) > 0);
    xPlot = xPlot(xPlot >= zoneDiv(i) & xPlot < zoneDiv(i+1));
    FRIz  = diff(xPlot);

    if max(FRIz) > 5000
        % FRIs too long to characterise — leave as -999
        continue
    end
    if length(FRIz) <= 1
        continue
    end

    [FRI_freq, FRI_bin] = histcounts(FRIz, binWidth:binWidth:1000, ...
                                     'Normalization','count');
    FRI_bin_centers = FRI_bin(1:end-1) + binWidth/2;

    % Weibull fit to binned data
    param = wblfit(FRI_bin_centers, [], [], FRI_freq);

    % One-sample KS goodness-of-fit against fitted Weibull CDF
    FRIBinKS  = 0:20:5000;
    wbl_cdf   = wblcdf(FRIBinKS, param(1), param(2));
    [~, pKS]  = kstest(FRIz, [FRIBinKS', wbl_cdf']);

    % Bootstrap CIs on Weibull parameters and mFRI
    mean_mFRI  = NaN(nBoot, 1);
    param_t    = NaN(nBoot, 2);
    for t = 1:nBoot
        FRI_t            = FRIz(randi(length(FRIz), length(FRIz), 1));
        [f_t, ~]         = histcounts(FRI_t, FRI_bin, 'Normalization','count');
        param_t(t,:)     = wblfit(FRI_bin_centers, [], [], f_t);
        mean_mFRI(t)     = mean(FRI_t);
    end
    wbl_a_ci       = prctile(param_t(:,1), [2.5 97.5]);
    wbl_b_ci       = prctile(param_t(:,2), [2.5 97.5]);
    mean_mFRI_ci   = prctile(mean_mFRI,    [2.5 97.5]);

    FRI_params_zone(i,:) = [length(FRIz), mean(FRIz), mean_mFRI_ci, ...
                             param(1), wbl_a_ci, param(2), wbl_b_ci];

    % Attach per-zone KS result and fitted curve to Post for plotting
    Post.zone(i).FRI             = FRIz;
    Post.zone(i).FRI_bin_centers = FRI_bin_centers;
    Post.zone(i).FRI_freq        = FRI_freq;
    Post.zone(i).param           = param;
    Post.zone(i).pKS             = pKS;
    Post.zone(i).wbl_est         = binWidth * wblpdf(1:1000, param(1), param(2));
    Post.zone(i).wbl_a_ci        = wbl_a_ci;
    Post.zone(i).wbl_b_ci        = wbl_b_ci;
    Post.zone(i).mean_mFRI       = mean(FRIz);
    Post.zone(i).mean_mFRI_ci    = mean_mFRI_ci;
    Post.zone(i).xPlot           = xPlot;
end

%% ?? 6. Update Charcoal struct ???????????????????????????????????????????
Charcoal.peakInsig        = zeros(size(Charcoal.ybpI));
Charcoal.peakInsig(peakScreenIn) = 1;
Charcoal.peakMagnitude    = peak_mag;
Charcoal.smoothedFRI      = smFRI;
Charcoal.smoothedFireFrequ = ff_sm;
Charcoal.peaksFrequ       = peaksFrequ;

%% ?? 7. Assemble output matrix ???????????????????????????????????????????
N = length(Charcoal.cmI);
charResults = NaN * ones(N, 33);
charResults(:, 1:22) = [Charcoal.cmI,  Charcoal.ybpI,  Charcoal.countI, ...
    Charcoal.volI,   Charcoal.conI,    Charcoal.accI,  Charcoal.accIS, ...
    Charcoal.peak,   CharThresh.pos,   CharThresh.neg(:,end), ...
    CharThresh.SNI,  CharThresh.GOF(:,1), CharcoalCharPeaks, ...
    Charcoal.peakInsig, Charcoal.peakMagnitude, Charcoal.smoothedFireFrequ];
if ~isnan(yis)
    charResults(1:length(yis), 23) = yis;
end
charResults(1:nZones, 24:33) = FRI_params_zone;

%% ?? 8. Save output ??????????????????????????????????????????????????????
if Results.save == 1
    if min(fileName(end-2:end) == 'xls') > 0
        writematrix(charResults, fileName, 'Sheet', 'charResults', ...
                    'Range', 'A2');
    else
        outputResults(charResults, fileName, site);
    end
end

%% ?? 9. Pack Post struct ?????????????????????????????????????????????????
Post.peakIn           = peakIn;
Post.peakScreenIn     = peakScreenIn;
Post.CharcoalCharPeaks = CharcoalCharPeaks;
Post.threshIn         = threshIn;
Post.peak_mag         = peak_mag;
Post.ff_sm            = ff_sm;
Post.FRIyr            = FRIyr;
Post.FRI              = FRI;
Post.smFRIyr          = smFRIyr;
Post.smFRI            = smFRI;
Post.smFRIci          = smFRIci;
Post.yis              = yis;
Post.FRI_params_zone  = FRI_params_zone;
Post.alpha            = alpha;
Post.nBoot            = nBoot;
Post.charResults      = charResults;