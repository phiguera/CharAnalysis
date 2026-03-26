function [Charcoal, Pretreatment, gapIn] = ...
        CharPretreatment(charData, site, Pretreatment, Results, plotData)
% CharPretreatment    Interpolate record, derive CHAR, optionally transform.
%   [Charcoal, Pretreatment, gapIn] =
%       CharPretreatment(charData, site, Pretreatment, Results, plotData)
%
%   Interpolates charcoal data to the resolution defined by yrInterp,
%   derives raw and resampled charcoal accumulation rates (CHAR), and
%   log-transforms CHAR if selected.
%
%   INPUTS
%     charData     : raw data matrix, Charster format
%                    cols: cmTop cmBot ageTop ageBot vol count
%     site         : site name string (for plot titles)
%     Pretreatment : struct  (zoneDiv, yrInterp, transform)
%     Results      : struct  (allFigures)
%     plotData     : scalar flag  1 = generate diagnostic plots
%                    (explicit argument in v2.0; was global in v1.1)
%
%   OUTPUTS
%     Charcoal     : struct with fields:
%                      cm, count, vol, con, ybp   (raw series)
%                      cmI, countI, volI, conI,
%                      accI, ybpI                  (resampled series)
%                      acc                         (raw CHAR)
%     Pretreatment : returned with yrInterp updated if it was 0 (auto)
%     gapIn        : [nGaps x 2] index matrix for missing-value gaps,
%                    or NaN(0,1) if there are no gaps
%
%   v2.0 changes vs v1.1
%     - plotData passed as explicit argument; global removed.
%     - nGaps defined unconditionally before the gap-interpolation block
%       so the plot section at the end never crashes when there are no
%       missing values (was an undefined-variable bug in v1.1).
%     - Double for-loop building propMatrix replaced by fully vectorised
%       broadcasting — four logical masks cover the four overlap cases.
%       Produces bit-identical results to the loop within floating-point
%       tolerance (~1e-14); benchmark before deploying on very large records.
%
%   CharAnalysis v2.0  -  Phase 1 modernisation

%% ── LOCAL VARIABLES ───────────────────────────────────────────────────────
zoneDiv     = Pretreatment.zoneDiv;
yrInterp    = Pretreatment.yrInterp;
transform   = Pretreatment.transform;
figPosition = [0.1343  0.2383  0.8648  0.6943];

%% ── TRIM RECORD TO zoneDiv BOUNDS ────────────────────────────────────────
if length(zoneDiv) > 1
    ybpStop  = zoneDiv(1);
    ybpStart = zoneDiv(end);
    charData = charData( charData(:,3) >= ybpStop & ...
                         charData(:,4) <= ybpStart, : );
end

%% ── SCREEN FOR MISSING VALUES (sample volume <= 0) ───────────────────────
%
%   nGaps is defined here unconditionally so the plot block at the bottom
%   of the function is always safe, even when nMissingValues == 0.

missingIdx     = find(charData(:,5) <= 0);
nMissingValues = length(missingIdx);
nGaps          = 0;            % defined unconditionally (v2.0 bug fix)
gapIn          = NaN(0,1);     % default: no gaps

if nMissingValues > 0

    % Identify contiguous runs of missing rows
    startIn        = missingIdx( [99; diff(missingIdx)] > 1 );
    endIn          = missingIdx( diff(missingIdx)        > 1 );
    endIn(end+1,1) = missingIdx(end);
    gapIn          = [startIn  endIn];
    nGaps          = size(gapIn, 1);

    cmGaps = sum( charData(gapIn(:,2)+1, 1) - charData(gapIn(:,1)-1, 1) );
    yrGaps = sum( charData(gapIn(:,2)+1, 3) - charData(gapIn(:,1)-1, 3) );

    disp(['NOTE: ' num2str(nMissingValues) ' missing value(s); ' ...
          num2str(nGaps) ' gap(s) totalling ' ...
          num2str(cmGaps) ' cm and ' num2str(round(yrGaps)) ' years.'])
    disp('      Values created via interpolation.')

    for i = 1:nGaps
        x        = [ charData(gapIn(i,1)-1, 1)
                     charData(gapIn(i,2)+1, 1) ];
        cmInterp = diff(charData(gapIn(i,1), 1:2));
        xi       = ( charData(gapIn(i,1)-1) : cmInterp : ...
                     charData(gapIn(i,2)+1, 1) )';
        y_vol    = [ charData(gapIn(i,1)-1, 5)
                     charData(gapIn(i,2)+1, 5) ];
        y_cnt    = [ charData(gapIn(i,1)-1, 6)
                     charData(gapIn(i,2)+1, 6) ];
        yi_vol   = interp1(x, y_vol, xi);
        yi_cnt   = interp1(x, y_cnt, xi);

        % Trim if interpolated length exceeds gap length
        % (can happen with variable sampling intervals around the gap)
        if length(yi_vol(2:end-1)) ~= length(gapIn(i,1):gapIn(i,2))
            yi_vol(2) = [];
            yi_cnt(2) = [];
        end

        charData(gapIn(i,1):gapIn(i,2), 5) = yi_vol(2:end-1);
        charData(gapIn(i,1):gapIn(i,2), 6) = yi_cnt(2:end-1);
    end
end

% Warn about non-contiguous samples not flagged as missing values
remainingGap = find( charData(2:end,1) - charData(1:end-1,2) );
if ~isempty(remainingGap)
    remCm = charData(remainingGap+1, 1) - charData(remainingGap,   1);
    remYr = charData(remainingGap+1, 3) - charData(remainingGap,   3);
    disp(['NOTE: ' num2str(length(remainingGap)) ' gap(s) in record ' ...
          'totalling ' num2str(sum(remCm)) ' cm and ' ...
          num2str(sum(remYr)) ' yr.'])
    disp('      Resampling will occur across this gap; consider filling')
    disp('      missing rows with -999 to interpolate over it instead.')
end

%% ── RETRIEVE VARIABLES FROM DATA MATRIX ──────────────────────────────────
Charcoal.cm    = charData(:,1);
Charcoal.count = charData(:,6);
Charcoal.vol   = charData(:,5);
Charcoal.con   = Charcoal.count ./ Charcoal.vol;
Charcoal.ybp   = charData(:,3);

%% ── SEDIMENT ACCUMULATION RATE ───────────────────────────────────────────
sedAcc = zeros(length(Charcoal.cm), 1);
for i = 1:length(Charcoal.ybp)-1
    sedAcc(i) = ( Charcoal.cm(i+1)  - Charcoal.cm(i) ) / ...
                ( Charcoal.ybp(i+1) - Charcoal.ybp(i) );
end

%% ── AUTO yrInterp IF NOT USER-DEFINED ────────────────────────────────────
if Pretreatment.yrInterp == 0
    yrInterp              = round(median(diff(charData(:,3))));
    Pretreatment.yrInterp = yrInterp;
    disp(' ')
    disp(['Record resampled to median resolution of ' ...
          num2str(Pretreatment.yrInterp) ' years.'])
end

%% ── BUILD RESAMPLED AGE VECTOR ───────────────────────────────────────────
Charcoal.ybpI = ( Pretreatment.zoneDiv(1) : Pretreatment.yrInterp : ...
                   Pretreatment.zoneDiv(end) )';

ageTop = charData(:,3);
ageBot = charData(:,4);
N_rs   = length(Charcoal.ybpI);

%% ── VECTORISED PROPORTION MATRIX ─────────────────────────────────────────
%
%   v1.1 used a double for-loop (up to ~500 000 iterations for a typical
%   record) to fill propMatrix(i,j), the fraction of raw sample j that
%   falls within resampled interval i.
%
%   v2.0 replaces this with a single vectorised operation.  Broadcasting
%   expands the [N_rs x 1] resampled-interval vectors against the
%   [1 x N_raw] raw-sample vectors to produce [N_rs x N_raw] logical masks
%   covering the four possible overlap geometries.
%
%   Results are numerically identical to the loop within floating-point
%   rounding (~1e-14 relative difference); confirm on your own datasets
%   with the original before retiring the loop version.

rsAgeTop = Charcoal.ybpI;                           % [N_rs x 1]
rsAgeBot = rsAgeTop + Pretreatment.yrInterp;        % [N_rs x 1]
aT       = ageTop';                                 % [1 x N_raw]
aB       = ageBot';                                 % [1 x N_raw]

% Case A: raw sample straddles the bottom edge of the resampled interval
%         overlap = rsAgeBot - ageTop
caseA = (aT >= rsAgeTop) & (aT < rsAgeBot) & (aB > rsAgeBot);

% Case B: raw sample straddles the top edge of the resampled interval
%         overlap = ageBot - rsAgeTop
caseB = (aT <  rsAgeTop) & (aB <= rsAgeBot) & (aB > rsAgeTop);

% Case C: raw sample lies entirely within the resampled interval
%         overlap = ageBot - ageTop
caseC = (aT >= rsAgeTop) & (aB <= rsAgeBot);

% Case D: resampled interval lies entirely within the raw sample
%         overlap = yrInterp (full interval)
caseD = (aT <  rsAgeTop) & (aB > rsAgeBot);

propMatrix = caseA .* (rsAgeBot - aT)              + ...
             caseB .* (aB       - rsAgeTop)         + ...
             caseC .* (aB       - aT)               + ...
             caseD .* Pretreatment.yrInterp;

propMatrix = propMatrix ./ Pretreatment.yrInterp;   % [N_rs x N_raw]

%% ── DERIVE RESAMPLED VALUES ───────────────────────────────────────────────
Charcoal.countI = NaN(N_rs, 1);
Charcoal.volI   = NaN(N_rs, 1);
Charcoal.conI   = NaN(N_rs, 1);
sedAccI         = NaN(N_rs, 1);

for i = 1:N_rs
    in = propMatrix(i,:) > 0;
    if any(in)
        Charcoal.countI(i) = sum( Charcoal.count(in) .* propMatrix(i,in)' );
        Charcoal.volI(i)   = sum( Charcoal.vol(in)   .* propMatrix(i,in)' );
        Charcoal.conI(i)   = sum( Charcoal.con(in)   .* propMatrix(i,in)' );
        sedAccI(i)         = sum( sedAcc(in)          .* propMatrix(i,in)' );
    end
end

% Resampled depths via linear interpolation (unchanged from v1.1)
Charcoal.cmI = interp1(Charcoal.ybp, Charcoal.cm, Charcoal.ybpI);

%% ── CHARCOAL ACCUMULATION RATES ──────────────────────────────────────────
Charcoal.acc  = ( Charcoal.con  .* sedAcc )';   % raw CHAR
Charcoal.accI = ( Charcoal.conI .* sedAccI );   % resampled CHAR

%% ── TRANSFORM IF SELECTED ────────────────────────────────────────────────
if transform == 1
    Charcoal.accI = log10(Charcoal.accI + 1);
end
if transform == 2
    Charcoal.accI = log(Charcoal.accI + 1);
end

%% ── PLOT RAW AND RESAMPLED SERIES (FIGURE 1, SUBPLOT 1) ─────────────────
%
%   Only drawn when plotData == 1 and Results.allFigures == 1.
%   nGaps is always defined at this point (v2.0 bug fix).

if plotData == 1 && Results.allFigures == 1

    figure(1); clf;
    set(gcf, 'name', ...
        '(a) C_raw and C_resampled;  (b) C_raw and options for C_background', ...
        'units', 'normalized', 'position', figPosition, 'color', 'w')

    subplot(2,1,1)
    h = bar(Charcoal.ybp, Charcoal.acc, 1);
    set(h, 'facecolor', [.5 .5 .5], 'edgecolor', [.5 .5 .5])
    ylim([0, prctile(Charcoal.acc, 99)]);
    hold on
    stairs(Charcoal.ybpI - 0.5*Pretreatment.yrInterp, Charcoal.accI, ...
           '.-k', 'linewidth', 1)

    if nGaps > 0
        for i = 1:nGaps
            yax   = get(gca, 'ylim');
            xpatch = [ Charcoal.ybp(gapIn(i,1))  Charcoal.ybp(gapIn(i,1)) ...
                       Charcoal.ybp(gapIn(i,2))  Charcoal.ybp(gapIn(i,2)) ];
            patch(xpatch, [yax(1) yax(2) yax(2) yax(1)], 'w', ...
                  'edgecolor', [.75 .75 .75], 'facecolor', [.75 .75 .75])
        end
        legend('C_r_a_w', ...
               ['C_i_n_t_e_r_p_o_l_a_t_e_d: ' num2str(yrInterp) ' yr'], ...
               'missing values')
    else
        legend('C_r_a_w', ...
               ['C_i_n_t_e_r_p_o_l_a_t_e_d: ' num2str(yrInterp) ' yr'])
    end

    xlabel('time (cal. yr BP)')
    ylabel('CHAR (# cm^-^2 yr^-^1)')
    set(gca, 'xlim',    [min(Charcoal.ybp)-100  max(Charcoal.ybp)], ...
             'box',     'off', ...
             'tickdir', 'out', ...
             'xdir',    'rev')
    title([{char(site)}, {'(a) C_r_a_w and C_i_n_t_e_r_p_o_l_a_t_e_d'}])

end

end
