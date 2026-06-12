% z_dump_bkgSens_refs.m
%
% Generate MATLAB v2.0 reference arrays for the R background-sensitivity
% regression test (CharAnalysis_2_0_R/tests/test_bkg_sensitivity.R).
%
% For the Code Lake (CO) dataset it runs the full pipeline (which draws
% Figure 2, needed by the global-threshold contour path) and then calls
% bkgCharSensitivity directly, capturing the three returned arrays:
%   z      peak-count surface/columns
%   GOF_i  goodness-of-fit (KS p-value)
%   SNI_i  signal-to-noise index
%
% Two configurations are written:
%   CO local   threshType = 2  (CO as shipped; 2x2 diagnostic path)
%   CO global  threshType = 1  (temporary variant; contour path)
%
% Output: ./bkgSens_refs/CO_<local|global>_<z|GOF_i|SNI_i>.csv
%
% Requirements: MATLAB R2022a+ (readlines/writelines). Run from the
%   CharAnalysis_2_0_MATLAB directory:
%       >> z_dump_bkgSens_refs
%
% AI assistance (Claude / Anthropic) was used to draft this harness.

clear; close all; clc
addpath('src')

outDir = 'bkgSens_refs';
if ~exist(outDir, 'dir')
    mkdir(outDir)
end

% ---- CO local (threshType = 2, as shipped) ------------------------------
disp('=== CO local (threshType = 2) ===')
dump_one('CO_charParams.csv', 'local', outDir)

% ---- CO global (threshType = 1, temporary variant) ----------------------
disp('=== CO global (threshType = 1) ===')
copyfile('CO_charData.csv', 'CO_global_charData.csv')

T = readlines('CO_charParams.csv');
for k = 1:numel(T)
    if contains(T(k), ',threshType,')
        T(k) = replace(T(k), ',threshType,2,', ',threshType,1,');
    end
end
writelines(T, 'CO_global_charParams.csv')

dump_one('CO_global_charParams.csv', 'global', outDir)

% Clean up the temporary global variant.
delete('CO_global_charParams.csv')
delete('CO_global_charData.csv')

disp(' ')
disp(['Reference arrays written to: ', fullfile(cd, outDir)])
disp('Copy this folder to CharAnalysis_2_0_R/tests/bkgSens_refs/ (or')
disp('point CHAR_MATLAB_REFS in test_bkg_sensitivity.R at it).')

% =========================================================================
% Local function: run one configuration and write the three arrays.
% =========================================================================
function dump_one(paramFile, tag, outDir)
    % 'standard' run draws Figures 1-9 and leaves Figure 2 available for the
    % global contour path inside bkgCharSensitivity.
    R = CharAnalysis(paramFile, 'standard');

    [z, GOF_i, SNI_i] = bkgCharSensitivity(R.Charcoal, R.CharThresh, ...
        R.PeakAnalysis, R.Pretreatment, R.Smoothing, R.Results, R.site);

    writematrix(z,     fullfile(outDir, ['CO_', tag, '_z.csv']))
    writematrix(GOF_i, fullfile(outDir, ['CO_', tag, '_GOF_i.csv']))
    writematrix(SNI_i, fullfile(outDir, ['CO_', tag, '_SNI_i.csv']))

    disp(['  wrote CO_', tag, '_{z,GOF_i,SNI_i}.csv  (z is ', ...
        num2str(size(z,1)), ' x ', num2str(size(z,2)), ')'])
end
