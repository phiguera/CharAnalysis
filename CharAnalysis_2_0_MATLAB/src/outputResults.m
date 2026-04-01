function outputResults(charResults, fileName, site)
% outputResults   Write charResults matrix to a .csv file with headers.
%   outputResults(charResults, fileName, site)
%
%   charResults : [N x 33] numeric matrix of output values
%   fileName    : input parameter file name (used only to detect path separator)
%   site        : site name string (used to build output file name)
%
%   Output file: <current directory>/<site>_charResults.csv
%
%   v2.0: replaces manual cell padding + fprintf with writecell(),
%   which handles quoting, delimiters, and line endings automatically.

%% ?? Column headers ???????????????????????????????????????????????????????
header = { ...
    'cm Top_i (cm)',             'age Top_i (yr BP)', ...
    'char Count_i (#)',          'char Vol_i (cm3)', ...
    'char Con_i (# cm-3)',       'char Acc_i (# cm-2 yr-1)', ...
    'charBkg (# cm-2 yr-1)',     'char Peak (# cm-2 yr-1)', ...
    'thresh 1 (# cm-2 yr-1)',    'thresh 2 (# cm-2 yr-1)', ...
    'thresh 3 (# cm-2 yr-1)',    'thresh FinalPos (# cm-2 yr-1)', ...
    'thresh FinalNeg (# cm-2 yr-1)', 'SNI (index)', ...
    'thresh GOF (p-val)',        'peaks 1', ...
    'peaks 2',                   'peaks 3', ...
    'peaks Final',               'peaks Insig.', ...
    'peak Mag (# cm-2 peak-1)',  'smPeak Frequ (peaks 1ka-1)', ...
    'smFRIs (yr fire-1)',        'nFRIs (#)', ...
    'mFRI (yr fire-1)',          'mFRI_uCI (yr fire-1)', ...
    'mFRI_lCI (yr fire-1)',      'WBLb (yr)', ...
    'WBLb_uCI (yr)',             'WBLb_lCI (yr)', ...
    'WBLc (unitless)',           'WBLc_uCI (unitless)', ...
    'WBLc_lCI (unitless)'};

%% ?? Build output file path ???????????????????????????????????????????????
i1 = regexp(fileName, '\');   % PC path separator
i2 = regexp(fileName, '/');   % Mac/Linux path separator
if length(i1) >= length(i2)
    fileNameResults = [cd '\' site '_charResults.csv'];
else
    fileNameResults = [cd '/' site '_charResults.csv'];
end

%% ?? Combine header and data into one cell array ??????????????????????????
% Convert numeric matrix to cell array of strings, NaN ? empty string.
[nRows, nCols] = size(charResults);
dataCells = cell(nRows, nCols);
for c = 1:nCols
    for r = 1:nRows
        v = charResults(r, c);
        if isnan(v)
            dataCells{r, c} = '';
        else
            dataCells{r, c} = num2str(v);
        end
    end
end

output = [header; dataCells];

%% ?? Write to file ????????????????????????????????????????????????????????
writecell(output, fileNameResults);
disp(['   Results saved to: ' fileNameResults])