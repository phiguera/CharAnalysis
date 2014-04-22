function outputResults(charResults,fileName,site)
% function outputResults(charResults,fileName)
% Place charResults into a .csv file, including column headers.

%% Get number of rows and columns for input data...
numrows = size(charResults,1);  % Number of rows in charResults.
numcols = size(charResults,2);  % Number of columns in charResults.

%% Create header row
header = {'cm Top_i (cm)','age Top_i (yr BP)','char Count_i (#)','char Vol_i (cm3)', ...
    'char Con_i (# cm-3)','char Acc_i (# cm-2 yr-1)','charBkg (# cm-2 yr-1)', ...
    'char Peak (# cm-2 yr-1)','thresh 1 (# cm-2 yr-1)','thresh 2 (# cm-2 yr-1)', ...
    'thresh 3 (# cm-2 yr-1)','thresh FinalPos (# cm-2 yr-1)','thresh FinalNeg (# cm-2 yr-1)', ...
    'SNI (index)','thresh GOF (p-val)','peaks 1','peaks 2','peaks 3','peaks Final', ...
    'peaks Insig.','peak Mag (# cm-2 peak-1)','smPeak Frequ (peaks 1ka-1)', ...
    'smFRIs (yr*fire-1)','nFRIs (#)','mFRI (yr fire-1)','mFRI_uCI (yr fire-1)', ...
    'mFRI_lCI (yr fire-1)','WBLb (yr)','WBLb_uCI (yr)','WBLb_lCI (yr fire-1)', ...
    'WBLc (unitless)','WBLc_uCI (unitless)','WBLc_lCI (unitless)'};

%% Preallocate temporary matrices
temp1 = cell(numrows,numcols);
temp2 = cell(numrows + 1,numcols);

%% Combine data and header into one cell array of strings
for i=1:numcols
    temp1(1,i) = header(1,i);
    for j = 1:numrows
        temp1(j+1,i) = {num2str(charResults(j,i))};
    end
end
numrows = numrows+1;

%% Now convert all the strings to characters, adding spaces to make cell
    % lengths uniform (32 chars long)
for i = 1:numrows
    for j = 1:numcols
        a = length(char(temp1(i,j))); 
        if j == 1 % 1st column is special--it doesn't start with ','
            temp2(i,j) = {[char(temp1(i,1)), char(32*ones(1,32-a))]};
        else
            temp2(i,j) = {[',', char(temp1(i,j)), char(32*ones(1,31-a))]};
        end
    end
end

%% Convert cell array to matrix
output = cell2mat(temp2);

%% Print to file
i1 = regexp(fileName,'\'); % Find where site directory starts in path (PC).
i2 = regexp(fileName,'/'); % Find where site directory starts in path (MAC).
    if length(i1) > length(i2)
        fileNameResults = [cd '\' site '_charResults.csv'];
    else
        fileNameResults = [cd '/' site '_charResults.csv'];
    end 
fid=fopen(fileNameResults, 'wt');
for i=1:numrows
    fprintf(fid,'%s\n',output(i,:));
end
fclose(fid);
