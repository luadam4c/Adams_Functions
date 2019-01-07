function vvecsNew = m3ha_correct_unbalanced_bridge(fileNames, vvecsOld, ivecsOld, initialSlopesFile)
%% Fix current pulse response traces that may have out-of-balance bridges
% Usage: vvecsNew = m3ha_correct_unbalanced_bridge(fileNames, vvecsOld, ivecsOld, initialSlopesFile)
% Arguments:    
%       fileNames   - file names for the vectors
%                   must be a string scalar or a character vector
%       vvecsOld    - Voltage vectors to correct
%                   must be a numeric vector
%       ivecsOld    - Corresponding current vectors
%                   must be a numeric vector
%       initialSlopesFile - path to initial slopes file
%                   must be a string scalar or a character vector
%
% Requires:
%       cd/correct_unbalanced_bridge.m
%       cd/find_ind_str_in_cell.m
%
% Used by:    
%       /media/adamX/m3ha/data_dclamp/dclampPassiveFitter.m
%       /media/adamX/m3ha/optimizer4gabab/import_rawtraces.m TODO

% File History:
% 2018-10-15 Adapted from /media/adamX/m3ha/optimizer4gabab/import_rawtraces.m
% TODO: Input parser
% TODO: Make vvecsOld, ivecsOld, initialSlopesFile, diaryFile optional arguments

initialSlopesFileDefault = ...
    fullfile('~', 'm3ha', 'data_dclamp', 'take4', ...
            'initial_slopes_nSamplesForPlot_2_threeStdMainComponent.mat');
fid = 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Count the number of sweeps
nSwps = size(vvecsOld, 2);

% Load the initial slopes output file
initialSlopes = matfile(initialSlopesFile);

% Extract necessary variables
initialSlopeFilenames = initialSlopes.filenamesSorted;
initialSlopeThreshold1IndexBalanced = initialSlopes.iThreshold1Balanced;
initialSlopeThreshold2IndexBalanced = initialSlopes.iThreshold2Balanced;

% Find the index of file in all files sorted by initial slope
%   in descending order
ftemp = @(x) find_ind_str_in_cell(x, initialSlopeFilenames);

% Determine whether the initial slopes exceed threshold
%   Note: These may have out-of-balance bridges
isOutOfBalance = ...
    cellfun(@(x) ftemp(x) < initialSlopeThreshold2IndexBalanced || ...
                ftemp(x) > initialSlopeThreshold1IndexBalanced, ...
                fileNames);

% Print out an appropriate message
if any(isOutOfBalance)
    fprintf(fid, ['The following current pulse responses will be ', ...
                    'corrected due to out-of-balance bridges:\n']);
    print_cellstr(fileNames(isOutOfBalance), ...
                  'FileID', fid, 'OmitBraces', true, ...
                  'OmitQuotes', true, 'Delimiter', '\n', 'Prefix', '\t');
else
    fprintf(fid, ['There are no current pulse response traces ', ...
                    'with out-of-balance bridges!\n']);
end

% Correct for traces that may have out-of-balance bridges
vvecsNew = zeros(size(vvecsOld));
parfor iSwp = 1:nSwps
    % Get the old voltage trace
    vvecOld = vvecsOld(:, iSwp);
    if isOutOfBalance(iSwp)
        % Get the old current trace
        ivecOld = ivecsOld(:, iSwp);

        % Correct any unbalanced bridge in the voltage trace
        vvecNew = correct_unbalanced_bridge(vvecOld, ivecOld);

        % Store the new trace
        vvecsNew(:, iSwp) = vvecNew;
    else
        % Store the old trace
        vvecsNew(:, iSwp) = vvecOld;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

load(initialSlopesFile, 'filenamesSorted', ...
        'iThreshold1Balanced', 'iThreshold2Balanced'); 

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%