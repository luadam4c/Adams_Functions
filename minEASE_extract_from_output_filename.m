function infoStruct = extract_from_minEASE_outputfilename (fileName)
%% Extracts information from a minEASE output filename
% Usage: infoStruct = extract_from_minEASE_outputfilename (fileName)
% Outputs:
%       TODO
% Arguments:
%       fileName    - file name for the minEASE output file
%                   Note: everything before the directionLabel
%                           or the sweepLabel with be considered
%                           part of the fileIdentifier
% Requires:
%       /home/Matlab/Adams_Functions/find_in_strings.m
%       /home/Matlab/Adams_Functions/sscanf_full.m
%
% Used by:
%       cd/minEASE_combine_events.m
%       cd/minEASE_load_output.m

% File History:
%   2018-07-26 Created by Adam Lu
%   2018-08-02 Fixed the combined output file condition
%   2018-08-02 Now uses find_in_strings.m in 'regexp' mode
%   2018-08-02 Now outputs everything in an infoStruct structure

%% Hard-coded parameters
% To be consistent with minEASE.m
individualOutputStr = '_output';
directionLabelRegexp = '[E|I]PSC';
sweepLabelRegexp = 'Swp\d*';
timeUnitsDefault = 'samples';

% To be consistent with minEASE_combine_events.m
timeUnitsRegexp = 'samples|ms';
combinedOutputStr = '_ALL_output';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'fileName', ...     % file name for the minEASE output file
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
%    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

%% Preparation
% Create an empty structure to be outputted if an error occurs
infoStruct = struct;

%% Do the job
% Get the file base
[~, fileBase, fileExt] = fileparts(fileName);

% Determine the file type
%   Note: Must check for combinedOutputStr before individualOutputStr
if any(strfind(fileBase, combinedOutputStr))
    fileType = 'combined';
elseif any(strfind(fileBase, individualOutputStr))
    fileType = 'individual';
else
    error('Cannot recognize this file name!');
end

% Remove the output string
if strcmp(fileType, 'combined')
    tempStr1 = strrep(fileBase, combinedOutputStr, '');
elseif strcmp(fileType, 'individual')
    tempStr1 = strrep(fileBase, individualOutputStr, '');
else
    error('Problem with code!');
end

% Split the resulting string with the underscore as the delimiter
partsToRead = strsplit(tempStr1, '_');

% Look for any matching directionLabel
[idxDirectionLabel, directionLabel] = ...
    find_in_strings(directionLabelRegexp, partsToRead, ...
                        'SearchMode', 'regexp', ...
                        'IgnoreCase', true);
if length(idxDirectionLabel) > 1
    % Too many direction labels are found, return error message
    fprintf('Too many direction labels are found for ''%s''!\n', fileName);
    return;
end

% Look for any matching sweepLabel
[idxSweepLabel, sweepLabel] = ...
    find_in_strings(sweepLabelRegexp, partsToRead, ...
                        'SearchMode', 'regexp', ...
                        'IgnoreCase', true);
if length(idxSweepLabel) > 1
    % Too many sweep labels are found, return error message
    fprintf('Too many sweep labels are found for ''%s''!\n', fileName);
    return;
end

% If there is a sweep label, get the sweep number.
%   Otherwise, use 'all'
if ~isempty(sweepLabel)
    sweepNumber = sscanf_full(sweepLabel, '%d');
else
    sweepNumber = 'all';
end

% Look for any matching timeUnits
idxTimeUnits = find_in_strings(timeUnitsRegexp, partsToRead, ...
                                    'SearchMode', 'regexp', ...
                                    'IgnoreCase', true);
if isempty(idxTimeUnits)
    % No time units are found, use default
    timeUnits = timeUnitsDefault;
elseif length(idxTimeUnits) > 1
    % Too many time units are found, return error message
    fprintf('Too many time units are found for ''%s''!\n', fileName);
    return;
else
    % Get the matching time unit
    timeUnits = partsToRead{idxTimeUnits};
end

% Get the file identifier:
%   Replace the direction label, sweep label and time units with spaces, 
%   then anything before the first space is considered the file identifier
if ~isempty(directionLabel)
    tempStr2 = strrep(tempStr1, ['_', directionLabel], ' ');
else
    tempStr2 = tempStr1;
end
if ~isempty(sweepLabel)
    tempStr3 = strrep(tempStr2, ['_', sweepLabel], ' ');
else
    tempStr3 = tempStr2;
end
if ~isempty(idxTimeUnits)
    tempStr4 = strrep(tempStr3, ['_', timeUnits], ' ');
else
    tempStr4 = tempStr3;
end
tempCell1 = strsplit(tempStr4, ' ');
fileIdentifier = tempCell1{1};

%% Output in structure
infoStruct.fileIdentifier = fileIdentifier;
infoStruct.directionLabel = directionLabel;
infoStruct.sweepLabel = sweepLabel;
infoStruct.sweepNumber = sweepNumber;
infoStruct.timeUnits = timeUnits;
infoStruct.fileExt = fileExt;
infoStruct.fileType = fileType;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fileExt = tempCell1{3};

if any(regexp(fileBase, [individualOutputSuffix, '$']))
    fileType = 1;
elseif any(regexp(fileBase, [combinedOutputSuffix, '$']))
    fileType = 2;
else
    error('Cannot recognize this file name!');
end

% Replace '_ALL_output_samples' or '_output' with a space
%   to get the full label for this output file
if fileType == 1
    fullLabel = strrep(fileBase, individualOutputSuffix, ' ');
elseif fileType == 2
    fullLabel = strrep(fileBase, combinedOutputSuffix, ' ');
else
    error('Problem with code!');
end

% Get the analysis label and possibly the sweep number
if fileType == 1
    % Replace '_Swp' by a space
    tempStr1 = strrep(fullLabel, '_Swp', ' ');

    % Use strsplit to get the analysis label and the sweep number
    tempCell1 = strsplit(tempStr1, ' ');
    analysisLabel = tempCell1{1};
    iSwp = str2double(tempCell1{2});
elseif fileType == 2
    % There is no sweep number for combined output files
    iSwp = 'all';

    % The analysis label is the same as the full label
    analysisLabel = fullLabel;
else
    error('Problem with code!');
end

% Apply strsplit again to get the direction label
tempCell2 = strsplit(analysisLabel, '_');
directionLabel = tempCell2{end};

% Use strrep to get the file identifier string
fileIdentifier = strrep(analysisLabel, ['_', directionLabel], '');

% Decide whether the file is an individual output file

%                   Note: the elements should be in the order:
%                       fileIdentifier, directionLabel, sweepLabel, 
%                           outputStr, timeUnits, fileExt
%                       where outputStr is either '_ALL_output' or '_output'

% Split the file base into many parts using the underscore as the delimiter
fileBaseParts = strsplit(fileBase, '_');

% Look for the direction label and its index
fileBaseParts

% Get the analysis label and possibly the sweep number
if fileType == 1
    % Replace '_Swp' by a space
    tempStr1 = strrep(fullLabel, '_Swp', ' ');

    % Use strsplit to get the analysis label and the sweep number
    tempCell1 = strsplit(tempStr1, ' ');
    analysisLabel = tempCell1{1};
    iSwp = str2double(tempCell1{2});
elseif fileType == 2
    % There is no sweep number for combined output files
    iSwp = 'all';

    % The analysis label is the same as the full label
    analysisLabel = fullLabel;
else
    error('Problem with code!');
end

% Apply strsplit again to get the direction label
tempCell2 = strsplit(analysisLabel, '_');
directionLabel = tempCell2{end};

% Use strrep to get the file identifier string
fileIdentifier = strrep(analysisLabel, ['_', directionLabel], '');

% Remove the file identifier
tempStr4 = strrep(tempStr3, fileIdentifier, '');

% Split this string again with the underscore as the delimiter
partsToRead2 = strsplit(tempStr4, '_');


%}
