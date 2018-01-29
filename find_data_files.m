function [dataType, allDataFiles, nDataFiles, message] = find_data_files (dataTypeUser, dataDirectory, possibleDataTypes, varargin)
%% Looks for data files in a dataDirectory according to either dataTypeUser or going through a list of possibleDataTypes
% Usage: [dataType, allDataFiles, nDataFiles, message] = find_data_files (dataTypeUser, dataDirectory, possibleDataTypes, varargin)
%
% Arguments:
%       TODO    
%       varargin    - 'FileIdentifier': data file identifier
%                   must be a string scalar or a character vector
%                   default == ''
% Used by:
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/Adams_Functions/combine_sweeps.m
%
% File History:
%   2018-01-29 Moved from /home/Matlab/minEASE/minEASE.m
%   2018-01-29 Added the case where dataTypeUser is not recognized 
%   2018-01-29 Added fileIdentifier as an optional argument
%   TODO: Make possibleDataTypes an optional argument? Default?
%

%% Default values for optional arguments
fileIdentifierDefault = '';     % no file identifier by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help find_data_files'' for usage']);
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'find_data_files';

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FileIdentifier', fileIdentifierDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));

% Read from the Input Parser
parse(iP, varargin{:});
fileIdentifier = iP.Results.FileIdentifier;

%% Perform job
% Extract the number of possible data types
nDataTypes = numel(possibleDataTypes);

% Find data files according to data type
switch dataTypeUser
case possibleDataTypes              % if user provided a possible data type
    dataType = dataTypeUser;
    allDataFiles = dir(fullfile(dataDirectory, ...
                    [fileIdentifier, '*.', dataTypeUser]));
    nDataFiles = length(allDataFiles);  % # of files in data subdirectory
    if nDataFiles == 0
        message = sprintf(['There are no %s files in this directory:\n', ...
                            '%s\n'], dataTypeUser, dataDirectory);
    end
case 'default'                      % if using default data types
    % Iterate through possibleDataTypes to look for possible data type
    for iDataType = 1:nDataTypes
        tempType = possibleDataTypes{iDataType};
        allDataFiles = dir(fullfile(dataDirectory, ...
                        [fileIdentifier, '*.', tempType]));
        nDataFiles = length(allDataFiles);% # of files in data subdirectory
        if nDataFiles > 0
            dataType = tempType;
            message = sprintf(['The .%s files in this directory ', ...
                                'will be used as data:\n', '%s\n\n'], ...
                                dataType, dataDirectory);
            break;
        end
    end
    if nDataFiles == 0
        dataType = '';
        message = sprintf(['There are no acceptable data files', ...
                            ' in this directory:\n', ...
                            '%s\n'], dataDirectory);
    end
otherwise
    dataType = '';
    allDataFiles = [];
    nDataFiles = 0;

    % Start message
    message = sprintf(['The data type %s is not recognized!'
                        'The possible data types are:\n'], dataTypeUser);
                        
    % Print possible data types in a line
    for iDataType = 1:nDataTypes
        message = [message, sprintf('%s', possibleDataTypes{iDataType})];
        if iDataType < nDataTypes
            message = [message, ', '];
        else
            message = [message, '\n'];
        end
    end
    % TODO: Make more general into print_cell.m
    %       function string = print_cell(cellArray, varargin)
    %       %% Prints a cell array of strings into a single line with each entry separated by a delimiter (default ',')
    %       % Usage: string = print_cell(cellArray, varargin)
    %       %   
    %       % Arguments:
    %       %       varargin    - 'Delimiter' - delimiter used to separate entries
    %       %                   default == ','
    %       %                   - 'OmitNewline' - whether to omit the newline character
    %       %                   default == false
    
    % End message
    message = [message, 'You could also say ''default''.\n'];


end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
