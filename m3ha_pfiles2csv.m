function sheetPaths = m3ha_pfiles2csv (varargin)
%% Converts the old m3ha .p files to spreadsheet files
% Usage: sheetPaths = m3ha_pfiles2csv (varargin)
% Explanation:
%       TODO
% Example(s):
%       m3ha_pfiles2csv;
%       m3ha_pfiles2csv('Recursive', true);
%
% Outputs:
%       sheetPaths  - full paths to the new spreasheets
%                   specified as a cell array of character vectors
% Arguments:
%       varargin    - 'Directory': the name of the directory containing 
%                                   the .p files, e.g. 'bestparams_E092910'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileNames': names of .p files to detect
%                   must be a characeter vector, a string array 
%                       or a cell array of character arrays
%                   default == detect from pwd
%                   - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'Recursive': whether to search recursively
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the save_params() function
%
% Requires:
%       cd/all_files.m
%       cd/create_error_for_nargin.m
%       cd/read_params.m
%       cd/save_params.m
%       cd/struct2arglist.m
%
% Used by:

% File History:
% 2019-08-14 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
directoryDefault = pwd;         % look for .p files in 
                                %   the present working directory by default
fileNamesDefault = {};          % detect from pwd by default
verboseDefault = true;          % print to standard output by default
recursiveDefault = false;       % don't search recursively by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Directory', directoryDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'FileNames', fileNamesDefault, ...
    @(x) isempty(x) || ischar(x) || iscellstr(x) || isstring(x));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Recursive', recursiveDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
directory = iP.Results.Directory;
fileNames = iP.Results.FileNames;
verbose = iP.Results.Verbose;
recursive = iP.Results.Recursive;

% Keep unmatched arguments for the save_params() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decide on the files to use
if isempty(fileNames)
    % Find all .p files in the directory
    [~, fileNames] = all_files('Directory', directory, 'Extension', '.p', ...
                                'Verbose', verbose, 'Recursive', recursive);

    % Return usage message if no .p files found
    if isempty(fileNames)
        fprintf('No .p files found!\n');
        return
    end
    
    % Print message
    if verbose
        fprintf('Parsing all .p files in %s ...\n', directory);
    end
else
    % Print message
    if verbose
        fprintf('Parsing all .p files ...\n');
    end
end

% Construct spreadsheet paths
sheetPaths = replace(fileNames, '.p', '.csv');

%% Do the job
% Load the parameter tables
paramsTables = cellfun(@read_params, fileNames, 'UniformOutput', false);

% Write the parameter tables
cellfun(@(x, y) save_params(x, 'FileName', y, otherArguments{:}), ...
                paramsTables, sheetPaths, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
