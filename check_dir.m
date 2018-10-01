function check_dir (directories, varargin)
%% Checks if needed directory(ies) exist and creates them if not
% Usage: check_dir (directories, varargin)
% Arguments:
%       directories - directory(ies) to check
%                   must be a cell array of character arrays
%                       or a scalar text
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
% Requires:
%       cd/construct_fullfilename.m
%
% Used by:
%       cd/check_subdir.m
%       cd/compute_and_plot_evoked_LFP.m
%       cd/create_input_file.m
%       cd/parse_all_abfs.m
%       cd/plot_all_abfs.m
%       cd/plot_traces_abf.m
%       /media/adamX/m3ha/data_dclamp/take4/find_initial_slopes.m
%
% File History:
% 2018-06-21 Modified from check_subdir.m
% 2018-09-18 Added input parser and verbose
% TODO: Use print_or_show_message

%% Hard-coded parameters
mtitle = 'New Directory Made';              % message title

% TODO:Make this an optional parameter
messageMode= 'show';

%% Default values for optional arguments
verboseDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'directories', ...
    @(x) iscellstr(directories) || ischar(directories) || ...
        isstring(directories));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, directories, varargin{:});
verbose = iP.Results.Verbose;

%% Check directory(ies)
if iscell(directories)
    for k = 1:numel(directories)        
        % Construct the full path to the directory
        directory = construct_fullfilename(directories{k});

        % Check the directory and create it if it doesn't already exist
        check_dir_helper(directory, mtitle, messageMode, verbose);
    end
else
    % Construct the full path to the directory
    directory = construct_fullfilename(directories);
    
    % Check the directory and create it if it doesn't already exist
    check_dir_helper(directory, mtitle, messageMode, verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function check_dir_helper(directory, mtitle, messageMode, verbose)

% Check if the directory exists
if exist(directory, 'dir') ~= 7
    % Create the directory
    mkdir(directory);

    % Show message and print to standard output
    msg = sprintf('New directory is made: %s\n\n', directory);
    print_or_show_message(msg, 'MTitle', mtitle, ...
                                'MessageMode', messageMode, ...
                                'Verbose', verbose);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{ 
OLD CODE:

if verbose
    fprintf('New directory is made: %s\n\n', directory);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
