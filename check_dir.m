function check_dir (directories, varargin)
%% Checks if needed directory(ies) exist and create them if not
% Usage: check_dir (directories, varargin)
% Arguments:
%       directories - directory(ies) to check
%                   must be a cell array of character arrays
%                       or a scalar text
%       varargin    - 'Verbose': whether to write to standard output
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Used by:
%       cd/plot_traces_abf.m
%       /media/adamX/m3ha/data_dclamp/take4/find_initial_slopes.m
%
% File History:
% 2018-06-21 Modified from check_subdir.m
% 2018-09-18 Added input parser and verbose
% TODO: Use print_or_show_message

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
        if exist(directories{k}, 'dir') ~= 7
            mkdir(directories{k});
            if verbose
                fprintf('New directory is made: %s\n\n', ...
                    directories{k});
            end
        end
    end
else
    if exist(directories, 'dir') ~= 7
        mkdir(directories);
        if verbose
            fprintf('New directory is made: %s\n\n', ...
                    directories);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
