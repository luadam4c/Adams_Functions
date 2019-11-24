function [statuses, messages, messageIds] = ...
                copy_into (sources, destFolder, varargin)
%% Copies source directories and files into a destination folder
% Usage: [statuses, messages, messageIds] = ...
%               copy_into (sources, destFolder, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
%
% Arguments:
%       sources     - source directories/files
%                   must be a character vector or a string vector
%                       or a cell array of character vectors
%       destFolder  - destination folder
%                   must be a string scalar or a character vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%
% Used by:
%       /media/adamX/m3ha/optimizer4compgabab/singleneuronfitting63.m

% File History:
% 2019-11-24 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'sources', ...
    @(x) assert(isempty(x) || ischar(x) || iscellstr(x) || isstring(x), ...
        ['sources must be empty or a character array or a string array ', ...
            'or cell array of character arrays!']));
addRequired(iP, 'destFolder', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, sources, destFolder, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
if iscell(sources)
    [statuses, messages, messageIds] = ...
        cellfun(@(x) copy_into_helper(x, destFolder), ...
                sources, 'UniformOutput', false);
elseif isstring(sources)
    [statuses, messages, messageIds] = ...
        arrayfun(@(x) copy_into_helper(x, destFolder), ...
                sources, 'UniformOutput', false);
else
    [statuses, messages, messageIds] = copy_into_helper(sources, destFolder);
end

%% Deal with outputs
% Force statuses as a numeric column
if iscell(statuses)
    statuses = vertcat(statuses{:});
end

% Force messages as column cell arrays
[messages, messageIds] = argfun(@force_column_cell, messages, messageIds);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [status, message, messageId] = copy_into_helper (source, destFolder)

% Decide on the destination acceptable by copyfile
if isfolder(source)
    % Grab the directory base
    dirBase = extract_fileparts(source, 'dirbase');

    % Construct a subdirectory with the same directory base inside destFolder
    destination = fullfile(destFolder, dirBase);
else
    destination = destFolder;
end

% Copy the source to destination
[status, message, messageId] = copyfile(source, destination);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%