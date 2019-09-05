function outTables = create_row_labels (inTables, varargin)
%% Creates row labels for table(s)
% Usage: outTables = create_row_labels (inTables, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       outTables   - output tables
%                   specified as a TODO
% Arguments:
%       inTables    - input tables
%                   must be a table or a cell array of tables
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair 
%                       for create_labels_from_numbers.m
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/combine_variables_across_tables.m

% File History:
% 2019-03-17 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'inTables', ...
    @(x) validateattributes(x, {'table', 'cell'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, inTables, varargin{:});
param1 = iP.Results.param1;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
if iscell(inTables)
    outTables = cellfun(@(x) create_row_labels_helper(x, otherArguments), ...
                        inTables, 'UniformOutput', false);
else
    outTables = create_row_labels_helper(inTables, otherArguments);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tble = create_row_labels_helper(tble, otherArguments)
%% Creates row labels from row numbers

% Count the number of rows
nRows = height(tble);

% Create row labels
tble.Properties.RowNames = ...
    create_labels_from_numbers(1:nRows, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%