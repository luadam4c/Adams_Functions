function swdTable = read_swd_sheet (filePath)
%% Read in an SWD table from a spreadsheet file
% Usage: swdTable = read_swd_sheet (filePath)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       swdTable    - SWD table loaded
%                   specified as a 2D table
% Arguments:
%       filePath    - SWD spreadsheet file to read
%                   must be a string scalar or a character vector
%
% Used by:
%       cd/load_swd_sheets.m

% File History:
% 2018-11-27 Created by Adam Lu
% 

%% Hard-coded parameters
% Must be consistent with parse_assyst_swd.m
dateTimePattern = 'M/d/yyyy HH:mm:ss.SSS';

%% Default values for optional arguments

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
addRequired(iP, 'filePath', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, filePath);

%% Preparation

%% Do the job
% Read in the SWD table
swdTable = readtable(filePath);

% Get all variables
varNames = swdTable.Properties.VariableNames;

% Process further if there is data
if height(swdTable) ~= 0
    % Convert strings to datetime
    if any(strcmpi(varNames, 'startTimeOrig'))
        swdTable.startTimeOrig = ...
            datetime(swdTable.startTimeOrig, 'InputFormat', dateTimePattern);    
    end
    if any(strcmpi(varNames, 'endTimeOrig'))
        swdTable.endTimeOrig = ...
            datetime(swdTable.endTimeOrig, 'InputFormat', dateTimePattern);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

swdTable.startTimeOrig = ...
    cellfun(@(x) textscan(x, dateTimeFormatSpec), swdTable.startTimeOrig);

startTimeOrigDateTime = ...
    cellfun(@(x) datetime(x, 'InputFormat', dateTimePattern), ...
            swdTable.startTimeOrig);

startTimeOrigDateTime = ...
    datetime(swdTable.startTimeOrig, 'InputFormat', dateTimePattern);

swdTable = addvars(swdTable, startTimeOrigDateTime, 'before', 'startTimeOrig');

swdTable = removevars(swdTable, 'startTimeOrig');

swdTable.Properties.VariableNames{'startTimeOrigDateTime'} = 'startTimeOrig';

if any(strcmpi(varNames, 'durationOrig'))
    swdTable.durationOrig = ...
        duration(swdTable.durationOrig, 'InputFormat', durationPattern);    
end

durationPattern = 'hh:mm:ss.SSS';

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%