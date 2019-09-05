function timeTable = read_timetable (sheetName, varargin)
%% Reads a time table from a spreadsheet file
% Usage: timeTable = read_timetable (sheetName, varargin)
% Explanation:
%       The built-in readtable function doesn't support timetables
%       So this is simply a wrapper function that converts a table
%           to a timetable after readtable
%
% Example(s):
%       load_examples;
%       write_timetable(myTimeTable1, 'myTimeTable1.csv');
%       myTimeTable1Copy = read_timetable('myTimeTable1.csv');
%
% Outputs:
%       timeTable   - time table read
%                   specified as a timetable
%
% Arguments:
%       sheetName   - spreadsheet file name for reading
%                   must be a string scalar or a character vector
%       varargin    - Any other parameter-value pair for writetable()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       TODO

% File History:
% 2019-08-27 Created by Adam Lu
% 

%% Hard-coded parameters

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'sheetName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, sheetName, varargin{:});

% Keep unmatched arguments for the writetable() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Read the table
tableRead = readtable(sheetName, otherArguments{:});

% Get all column names
columnNames = tableRead.Properties.VariableNames;

% Convert time vector to a duration or datetime array
if any(ismatch(columnNames, 'Time'))
    % Extract time vector
    timeVector = tableRead.Time;

    % Remove the time variable
    removevars(tableRead, 'Time');

    % If the time vector is text, convert it to duration or datetime
    tableRead.Time = text2time(timeVector);
else
    disp("A time table must have a 'Time' column!");
    timeTable = [];
    return;
end

% Convert to a timetable
timeTable = table2timetable(tableRead);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function timeVector = text2time (textVector)
%% Converts text to duration or datetime
% TODO: Pull out as its own function

if isduration(textVector) || isdatetime(textVector)
    timeVector = textVector;
    return
end

% TODO: make more general
timeVectorBeforeMin = extractBefore(textVector, 'min');

% TODO
timeVectorNumbers = str2double(timeVectorBeforeMin);

% TODO: make more general
timeVector = minutes(timeVectorNumbers);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Rename 'rowTimes' -> 'Time' if exists
table = renamevars(table, 'rowTimes', 'Time');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%