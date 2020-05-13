function table = write_timetable (timeTable, varargin)
%% Writes a time table to a spreadsheet file
% Usage: table = write_timetable (timeTable, sheetName (opt), varargin)
% Explanation:
%       The built-in writetable function doesn't support timetables
%       So this is simply a wrapper function that converts a timetable
%           to a table, then save
%
% Example(s):
%       load_examples;
%       write_timetable(myTimeTable1, 'myTimeTable1.csv');
%
% Outputs:
%       table       - converted table
%                   specified as a table
%
% Arguments:
%       timeTable   - time table to save
%                   must be a time table
%       sheetName   - (opt) spreadsheet file name for saving
%                   must be a string scalar or a character vector
%                   default == specified by writetable()
%       varargin    - Any other parameter-value pair for writetable()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_population_average.m
%       cd/plot_measures.m

% File History:
% 2019-08-27 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
sheetNameDefault = '';

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
addRequired(iP, 'timeTable', ...
    @(x) validateattributes(x, {'timetable'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'sheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, timeTable, varargin{:});
sheetName = iP.Results.sheetName;

% Keep unmatched arguments for the writetable() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Convert to a table
table = timetable2table(timeTable);

% Save the table
if isempty(sheetName)
    writetable(table, otherArguments{:});
else
    writetable(table, sheetName, otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Rename 'rowTimes' -> 'Time' if exists
table = renamevars_custom(table, 'rowTimes', 'Time');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
