function lineStrs = read_lines_from_file (filePath, varargin)
%% Reads line(s) from a file
% Usage: lineStrs = read_lines_from_file (filePath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       lineStrs    - line(s) read
%                   specified as a character vector
%
% Arguments:
%       filePath    - path to file to read
%                   must be a string scalar or a character vector
%       varargin    - 'LineNumber': line number to read
%                   must be empty or a positive integer scalar
%                   default == []
%                   - 'Keyword': keyword for the line to read
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'IncludeNewLine': whether to include the newline character
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNum': maximum number of lines to read
%                   must be empty or a positive integer scalar
%                   default == read all lines of the file
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/ispositiveintegerscalar.m
%       cd/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/parse_atf_swd.m

% File History:
% 2019-09-08 Created by Adam Lu
% TODO FOR UNDERGRAD: complete 'IncludeNewLine' as an optional argument
% TODO: Expand to multiple lines
% TODO: Complete 'MaxNum' as an optional argument

%% Hard-coded parameters

%% Default values for optional arguments
lineNumberDefault = [];
keywordDefault = '';
includeNewLineDefault = [];             % default TODO: Description of includeNewLine
maxNumDefault = [];                     % TODO by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'filePath');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'LineNumber', lineNumberDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));
addParameter(iP, 'Keyword', keywordDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'IncludeNewLine', includeNewLineDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNum', maxNumDefault, ...
    @(x) assert(isempty(x) || ispositiveintegerscalar(x), ...
                'MaxNum must be either empty or a positive integer scalar!'));

% Read from the Input Parser
parse(iP, filePath, varargin{:});
lineNumber = iP.Results.LineNumber;
keyword = iP.Results.Keyword;
includeNewLine = iP.Results.IncludeNewLine;
maxNum = iP.Results.MaxNum;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

% Check relationships between arguments
% TODO

%% Preparation
% Decide on the function to use based on includeNewLine
% TODO: Hint: fgetl vs. fgets


%% Do the job
% Open the file for reading
fid = fopen(filePath, 'r');

% Read a line 
lineStrs = fgetl(fid);

% Ignore lines that are not needed
if ~isempty(lineNumber)
    % Read lines until a line number is reached
    lineCount = 1;
    while lineCount < lineNumber
        lineStrs = fgetl(fid);
        lineCount = lineCount + 1;
    end
elseif ~isempty(keyword)
    % Read lines until a line containing keyword is read
    while ~contains(lineStrs, keyword)
        lineStrs = fgetl(fid);
    end
end

% Close the trace file
fclose(fid);

%% Output results
% TODO

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%