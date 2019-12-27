function lineStrs = read_lines_from_file (filePath, varargin)
%% Reads line(s) from a file
% Usage: lineStrs = read_lines_from_file (filePath, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       path = '/home/Matlab/Adams_Functions/read_lines_from_file.m'
%       lineStrs = read_lines_from_file(path, 'LineNumber', 2)
%       lineStrs = read_lines_from_file(path, 'Keyword', 'unix')
%       lineStrs = read_lines_from_file(path, 'Keyword', 'unix', 'MaxNum', 1)
%
% Outputs:
%       lineStrs    - line(s) read
%                   specified as a character vector 
%                       or a cell array of character vectors
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
%                   default == Inf
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/count_strings.m
%       cd/force_column_cell.m
%       cd/ispositiveintegerscalar.m
%       cd/struct2arglist.m
%       /TODO:dir/TODO:file
%
% Used by:
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_atf_swd.m
%       cd/parse_iox.m

% File History:
% 2019-09-08 Created by Adam Lu
% 2019-09-10 Added 'MaxNum' as an optional argument
% TODO FOR UNDERGRAD: complete 'IncludeNewLine' as an optional argument
% TODO: Expand to multiple lines

%% Hard-coded parameters

%% Default values for optional arguments
lineNumberDefault = [];
keywordDefault = '';
includeNewLineDefault = [];             % default TODO: Description of includeNewLine
maxNumDefault = Inf;

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
% Initialize output
lineStrs = {};

% Decide on the function to use based on includeNewLine
% TODO: Hint: fgetl vs. fgets

% Return if file does not exist
if ~isfile(filePath)
    fprintf('The file %s doesn''t exist!\n', filePath);
    return
end

% Update maxNum if lineNumber provided
if ~isempty(lineNumber)
    if ~isinf(maxNum) && maxNum <= 1
        fprintf('One can only read at most 1 line with a specific line number!\n');
        return;
    else
        maxNum = 1;
    end
end

%% Do the job
% if isunix
if false
    if ~isempty(lineNumber)
        [~, unixOut] = unix(sprintf('sed -n ''%dp'' %s', lineNumber, filePath));
    else
        % Use grep
        if isinf(maxNum)
            [~, unixOut] = unix(sprintf('grep "%s" %s', keyword, filePath));
        else
            [~, unixOut] = unix(sprintf('grep -m %d "%s" %s', ...
                                            maxNum, keyword, filePath));
        end
    end

    % Split output by endline character
    lineStrs = strsplit(unixOut, '\n');

    % Remove the last string if it's empty
    if isempty(lineStrs{end})
        lineStrs = lineStrs(1:end-1);
    end
else
    lineStrs = read_lines_from_file_alt(filePath, lineNumber, keyword, maxNum);
%     lineStrs = read_lines_from_files_slow(filePath, lineNumber, keyword, maxNum);
end

%% Output results
% Output as a character vector if there is only one line
if iscell(lineStrs) && numel(lineStrs) == 1
    lineStrs = lineStrs{1};
else
    lineStrs = force_column_cell(lineStrs);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lineStrs = read_lines_from_file_alt (filePath, lineNumber, ...
                                                keyword, maxNum)

% Read the file as a table
tableOfStrs = readtable(filePath, 'FileType', 'text', ...
                        'ReadVariableNames', false, 'Delimiter', '\n', ...
                        'HeaderLines', 0);

% Extract strings
if ~isempty(lineNumber)
    lineStrs = tableOfStrs{lineNumber, 'Var1'};
else
    lineStrs = tableOfStrs{:, 'Var1'};
end

% Restrict to specific keyword(s) if requested
if ~isempty(keyword)
    % Test whether each string contains a keyword
    containsKeyword = contains(lineStrs, keyword);

    % Restrict to those strings
    lineStrs = lineStrs(containsKeyword);
end

% Count how many lines are in lineStrs
nLines = count_strings(lineStrs);

% If there are more lines than maxNum, restrict to the first maxNum lines
if nLines > maxNum
    lineStrs = lineStrs(1:maxNum);    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function lineStrs = read_lines_from_files_slow (filePath, lineNumber, ...
                                                keyword, maxNum)
%% Read lines with MATLAB functions

% Open the file for reading
fid = fopen(filePath, 'r');

% Initialize a counter for stored lines
nStored = 0;

% Initialize a counter for read lines
nRead = 0;

% Read a line and increment counter
[lineStrThis, fid, nRead] = read_and_increment(fid, nRead);

% Read until maxNum or end of file is reached
while nStored <= maxNum && ischar(lineStrThis)
    % Ignore lines that are not needed
    if ~isempty(lineNumber)
        % Read lines until a line number is reached
        while nRead < lineNumber
            [lineStrThis, fid, nRead] = read_and_increment(fid, nRead);
        end
    elseif ~isempty(keyword)
        % Read lines until a line containing keyword is read
        while ~contains(lineStrThis, keyword)
            [lineStrThis, fid, nRead] = read_and_increment(fid, nRead);
        end
    end

    % Store the current line
    nStored = nStored + 1;
    lineStrs{nStored} = lineStrThis;
end

% Close the trace file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [lineStrThis, fid, nRead] = read_and_increment (fid, nRead)
%% Read line and increment counter

% Read new line
lineStrThis = fgetl(fid);

% Increment counter
nRead = nRead + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%