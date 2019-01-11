function string = print_cellstr (cellStr, varargin)
%% Prints and returns a string for the contents stored in a cell array
% Usage: string = print_cellstr (cellStr, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Side Effects:
%       Prints to standard output or a file with given FileID
% Outputs:
%       string      - a character array that includes elements of a cell array
% Arguments:
%       cellStr   - a cell array of character arrays to be printed
%                       or a character array
%                   must be a cell array of character arrays or a character array
%       varargin    - 'Delimiter': used to delimit separate entries
%                   must be a character array
%                   default == ''
%                   - 'ColDelimiter': used to delimit separate columns
%                   must be a character array
%                   default == ', '
%                   - 'RowDelimiter': used to delimit separate rows
%                   must be a character array
%                   default == '; '
%                   - 'Prefix': prefix to prepend to each element
%                   must be a character array
%                   default == ''
%                   - 'Suffix': suffix to append to each element
%                   must be a character array
%                   default == ''
%                   - 'OmitQuotes': whether to omit quotes
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OmitBraces': whether to omit braces
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OmitNewline': whether to omit the newline character
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'ToPrint': whether to actually print the string
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'FileID': file ID returned by fopen()
%                   must be an integer
%                   default == 1 (standard output)
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_parse_mat.m
%       cd/m3ha_select_cells.m
%       cd/m3ha_xolotl_create_neuron.m
%       cd/print_or_show_message.m
%       cd/print_and_show_message.m
%       cd/print_structure.m
%
% File History:
% 2018-01-31 Modified from code in find_data_files.m by Mira
% 2018-02-01 AL - Modified descriptions
% 2018-02-02 AL - This function will now be used by print_or_show_message.m
% 2018-02-05 MD - Shortened code and modified optional arguments
% 2018-02-05 MD - Allowed first argument to be just a character array as well
% 2018-06-21 AL - Changed default delimiter from ',' to ', '
% 2018-06-21 AL - Changed iP.KeepUnmatched to false
% 2018-06-21 AL - Added 'ToPrint', 'FileID'
% 2018-06-21 AL - Now also prints to standard output by default
% TODO: Consider 3-D cell arrays of strings

%% Hard-coded parameters
newLineChar = '\n';

%% Default values for optional arguments
delimiterDefault = '';
rowDelimiterDefault = '';               % default: separate rows with '; '
colDelimiterDefault = '';               % default: separate columns with ', '
prefixDefault = '';                     % default: prepend nothing
suffixDefault = '';                     % default: append nothing
omitQuotesDefault = false;              % whether to omit quotes by default
omitBracesDefault = false;              % whether to omit braces by default
omitNewlineDefault = false;             % whether to omit a new line by default
toPrintDefault = true;                  % default: print the string
fileIdDefault = 1;                      % default: print to standard output

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
addRequired(iP, 'cellStr', ...
    @(x) iscellstr(x) || ischar(x) || isstring(x));     

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RowDelimiter', rowDelimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'ColDelimiter', colDelimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OmitQuotes', omitQuotesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OmitBraces', omitBracesDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OmitNewline', omitNewlineDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'ToPrint', toPrintDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'FileId', fileIdDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'integer'}));

% Read from the Input Parser
parse(iP, cellStr, varargin{:});
delimiter = iP.Results.Delimiter;               % Examples: ',' '/'
rowDelimiter = iP.Results.RowDelimiter;
colDelimiter = iP.Results.ColDelimiter;
prefix = iP.Results.Prefix;
suffix = iP.Results.Suffix;
omitQuotes = iP.Results.OmitQuotes;             % omitQuotes is true/false
omitBraces = iP.Results.OmitBraces;             % omitBraces is true/false
omitNewline = iP.Results.OmitNewline;           % omitNewline is true/false
toPrint = iP.Results.ToPrint;
fileId = iP.Results.FileId;

% Set dependent argument defaults
if ~isempty(delimiter)
    if ~isempty(rowDelimiter)
        fprintf('Warning: Delimiter ignored as RowDelimiter is defined!\n');
    else
        rowDelimiter = delimiter;
    end
    if ~isempty(colDelimiter)
        fprintf('Warning: Delimiter ignored as ColDelimiter is defined!\n');
    else
        colDelimiter = delimiter;
    end
else
    if isempty(rowDelimiter)
        rowDelimiter = '; ';
    end
    if isempty(colDelimiter)
        colDelimiter = ', ';
    end
end

%% Perform task
% Prepend prefices and append suffices
cellStr = strcat(prefix, cellStr, suffix);

% Print elements of cellStr in a single string delimited by delimiter
if ischar(cellStr)                     % a character vector
    if ~omitQuotes
        string = ['''', cellStr, ''''];
    else
        string = cellStr;
    end
elseif isstring(cellStr)                % a string array
    % TODO: deal with string arrays with numel > 1
    if ~omitQuotes
        string = ['"', cellStr, '"'];
    else
        string = cellStr;
    end
else                                    % a cell array of character vectors
    % Add single quotes to each string if requested
    if ~omitQuotes
        cellStr = cellfun(@(x) ['''', x, ''''], cellStr, ...
                'UniformOutput', false);
    end

    % Count the number of rows
    nRows = size(cellStr, 1);

    % Count the number of columns
    nColumns = size(cellStr, 2);

    % Initialize the final string
    string = '';

    % Append each element of cellStr and necessary delimiters
    for iRow = 1:nRows                  % for each row
        for iCol = 1:nColumns           % for each column
            % Append the element
            string = [string, cellStr{iRow, iCol}];  

            % If we are not at the last column, add a column delimiter
            if iCol < nColumns
                string = [string, colDelimiter];
            end

            % If we are at the last column but not at the last row,
            %   add a row delimiter
            if iCol == nColumns && iRow < nRows
                string = [string, rowDelimiter];
            end
        end
    end
end

% Close the string out with braces if requested
if ~omitBraces
    string = ['{', string, '}'];
end

% Append a newline character at the end if requested
if ~omitNewline                         % if not omitNewline
    string = [string, newLineChar];
end

% Print the single string to standard output or a file
if toPrint
    fprintf(fileId, string);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

%       function string = print_cellstr(cellStr, varargin)
%       %% Prints a cell array of strings into a single line with each entry separated by a delimiter (default ',')
%       % Usage: string = print_cellstr(cellStr, varargin)
%       %   
%       % Arguments:
%       %       varargin    - 'Delimiter' - delimiter used to separate entries
%       %                   default == ','
%       %                   - 'OmitNewline' - whether to omit the newline character
%       %                   default == false

%string = cellStr{1};%taking the first element of the cellStr 
%string = [string, ', '];
%for iElement = 1:nElements          % for all elements except 1 in cellStr

string = '';
for iElement = 1:nElements              % for all elements in cellStr
    string = [string, cellStr{iElement}];   % concatenate strings  
    if iElement < nElements             % if we are not at the last element
        string = [string, ', '];        % put a comma after
    else                                % if we are at the last element
        string = [string, '\n'];        % put a newline character after
    end
end

string = strjoin(cellStr, ', ');    % concatenate strings  
string = [string, '\n'];            % put a newline character after

delimiterWithSpace = [delimiter, ' '];

if omitNewline == false

iP.KeepUnmatched = true;                        % allow extraneous options
% Display warning message if some inputs are unmatched
if ~isempty(fieldnames(iP.Unmatched))
    fprintf('WARNING: The following name-value pairs could not be parsed: \n');
    disp(iP.Unmatched);
end

@(x) iscellstr(x) || ischar(x));     

string = strjoin(cellStr, delimiter);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%