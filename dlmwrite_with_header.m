function dlmwrite_with_header (filename, M, varargin)
%% Write a comma-separated value file with given header
% Usage: dlmwrite_with_header (filename, M, varargin)
% Explanation:
%   TODO
%
% Example(s):
%       dlmwrite_with_header('magic3.txt', magic(3));
%       dlmwrite_with_header('magic3_h.txt', magic(3), 'ColumnHeader', {'A', 'B', 'C'});
%       dlmwrite_with_header('magic4.csv', magic(4));
%       dlmwrite_with_header('magic3_tab.txt', magic(3), 'Delimiter', '\t');
%       dlmwrite_with_header('magic3_p2.txt', magic(3), 'Precision', 2);
%
% Arguments:    
%       filename    - file name
%                   must be a string scalar or a character vector 
%       varargin    - 'Delimiter': delimiter used to separate entries
%                   must be a string scalar or a character vector
%                   default == ','
%                   - 'Precision': precision for the number
%                   must be a number or a string scalar or a character vector
%                   default == '%g'
%                   - 'RowHeader': row header
%                   must be a string vector or cell array of character vectors
%                   default == none
%                   - 'ColumnHeader': row header
%                   must be a string vector or cell array of character vectors
%                   default == none
%                   - 'AppendToFile': whether to append to the file
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Used by:    
%       cd/m3ha_autocorrelogram.m
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/minEASE/combine_eventInfo.m
%       /home/Matlab/minEASE/filter_minEASE_output.m
%       /home/Matlab/plethRO1/spike2Loader.m
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%
% File History:
% 2017-07-24 Moved from examine_gapfree_events
% 2019-09-02 Renamed csvwrite_with_header -> dlmwrite_with_header.m
% 2019-09-02 Added 'Delimiter' and 'Precision' as optional arguments
% 2019-09-02 Added 'AppendToFile'
% TODO: Rename as dlmwrite_custom
% TODO: Use dlmwrite.m with the '-append' option and '%g' precision mode by default
% TODO: Use print_cellstr.m
%

%% Default values for optional arguments
delimiterDefault = ',';
precisionDefault = '';      % set later
rowHeaderDefault = {};
columnHeaderDefault = {};
appendToFileDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end


% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'filename', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addRequired(iP, 'M', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Delimiter', delimiterDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Precision', precisionDefault);
addParameter(iP, 'RowHeader', rowHeaderDefault, ...
    @(x) assert(iscell(x) || isstring(x) || iscellstr(x), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));
addParameter(iP, 'ColumnHeader', columnHeaderDefault, ...
    @(x) assert(iscell(x) || isstring(x) || iscellstr(x), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));
addParameter(iP, 'AppendToFile', appendToFileDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, filename, M, varargin{:});
delimiter = iP.Results.Delimiter;               % Examples: ',' '\t'
precision = iP.Results.Precision;
rowHeader = iP.Results.RowHeader;
columnHeader = iP.Results.ColumnHeader;
appendToFile = iP.Results.AppendToFile;             % appendToFile is true/false

%% Preparation
% Decide on the precision
precision = decide_on_precision(precision);

%% Do the job
% Open file for writing
if appendToFile
    fid = fopen(filename, 'a');
else
    fid = fopen(filename, 'w');
end

% Print header
if ~isempty(columnHeader)
    % Print a column for the row header
    if ~isempty(rowHeader)
        fprintf(fid, ['%s', delimiter], 'RowName');
    end

    % Count the number of columns
    nColumns = numel(columnHeader);

    % Print all headers except the last, followed by the delimiter
    if nColumns >= 2
        for iColumn = 1:(nColumns-1)
            fprintf(fid, ['%s', delimiter], columnHeader{iColumn});
        end
    end

    % Print the last header, followed by the newline character
    fprintf(fid, '%s\n', columnHeader{nColumns});

end

% Print matrix contents
if ~isempty(M)
    % Make sure header is the correct length
    if numel(columnHeader) ~= size(M, 2) && ~isempty(columnHeader)
        error('Column header does not have the correct length!');
    end
    if numel(rowHeader) ~= size(M, 1) && ~isempty(rowHeader)
        error('Row header does not have the correct length!');
    end

    % Count the number of rows
    nRows = size(M, 1);

    % Count the number of columns
    nColumns = size(M, 2);

    % Print matrix contents
    for iRow = 1:nRows
        % Print the row header
        if ~isempty(rowHeader)
            fprintf(fid, '%s, ', rowHeader{iRow});
        end

        % Print all columns of this row except the last, 
        %       followed by the delimiter
        if nColumns >= 2
            for iColumn = 1:(nColumns - 1)
                fprintf(fid, [precision, delimiter], M(iRow, iColumn));
            end
        end

        % Print the last column, followed by the newline character
        fprintf(fid, [precision, '\n'], M(iRow, nColumns));
    end
end

% Close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function precision = decide_on_precision(precision)

if isempty(precision)
    precision = '%g';
elseif isnumeric(precision)
    precision = ['%.', num2str(precision), 'g'];
else
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
