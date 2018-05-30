function csvwrite_with_header (filename, M, varargin)
%% Write a comma-separated value file with given header
% Usage: csvwrite_with_header (filename, M, varargin)
% Explanation:
%   TODO
% Example(s):
%       TODO
% Arguments:    
%       TODO
%
% Used by:    
%       /home/Matlab/minEASE/minEASE.m
%       /home/Matlab/minEASE/combine_eventInfo.m
%
% File History:
% 2017-07-24 Moved from examine_gapfree_events.m
% 2018-04-13 BT - Implemented Input Parser and added headers
% 2018-05-29 AL - Improved input parser
%

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
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty', 'scalartext'}));
addRequired(iP, 'M', ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'RowHeader', {}, ...
    @(x) assert(iscell(x) && (all(cellfun(@ischar, x)) ...
                || all(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));
addParameter(iP, 'ColumnHeader', {}, ...
    @(x) assert(iscell(x) && (all(cellfun(@ischar, x)) ...
                || all(cellfun(@isstring, x))), ...
                ['Second input must be a cell array ', ...
                'of strings or character arrays!']));

% Read from the Input Parser
parse(iP, filename, M, varargin{:});
filename = iP.Results.filename;
M = iP.Results.M;
rowHeader = iP.Results.RowHeader;
columnHeader = iP.Results.ColumnHeader;

%% Do the job
% Open file for writing
fid = fopen(filename, 'w');

% Print header
if ~isempty(columnHeader)
    if ~isempty(rowHeader)
        fprintf(fid, '%s, ', '');
    end
    for j = 1:numel(columnHeader)
        fprintf(fid, '%s, ', columnHeader{j});
    end
    fprintf(fid, '\n');
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

    % Print matrix contents
    for i = 1:size(M, 1)
        if ~isempty(rowHeader)
            fprintf(fid, '%s, ', rowHeader{i});
        end
        for j = 1:size(M, 2)
            fprintf(fid, '%g, ', M(i, j));
        end
        fprintf(fid, '\n');
    end
end

% Close file
fclose(fid);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}