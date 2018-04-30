function csvwrite_with_header (filename, M, varargin)
%% Write a comma-separated value file with given header
% Usage: csvwrite_with_header (filename, M, header)
% TODO
%
% Used by:    
%       /home/Matlab/minEASE/examine_gapfree_events.m
%
% File History:
% 2017-07-24 Moved from examine_gapfree_events.m
% 2018-04-13 BT - Implemented Input Parser and added headers
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

iP = inputParser;
addRequired(iP, 'filename', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty', 'scalartext'}));
addRequired(iP, 'M', ...
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addParameter(iP, 'RowHeader', {}, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
addParameter(iP, 'ColumnHeader', {}, ...
    @(x) validateattributes(x, {'cell'}, {'nonempty'}));
parse(iP, filename, M, varargin{:});
filename = iP.Results.filename;
M = iP.Results.M;
rowHeader = iP.Results.RowHeader;
columnHeader = iP.Results.ColumnHeader;

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

