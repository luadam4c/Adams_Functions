function csvwrite_with_header (filename, M, header)
%% Write a comma-separated value file with given header
% Usage: csvwrite_with_header (filename, M, header)
% TODO
%
% Used by:    
%       /home/Matlab/minEASE/examine_gapfree_events.m
%
% File History:
% 2017-07-24 Moved from examine_gapfree_events.m
% TODO: Input Parser
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Open file for writing
fid = fopen(filename, 'w');

% Print header
for j = 1:numel(header)
    fprintf(fid, '%s, ', header{j});
end
fprintf(fid, '\n');

% Print matrix contents
if ~isempty(M)
    % Make sure header is the correct length
    if numel(header) ~= size(M, 2)
        error('Header does not have the correct length!');
    end

    % Print matrix contents
    for i = 1:size(M, 1)
        for j = 1:size(M, 2)
            fprintf(fid, '%g, ', M(i, j));
        end
        fprintf(fid, '\n');
    end
end

% Close file
fclose(fid);

