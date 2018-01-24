function print_next_in_csv(fid, i, n)
% What to print next in a csv file
% 
% Used by:
%       /home/Matlab/Adams_Functions/log_matfile.m
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

if i == n
    fprintf(fid, '\n');
else
    fprintf(fid, ', ');
end        

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

