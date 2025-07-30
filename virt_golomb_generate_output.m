% Hard-coded parameters
folder1 = 'pb_v_rast_intra';
folder2 = 'pb_v_rast_no_intra';

% Save parent directory
parentDir = pwd;

% Check if directory contains datfig
if ~isfolder('datfig')
    disp("Please run this code with the 'datfig' in current directory");
    return
end

% Change directory to datfig
cd datfig

% Generate input files
command0 = sprintf('julia %s\\genfig\\gen_input_files.jl', parentDir);
system(command0);

% Change directory to dir_here
cd dir_here

% Data for generating figures 6A,B
command1 = sprintf('julia %s\\prog\\irt.jl %s', parentDir, folder1);
command2 = sprintf('julia %s\\prog\\irt.jl %s', parentDir, folder2);
system(command1);
system(command2);


