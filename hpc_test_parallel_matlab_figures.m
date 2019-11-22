function hpc_test_parallel_matlab_figures(jobnum)
%% Set up parallel pool and run hpc_test_matlab_figures.m on Rivanna
%
% Requires:
%       cd/hpc_test_matlab_figures.m

% File History:
% 2018-03-08 Modified from parallel_test_Matlab_unix_NEURON.m

% Set up MATLAB compiler runtime profile
setmcruserdata('ParallelProfile', 'local.settings');

% Set up Parallel Cluster Profile
pc = parcluster('local');

% Explicitly set the JobStorageLocation to the temporary directory
%	that was created in the SLURM script
pc.JobStorageLocation = ['/scratch/', getenv('USER'), '/', ...
                          getenv('slurmArrayID')];

% Start the parallel pool with the number of workers defined in the SLURM script
parpool(pc, str2num(getenv('numWorkers')));

% Make sure arguments are numbers
if ~isnumeric(jobnum) && ischar(jobnum)
    jobnum = str2num(jobnum);
end

% Run hpc_test_matlab_figures.m
hpc_test_matlab_figures(jobnum);
