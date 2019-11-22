function hpc_test_parallel_matlab_unix_neuron(jobnum)
%% Tests whether MATLAB can call NEURON with the unix command under parallel loop on a high performance computing server
%
% Requires:
%       cd/hpc_test_matlab_unix_neuron.m

% File History:
% 2018-03-02 Modified from parallel_test_Matlab_unix_NEURON_test20_Rivanna.m

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

% Run hpc_test_matlab_unix_neuron.m
hpc_test_matlab_unix_neuron(jobnum);
