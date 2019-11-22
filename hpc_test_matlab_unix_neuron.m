function hpc_test_matlab_unix_neuron(jobnum)
%% Tests whether MATLAB can call NEURON with the unix command
% 
% Used by:
%       cd/hpc_test_parallel_matlab_unix_neuron.m

% File History:
%   2018-03-02 Created

%% Hard-coded parameters
% The number of simulations to perform
numswps = 19;

% The command to simulate
simCommand = sprintf([ ...
    'build("passive", 38.42, 97.16, 0.2272, 50)\n', ...
    'adjust_globalpas(0.88, 173, 7.954)\n', ...
    'adjust_leak(1e-05, -80, 7.954)\n', ...
    'sim("passive", "test_here.out", 2260, -65.9692, ', ...
        '-0.0506104, 0.016, 52, 90.1, 1073.2, 0.952)\n']);

% The command for running NEURON
runNeuronCommand = 'x86_64/special singleneuron4compgabab.hoc';

% The command on a high performance computing server 
%   for loading modules required for NEURON to work
moduleLoadCommands = sprintf([ ...
    '# Refresh bash shell to load environmental variables\n', ...    
    'source /etc/bashrc\n', ...
    '# Purge all loaded modules\n', ...
    'module purge\n', ...
    '# Load the newest version of the Intel compiler (16.0)\n', ...
    '#   This is necessary before loading the Open MPI library\n', ...
    'module load intel\n', ...
    '# Load the newest version of the Open MPI library (2.1.1)\n', ...
    'module load openmpi\n', ...
    '# Load the newest version of NEURON (7.4)\n', ...
    'module load neuron\n']);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Print a starting message
fprintf('Starting job number #%d ...\n\n', jobnum);

%% RUN NEURON
results = cell(1, numswps);            % stores simulation standard outputs
tstart_NEURON = tic();
%##########
%##############
parfor k = 1:numswps
    [status, results{k}] = ...
        unix(sprintf(['%s\n', '%s - << here\n', ...
                      '%s\n', 'print "No_Errors!"\n', 'here'], ...
                     moduleLoadCommands, runNeuronCommand, simCommand));

    % Save simulation standard output in a text file
    fid = fopen([mfilename, '_job', num2str(jobnum), ...
                '_swp', num2str(k), '_output.txt'], 'w');
    fprintf(fid, ['Return status was: %d\n\n', ...
                  'Simulation output was:\n\n%s\n'], ...
                  status, results{k});
    fclose(fid);
end
%##############
%##########

%% Print a completion message
fprintf('Finished job number #%d!\n\n', jobnum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% moduleLoadCommands = '';

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%