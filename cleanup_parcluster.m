%% Cleans up parallel cluster, removing all jobs that contain crash dump files
% Shutdown parallel pool
delete(gcp('nocreate'));

% Create parallel cluster
myCluster = parcluster('local');

% Remove all jobs created with profile local
delete(myCluster.Jobs);