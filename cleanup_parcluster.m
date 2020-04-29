%% Cleans up parallel cluster, removing all jobs that contain crash dump files
% Used by:
%       cd/array_fun.m

% File History:
% 2019-03-13 Created

% Shutdown parallel pool
delete(gcp('nocreate'));

% Create parallel cluster
myCluster = parcluster('local');

% Remove all jobs created with profile local
delete(myCluster.Jobs);