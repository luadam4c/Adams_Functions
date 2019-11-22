function hpc_test_matlab_figures(jobnum)
%% Test whether MATLAB figures can be suppressed on a high performance computing server
% 
% Used by:
%       cd/hpc_test_parallel_matlab_figures.m

% File History:
% File History:
%   2018-03-08 Created
%   2018-03-08 Changed nFigs to 1
%   2018-03-08 Changed nFigs back to 10
%   2018-03-08 Changed title to suptitle (what seems to be causing the problem)

% Hard-coded parameters
nFigs = 10;
nNums = 10;
outFolder = 'hpc_test_matlab_figures_output';
titleBase = 'Random figure';
fignameBase = 'randfig';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Print a starting message
fprintf('Starting job number #%d ...\n\n', jobnum);

%% Create output folder if it doesn't exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
end

%% Plot figures
h = cell(nFigs, 1);
%##########
%##############
parfor k = 1:nFigs
    h{k} = figure('Visible', 'off');
    hold on;
    plot(1:nNums, rand(1, nNums), 'r');
end
%##############
%##########

%% Plot more and save figures
parfor k = 1:nFigs
    id = [num2str(jobnum), '-', num2str(k)];
    figure(h{k});
    plot(1:nNums, rand(1, nNums), 'g');
    suptitle([titleBase, ' #', id]);
    saveas(h{k}, fullfile(outFolder, [fignameBase, '_', id]), 'png');
end

%% Close all figures
close all force hidden

%% Print a completion message
fprintf('Finished job number #%d!\n\n', jobnum);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

nFigs = 10;
title([titleBase, ' #', id]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%