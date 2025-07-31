function functionsDirectory = locate_functionsdir
%% Locate the first shared functions directory that exists
% Usage: functionsDirectory = locate_functionsdir
% Outputs:
%       functionsDirectory  - the first functions directory that exists
%                               and contains Adams_Functions
%                           specified as a character vector
% Arguments:    
%
% Requires:
%       cd/locate_dir.m
%
% Used by:
%       cd/char2rgb.m
%       cd/combine_sweeps.m
%       cd/create_colormap.m
%       cd/find_passive_params.m
%       cd/minEASE.m
%       cd/parse_abf.m
%       cd/plot_grouped_jitter.m
%       cd/plot_violin.m
%       ~/m3ha/data_dclamp/dclampPassiveFitter.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m and beyond
%       ~/m3ha/optimizer4gabab/m3ha_optimizer_4compgabab.m
%       ~/FluoroSNNAP/FluroSNNAP.m
%

% File History:
% 2018-10-04 Created by Adam Lu
% 2019-05-16 Added MatlabFishFish to functions directory candidates
% 

%% Hard-coded parameters
functionsDirectoryCandidates = {pwd, '/home/Matlab/', ...
    '/scratch/al4ng/Matlab/', ...
    '/scratch/kas5dv/Matlab/', ...
    '/sfs/lustre/scratch/al4ng/Matlab/', ...
    };
directoryType = 'functions directory';
containedSubdir = 'Adams_Functions';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Locate functions directory
functionsDirectory = locate_dir(functionsDirectoryCandidates, ...
                          'DirectoryType', directoryType, ...
                          'ContainedSubdir', containedSubdir);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
