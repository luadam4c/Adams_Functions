function homeDirectory = m3ha_locate_homedir
%% Locate the first home directory that exists for the GAT blockade project
% Usage: homeDirectory = m3ha_locate_homedir
% Outputs:
%       homeDirectory   - the first home directory that exists
%                       specified as a character vector
% Arguments:    
%
% Requires:
%       cd/locate_dir.m
%
% Used by:
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_load_sweep_info.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/parse_lts.m
%       ~/m3ha/data_dclamp/dclampDataExtractor.m
%       ~/m3ha/data_dclamp/dclampdatalog_analyze.m
%       ~/m3ha/data_dclamp/dclampPassiveFitter.m
%       ~/m3ha/data_dclamp/PlotHistogramsRefineThreshold.m
%       ~/m3ha/optimizer4gabab/singleneuronfitting42.m and beyond

% File History:
% 2018-10-04 Created by Adam Lu
% 

%% Hard-coded parameters
homeDirectoryCandidates = {'/tmp/data/m3ha/', '/media/adamX/m3ha/', ...
                          '/home/adam/m3ha/', '/scratch/al4ng/m3ha/', ...
                          '/sfs/lustre/scratch/al4ng/m3ha/', ...
                          '/home/shared/mcm/'};
directoryType = 'home directory';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Locate home directory
homeDirectory = locate_dir(homeDirectoryCandidates, ...
                          'DirectoryType', directoryType);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
