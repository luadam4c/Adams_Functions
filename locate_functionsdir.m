function functionsDirectory = locate_functionsdir
%% Locate the first functions directory that exists for the GAT blockade project
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
%       cd/combine_sweeps.m
%       cd/find_passive_params.m
%       cd/parse_abf.m
%

% File History:
% 2018-10-04 Created by Adam Lu
% 

%% Hard-coded parameters
functionsDirectoryCandidates = {pwd, '/home/Matlab/', '/scratch/al4ng/Matlab/'};
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