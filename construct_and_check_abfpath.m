function [abfPath, pathExists] = construct_and_check_abfpath (abfFileName, varargin)
%% Constructs the full path to a .abf file and checks whether it exists
% Usage: [abfPath, pathExists] = construct_and_check_abfpath (abfFileName, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       abfPath     - the full path(s) to .abf file(s) constructed
%                   specified as a character vector 
%                       or a column cell array or character vectors
%       pathExists  - whether the path(s) exists
%                   specified as a column logical array
% Arguments:
%       abfFileName - file or directory name(s)
%                       e.g. 'A100110_0008_18.mat'
%                       e.g. {'folder1', 'folder2'}
%                   must be a string/character array or a cell array 
%                       of strings/character arrays
%       varargin    - whatever accepted by construct_and_check_fullpath.m
%
% Requires:
%       cd/construct_and_check_fullpath.m
%
% Used by:
%       cd/parse_abf.m
%       cd/parse_assyst_swd.m
%       cd/parse_current_family.m

% File History:
% 2018-11-21 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'abfFileName', ...
    @(x) assert(ischar(x) || iscellstr(x) || isstring(x), ...
                ['pathName must be either a string/character array ', ...
                    'or a cell array of strings/character arrays!']));

% Read from the Input Parser
parse(iP, abfFileName);

%% Do the job
[abfPath, pathExists] = ...
    construct_and_check_fullpath(abfFileName, 'Extension', '.abf', varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
