function timeStamp = create_time_stamp (varargin)
%% Creates a time stamp (default format yyyymmddTHHMM)
% Usage: timeStamp = create_time_stamp (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       timeStamp   - a string for the current date and time
%                   specified as a character array
%
% Arguments:
%       varargin    - 'FormatOut': format of the output text
%                   must be a string scalar or character vector
%                   default == 'yyyymmddTHHMM'
%
%
% Used by:    
%       cd/archive_dependent_scripts.m
%       cd/atfwrite.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_rank_neurons.m
%       cd/m3ha_simulate_population.m
%       cd/parse_all_abfs.m
%       cd/parse_multiunit.m
%       cd/plot_psth.m
%       cd/save_params.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting54.m

% File History:
% 2018-10-21 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
formatOutDefault = 'yyyymmddTHHMM'; % default output format

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FormatOut', formatOutDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, varargin{:});
formatOut = iP.Results.FormatOut;

%% Do the job
timeStamp = datestr(now, formatOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
