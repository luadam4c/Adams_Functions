function timeStamp = create_time_stamp (varargin)
%% Creates a time stamp (default format yyyyMMddTHHmm)
% Usage: timeStamp = create_time_stamp (varargin)
% Explanation:
%       Creates a time stamp by passing argument to either datetime.m (if
%       R2014b or later) versus datestr.m
%
% Example(s):
%       timeStamp = create_time_stamp
%       timeStamp = create_time_stamp('FormatOut', 'yyyyMMdd''T''HHmm')
%       timeStamp = create_time_stamp('FormatOut', 'yyyymmdd')
%       timeStamp = create_time_stamp('FormatOut', 'yyyymmddTHHMM')
%
% Outputs:
%       timeStamp   - a string for the current date and time
%                   specified as a character array
%
% Arguments:
%       varargin    - 'FormatOut': format of the output text
%                   must be a string scalar or character vector
%                   default == 'yyyymmddTHHMM'
%                   - 'UseDateStr': forces use of datestr function instead
%                   default == false
%
%
% Used by:    
%       cd/archive_dependent_scripts.m
%       cd/write_data_atf.m
%       cd/m3ha_network_launch.m
%       cd/m3ha_rank_neurons.m
%       cd/m3ha_simulate_population.m
%       cd/parse_all_abfs.m
%       cd/parse_multiunit.m
%       cd/plot_psth.m
%       cd/save_params.m
%       /media/adamX/m3ha/optimizer4gabab/singleneuronfitting54.m
%       \Shared\Code\vIRt\virt_moore.m
%       \Shared\Code\vIRt\virt_plot_whisk_analysis.m

% File History:
% 2018-10-21 Created by Adam Lu
% 2025-08-01 Update to use datetime.m and added 'UseDateStr' as parameter
%               to ensure backwards compatibility
% 2025-08-21 Fixed to output a character array

%% Hard-coded parameters

%% Default values for optional arguments
formatOutDefault = 'yyyyMMdd''T''HHmm'; % default output format
useDateStrDefault = false;              % use datetime.m by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FormatOut', formatOutDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
addParameter(iP, 'UseDateStr', useDateStrDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));                                                % introduced after R2016b

% Read from the Input Parser
parse(iP, varargin{:});
formatOut = iP.Results.FormatOut;
useDateStr = iP.Results.UseDateStr;

%% Check if format is for datestr function
if contains(formatOut, 'yyyymmdd')
    % Use datestr function
    useDateStr = true;
end

%% Do the job
if ~useDateStr && exist('datetime.m', 'file') == 2
    % Since R2014b
    timeStamp = char(datetime('now', 'Format', formatOut));
else            
    % Introduced before R2006a, not recommended since R2022b
    timeStamp = datestr(now, formatOut);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
