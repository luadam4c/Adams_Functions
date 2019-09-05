function allData = combine_multiunit_data (varargin)
%% Combines data across multiple .abf files (using multiunit data defaults) for each slice in the input folder (or for a particular slice)
% Usage: allData = combine_multiunit_data (varargin)
% Explanation:
%       This is the same as combine_data_from_same_slice.m
%           except for the following defaults:
%               'ChannelTypes' - {'voltage', 'current'}
%               'ChannelUnits' - {'uV', 'arb'}
%
% Example(s):
%       allData = combine_multiunit_data;
%       allData = combine_multiunit_data('SliceBase', 'slice4');
%
% Outputs:
%       allData     - combined data for all slices (or for a particular slice)
%                   specified as a table (or struct)
%
% Arguments:
%       varargin    - 'ChannelTypes': the channel types
%                   must be a cellstr with nChannels elements
%                       each being one of the following:
%                           'Voltage'       - voltage
%                           'Current'       - current
%                           'Conductance'   - conductance
%                           'Other'         - other un-identified types
%                   default == detected with identify_channels()
%                   - 'ChannelUnits': the channel units
%                   must be a cellstr with nChannels elements
%                   default == detected with identify_channels()
%                   - Any other parameter-value pair for 
%                           combine_data_from_same_slice() or 
%                           combine_abf_data() or parse_all_abfs()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/struct2arglist.m
%       cd/combine_data_from_same_slice.m
%
% Used by:
%       cd/parse_all_multiunit.m
%       cd/parse_multiunit.m

% File History:
% 2019-08-23 Created by Adam Lu
% 

%% Hard-coded parameters
validChannelTypes = {'Voltage', 'Current', 'Conductance', 'Other'};

%% Default values for optional arguments
channelTypesDefault = {'voltage', 'current'};
channelUnitsDefault = {'uV', 'arb'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ChannelTypes', channelTypesDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));
addParameter(iP, 'ChannelUnits', channelUnitsDefault, ...
    @(x) isempty(x) || isstring(x) || iscellstr(x));

% Read from the Input Parser
parse(iP, varargin{:});
channelTypes = iP.Results.ChannelTypes;
channelUnits = iP.Results.ChannelUnits;

% Keep unmatched arguments for combine_data_from_same_slice() or 
%                           combine_abf_data() or parse_all_abfs()
otherArguments = struct2arglist(iP.Unmatched);

% Validate channel types
if ~isempty(channelTypes)
    channelTypes = cellfun(@(x) validatestring(x, validChannelTypes), ...
                            channelTypes, 'UniformOutput', false);
end

%% Do the job
allData = combine_data_from_same_slice('ChannelTypes', channelTypes, ...
                            'ChannelUnits', channelUnits, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%