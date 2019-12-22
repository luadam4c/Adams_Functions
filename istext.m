function isText = istext (x, varargin)
%% Returns whether the input is a character array, a string array or a cell array of character arrays
% Usage: isText = istext (x, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       isText      - whether the input is a text
%                   specified as a logical scalar
% Arguments:
%       x           - input to test
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/arglist2struct.m
%       cd/combine_data_from_same_slice.m
%       cd/ismatch.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_neuron_create_initial_params.m
%       cd/m3ha_select_raw_traces.m

% File History:
% 2019-01-11 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

%% Do the job
isText = ischar(x) || isstring(x) || iscellstr(x);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%