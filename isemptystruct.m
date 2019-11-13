function isEmptyStruct = isemptystruct (structArray)
%% Returns whether a structure has no fields
% Usage: isEmptyStruct = isemptystruct (structArray)
% Explanation:
%       TODO
%
% Example(s):
%       isemptystruct(struct)
%
% Outputs:
%       isEmptyStruct   - whether the structure array has no fields
%                       specified as a logical scalar
%
% Arguments:
%       structArray - a structure array
%                   must be a struct array
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/m3ha_optimizer_4compgabab.m

% File History:
% 2019-11-13 Created by Adam Lu
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'structArray');

% Read from the Input Parser
parse(iP, structArray, varargin{:});

%% Do the job
isEmptyStruct = isempty(fieldnames(structArray));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%