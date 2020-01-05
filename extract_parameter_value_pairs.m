function [paramsStruct, argList] = extract_parameter_value_pairs (argList)
%% Extracts parameter-value pairs from an argument list and return remaining
% Usage: [paramsStruct, argList] = extract_parameter_value_pairs (argList)
% Explanation:
%       TODO
%
% Example(s):
%       [a, b] = extract_parameter_value_pairs({1:5, 2:6, 'Name', 3})
%       [params, varargin] = extract_parameter_value_pairs(varargin)
%
% Outputs:
%       paramsStruct    - structure of extract parameters
%                       specified as a scalar structure
%       argList     - a cell array of remaining arguments
%                   must be a cell array
%
% Arguments:
%       argList     - a cell array of arguments
%                   must be a cell array
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/array_fun.m
%       cd/match_format_vectors.m

% File History:
% 2020-01-04 Moved from array_fun.m
% TODO: Right now this assumes parameters are char but not string

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
addRequired(iP, 'argList', ...
    @(x) validateattributes(x, {'cell'}, {'2d'}));

% Read from the Input Parser
parse(iP, argList);

%% Do the job
% Initialize starting index of parameter-value pairs
idxParamsStart = NaN;

% Initialize content tested
contentTested = argList;

% While there is still at least two element left, examine in pairs
while numel(contentTested) >= 2
    if ischar(contentTested{end - 1})
        idxParamsStart = numel(contentTested) - 1;
    end

    contentTested = contentTested(1:end-2);
end

% Divide up the content
if ~isnan(idxParamsStart)
    % Create a parameters list
    paramsList = argList(idxParamsStart:end);

    % Output as a parameters structure
    paramsStruct = arglist2struct(paramsList);

    % Truncate the original argument list
    argList = argList(1:idxParamsStart - 1);
else
    paramsStruct = struct.empty;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%