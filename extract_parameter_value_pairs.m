function [paramsOutput, argList] = extract_parameter_value_pairs (argList, varargin)
%% Extracts parameter-value pairs from an argument list and return remaining
% Usage: [paramsOutput, argList] = extract_parameter_value_pairs (argList, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [a, b] = extract_parameter_value_pairs({1:5, 2:6, 'Name', 3})
%       [params, varargin] = extract_parameter_value_pairs(varargin)
%
% Outputs:
%       paramsOutput    - structure or cell array of extract parameters
%                       specified as a scalar structure or cell array
%       argList     - a cell array of remaining arguments
%                   specified as a cell array
%
% Arguments:
%       argList     - a cell array of arguments
%                   must be a cell array
%       varargin    - 'OutputMode': whether to output params as 
%                                   a structure or a cell array
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'struct' - scalar structure
%                       'cell'   - cell array
%                   default == 'struct'
%
% Requires:
%       cd/arglist2struct.m
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/array_fun.m
%       cd/match_format_vectors.m
%       cd/write_table.m

% File History:
% 2020-01-04 Moved from array_fun.m
% 2025-09-12 Added 'OutputMode' as an optional argument
% TODO: Right now this assumes parameters are char but not string

%% Hard-coded parameters
validOutputModes = {'struct', 'cell'};

%% Default values for optional arguments
outputModeDefault = 'struct';       % default: returns params as a struct

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

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputMode', outputModeDefault, ...
    @(x) any(validatestring(x, validOutputModes)));

% Read from the Input Parser
parse(iP, argList, varargin{:});
outputMode = validatestring(iP.Results.OutputMode, validOutputModes);

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

    % Output parameters
    switch outputMode
        case 'struct'
            paramsOutput = arglist2struct(paramsList);
        case 'cell'
            paramsOutput = paramsList;
    end

    % Truncate the original argument list
    argList = argList(1:idxParamsStart - 1);
else
    switch outputMode
        case 'struct'
            paramsOutput = struct.empty;
        case 'cell'
            paramsOutput = {};
    end    
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%