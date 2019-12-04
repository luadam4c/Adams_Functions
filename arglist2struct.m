function structure = arglist2struct (argList, varargin)
%% Converts an argument list tp a scalar structure
% Usage: structure = arglist2struct (argList, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       arglist2struct({'a', 1, 'b', 2})
%       arglist2struct({'a', 'b', 'c', 'd'})
%
% Outputs:
%       structure   - structure with field names as parameters
%                   must be a scalar structure
%
% Arguments:
%       argList     - argument list (parameter value pairs)
%                   specified as a cell array of character vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/istext.m
%
% Used by:
%       cd/extract_vars.m

% File History:
% 2018-12-28 Moved from annotation_in_plot.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, argList, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Initialize structure
structure = struct;

% If empty, return empty structure
if isempty(argList)
    return
end

% Make sure there are an even number of elements
if mod(numel(argList), 2) ~= 0
    error('Argument lists must have an even number of elements!');
end

% Reshape as two rows, then transpose, so that 
%   names are in the 1st column and values in the 2nd column
nameValuePairs = transpose(reshape(argList, 2, []));

% The first column are the names
fieldNames = nameValuePairs(:, 1);

% Make sure all names are text
if ~istext(fieldNames)
    error('Names must be text!');
end

% The second column are the values
fieldValues = nameValuePairs(:, 2);

% Convert to a structure
structure = cell2struct(fieldValues, fieldNames);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%