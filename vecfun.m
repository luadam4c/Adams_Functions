function outVecs = vecfun (myFunction, inVecs, varargin)
%% Apply a function to each vector (each column of an array or each element of a cell array of vectors)
% Usage: outVecs = vecfun (myFunction, inVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       vecfun(@min, {1:5, 3:6})
%       vecfun(@smooth, magic(5))
%
% Outputs:
%       outVecs    - output array
%                   specified as an array
%
% Arguments:
%       myFunction  - a custom function
%                   must be a function handle
%       inVecs     - input vectors
%                   must be an array or a cell array of vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for cellfun() or arrayfun()
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/medianfilter.m
%       cd/movingaveragefilter.m

% File History:
% 2019-08-29 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction', ...               % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'2d'}));
addRequired(iP, 'inVecs');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, myFunction, inVecs, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the cellfun() or arrayfun() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
if iscellvector(inVecs)
    % Apply myFunction to each cell and always return a cell array
    outVecs = array_fun(myFunction, inVecs, 'UniformOutput', false, ...
                        otherArguments{:});
else
    % Count the number of columns
    nColumns = size(inVecs, 2);

    % Apply the myFunction to each column and return things in a cell array
    outVecsCell = array_fun(@(x) myFunction(inVecs(:, x)), ...
                        transpose(1:nColumns), 'UniformOutput', false, ...
                        otherArguments{:});

    % Force as a non-cell matrix
    outVecs = force_matrix(outVecsCell, 'TreatCellAsArray', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%