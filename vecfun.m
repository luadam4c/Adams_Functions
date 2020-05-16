function varargout = vecfun (myFunction, inVecs, varargin)
%% Apply a function to each vector (each column of an array or each element of a cell array of vectors)
% Usage: varargout = vecfun (myFunction, inVecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       vecfun(@min, {1:5, 3:6})
%       vecfun(@smooth, magic(5))
%       vecfun(@max, magic(5))
%       vecfun(@max, 1:5)
%
% Outputs:
%       outVecs    - output array
%                   specified as an array
%
% Arguments:
%       myFunction  - a custom function
%                   must be a function handle
%       inVecs      - input vectors
%                   must be an array or a cell array of vectors
%       varargin    - 'UniformOutput': whether outputs can be concatenated
%                                       as non-cell arrays if the input
%                                       is a cell array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for array_fun()
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/iscellvector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_pairwise_differences.m
%       cd/compute_sampling_interval.m
%       cd/medianfilter.m
%       cd/movingaveragefilter.m
%       cd/parse_ipsc.m

% File History:
% 2019-08-29 Created by Adam Lu
% 2020-04-16 Now forces array outputs as column vectors by default
% 2020-04-16 Now tries uniform output first
% 2020-04-16 Now allows multiple outputs
% 2020-05-16 Default UniformOutput is now true

%% Hard-coded parameters

%% Default values for optional arguments
uniformOutputDefault = true;

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
addParameter(iP, 'UniformOutput', uniformOutputDefault);

% Read from the Input Parser
parse(iP, myFunction, inVecs, varargin{:});
uniformOutput = iP.Results.UniformOutput;

% Keep unmatched arguments for the cellfun() or arrayfun() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
if ~iscell(inVecs)

end

%% Do the job
if iscellvector(inVecs)
    % Apply myFunction to each cell
    [varargout{1:nargout}] = ...
        array_fun(myFunction, inVecs, 'UniformOutput', uniformOutput, ...
                    otherArguments{:});
else
    % Count the number of columns
    nColumns = size(inVecs, 2);

    % Apply myFunction to each column and return things in a cell array
    [varargout{1:nargout}] = ...
        array_fun(@(x) myFunction(inVecs(:, x)), ...
                        transpose(1:nColumns), 'UniformOutput', false, ...
                        otherArguments{:});

    % Force as a non-cell matrix
    varargout = cellfun(@(x) force_matrix(x, 'TreatCellAsArray', false), ...
                        varargout, 'UniformOutput', false);    

    % Force as a column vector and ignore matrices
    varargout = ...
        cellfun(@(x) force_column_vector(x, 'IgnoreNonvectors', true), ...
                varargout, 'UniformOutput', false); 
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%