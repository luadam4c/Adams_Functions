function scalar = apply_iteratively (myFunction, array, varargin)
%% Applies a function iteratively to an array until it becomes a non-cell array scalar
% Usage: scalar = apply_iteratively (myFunction, array, varargin)
% Explanation:
%       Applies a function iteratively to an array.
%           The function must be able to take elements of the array as an argument
%               and return outputs that can be retaken as input
% Example(s):
%       a = apply_iteratively(@max, magic(3))
%       b = apply_iteratively(@min, {1:10, -10:5, 5:30})
%       c = apply_iteratively(@max, {1:10, -10:5, 5:30})
%       d = apply_iteratively(@max, {magic(3), -10:5})
%       e = apply_iteratively(@max, {{magic(3), magic(3)}, {magic(3)}})
%       f = apply_iteratively(@max, {{magic(3), magic(3)}, {[], []}})
% Outputs:
%       scalar      - the resulting scalar
%                   specified as a scalar
% Arguments:
%       myFunction  - a custom function
%                   must be a function handle
%       array       - an array to apply the function iteratively
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_axis_limits.m
%       cd/identify_channels.m
%       cd/plot_traces.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-28 Fixed code logic for case e in examples
% TODO: Add 'TreatCellAsArray' as an optional argument
% 
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

% Add required inputs to the Input Parser
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, myFunction, array, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO

%% Do the job
% Apply myFunction to array until it becomes a non-cell scalar
while iscell(array) || ~isscalar(array)
    if isempty(array)
        % Change to NaN
        array = NaN;
    elseif iscell(array)
        % If any element is empty, remove it
        array = array(~cellfun(@isempty, array));

        % Reduce the cell array to an ordinary array
        try
            % First try if one can reduce a cell array to an ordinary array
            %   by applying myFunction to each cell
            array = cellfun(myFunction, array);
        catch
            % If not, apply myFunction iteratively in each cell
            array = cellfun(@(x) apply_iteratively(myFunction, x), array, ...
                            'UniformOutput', false);
        end
    else
        % Apply myFunction to array until it becomes a scalar
        array = myFunction(array);
    end
end

% The final array should be a non-cell scalar
scalar = array;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

array = cellfun(myFunction, array, 'UniformOutput', false);

% If the final array is empty, change it to NaN
if isempty(array)
    scalar = NaN;
else
    scalar = array;
end

% Apply myFunction to array until it becomes a non-cell scalar or empty
while iscell(array) || ~(isscalar(array) || isempty(array))

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%