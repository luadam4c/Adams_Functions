function result = apply_iteratively (myFunction, array, varargin)
%% Applies a function iteratively to an array until it becomes a non-cell array result
% Usage: result = apply_iteratively (myFunction, array, varargin)
% Explanation:
%       Applies a function iteratively to an array.
%           The function must be able to take elements of the array as an argument
%               and return outputs that can be retaken as input
%
% Example(s):
%       a = apply_iteratively(@max, magic(3))
%       b = apply_iteratively(@min, {1:10, -10:5, 5:30})
%       c = apply_iteratively(@max, {1:10, -10:5, 5:30})
%       d = apply_iteratively(@max, {magic(3), -10:5})
%       e = apply_iteratively(@max, {{magic(3), magic(3)}, {magic(3)}})
%       f = apply_iteratively(@max, {{magic(3), magic(3)}, {[], []}})
%       g = apply_iteratively(@unique, {{magic(3), magic(3)}, {[], []}})
%       h = apply_iteratively(@unique_custom, {{magic(3), magic(3)}, {[], []}}, 'IgnoreNaN', true)
%       i = apply_iteratively(@union_over_cells, {{1:2, 3:5}, {1:3}})
%
% Outputs:
%       result      - the resulting vector
%                   specified as a vector
% Arguments:
%       myFunction  - a custom function
%                   must be a function handle
%       array       - an array to apply the function iteratively
%                   must be an array
%       varargin    - 'OptArg': optional argument that's 
%                                   not a parameter-value pair 
%                   default == []
%                   - optional parameter-value pairs for myFunction
%
% Requires:
%       cd/create_error_for_nargin.m
%
% Used by:
%       cd/compute_axis_limits.m
%       cd/create_default_grouping.m
%       cd/extract_subvectors.m
%       cd/identify_channels.m
%       cd/m3ha_network_show_net.m
%       cd/plot_raster.m
%       cd/plot_relative_events.m
%       cd/plot_traces.m
%       cd/plot_swd_histogram.m

% File History:
% 2018-12-19 Created by Adam Lu
% 2018-12-28 Fixed code logic for case e in examples
% 2019-04-24 Allowed the result to be a non-scalar vector
% TODO: Add 'TreatCellAsArray' as an optional argument
% 
% 

%% Hard-coded parameters
debugFlag = false;

%% Default values for optional arguments
optArgDefault = '';         % Not provided

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
addRequired(iP, 'myFunction', ...           % a custom function
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));
addRequired(iP, 'array');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OptArg', optArgDefault);

% Read from the Input Parser
parse(iP, myFunction, array, varargin{:});
optArg = iP.Results.OptArg;

% Keep unmatched arguments for myFunction
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Apply myFunction to array until it becomes a non-cell scalar
%   or if it does not change from the previous iteration
arrayOld = NaN;
while iscell(array) || ...
        ~(isscalar(array) || isequaln(array, arrayOld))
    % Update arrayOld
    arrayOld = array;

    % Apply function to array
    if isempty(array)
        % Change to NaN
        array = NaN;
    elseif iscell(array)
        % If any element is empty, remove it
        array = array(~cellfun(@isempty, array));

        % Save the old array
        arrayOld = array;

        % Reduce the cell array to an ordinary array
        try
            % First try if one can reduce a cell array to an ordinary array
            %   by applying myFunction to each cell
            if ~isempty(optArg)
                array = cellfun(@(x) myFunction(x, optArg, ...
                                                otherArguments{:}), array);
            else
                array = cellfun(@(x) myFunction(x, otherArguments{:}), array);
            end
            if debugFlag
                disp('Applied to each cell')
            end
        catch
            % If not, apply myFunction iteratively in each cell
            array = cellfun(@(x) apply_iteratively(myFunction, x, ...
                                    'OptArg', optArg, otherArguments{:}), ...
                            array, 'UniformOutput', false);
            if debugFlag
                disp('Applied iteratively to each cell')
            end
        end

        % If the array still hasn't changed, concatenate the vectors
        if isequaln(array, arrayOld)
            try
                array = horzcat(array{:});
            catch
                array = vertcat(array{:});
            end
        end
    else
        % Apply myFunction to array until it becomes a vector
        if ~isempty(optArg)
            array = myFunction(array, optArg, otherArguments{:});
        else
            array = myFunction(array, otherArguments{:});
        end
    end
    
    if debugFlag
        disp('Old array: ');
        disp(arrayOld);
        disp('New array: ');
        disp(array);
    end
end

% The final array should be a non-cell vector
result = array;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

while iscell(array) || ~isscalar(array)
array = myFunction(array);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%