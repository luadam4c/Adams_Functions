function vecsFilt = movingaveragefilter (vecs, varargin)
%% Applies a moving average filter to vectors
% Usage: vecsFilt = movingaveragefilter (vecs, filtWidth, si, varargin)
% Explanation:
%       Same as smooth() but with option of using a window in time units
%
% Example(s):
%       movingaveragefilter(magic(5))
%       movingaveragefilter(magic(5), 10, 10)
%       movingaveragefilter(magic(5), 10, 2)
%       movingaveragefilter(magic(5), 10, 1)
%
% Outputs:
%       vecsFilt    - filtered vector(s)
%                   specified as a numeric array
%
% Arguments:
%       vecs        - vector(s) to moving average filter
%                   must be a numeric array
%       filtWidth   - (opt) filter window width in time units
%                           must be a numeric vector
%       si          - (opt) sampling interval in time units
%                           must be a numeric vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for the smooth() function
%
% Requires:
%       cd/apply_to_nonnan_part.m
%       cd/create_error_for_nargin.m
%       cd/find_nearest_odd.m
%       cd/force_column_cell.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%       cd/vecfun.m
%
% Used by:
%       cd/compute_autocorrelogram.m
%       cd/m3ha_plot_simulated_traces.m
%       cd/parse_ipsc.m
%       cd/parse_lts.m

% File History:
% 2019-01-14 Created by Adam Lu
% 2019-08-29 Now uses vecfun.m
% 2019-11-18 Fixed the case with leading and trailing NaNs
% 

%% Hard-coded parameters

%% Default values for optional arguments
filtWidthDefault = 3;       % 3 sample points by default
siDefault = 1;              % treat filtWidth as samples by default
% param1Default = [];       % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vec', ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));

% Add optional inputs to the Input Parser
addOptional(iP, 'filtWidth', filtWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addOptional(iP, 'si', siDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vecs, varargin{:});
filtWidth = iP.Results.filtWidth;
si = iP.Results.si;
% param1 = iP.Results.param1;

% Keep unmatched arguments for the smooth() function
otherArguments = struct2arglist(iP.Unmatched);

%% Do the job
% Calculate the moving average filter window width in samples
%   Note: Round down to the nearest odd integer to preserve values!!
%           However, must be >= 1
filtWidthSamples = find_nearest_odd(filtWidth ./ si, 'Direction', 'down');

% Moving average filter vectors
if isscalar(filtWidthSamples)
    vecsFilt = vecfun(@(x) movingaveragefilter_helper(x, filtWidthSamples, ...
                                                    otherArguments), ...
                        vecs, 'UniformOutput', false);
else
    % Match formats
    vecs = force_column_cell(vecs);
    filtWidthSamples = force_column_vector(filtWidthSamples);

    % Apply a different filter window to each vector
    vecsFilt = cellfun(@(x, y) movingaveragefilter_helper(x, y, ...
                                                    otherArguments), ...
                        vecs, num2cell(filtWidthSamples), ...
                        'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function vecFilt = movingaveragefilter_helper(vec, filtWidthSamples, ...
                                                otherArguments)
%% Applies a moving average filter to one vector

% Create a custom smooth function
myFun = @(x) smooth(x, filtWidthSamples, otherArguments{:});

% Apply the custom smooth function to just the non-NaN part
%   Note: If this filter is applied to leading and trailing NaNs,
%           values become extrapolated, which is most likely not desired
vecFilt = apply_to_nonnan_part(myFun, vec);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%