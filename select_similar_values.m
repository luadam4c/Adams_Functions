function [valSelected, indSelected] = select_similar_values (values, varargin)
%% Selects values that are within a certain percentage range of the mean
% Usage: [valSelected, indSelected] = select_similar_values (values, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       valSelected - selected values
%                   specified as a numeric vector
%       indSelected - indices of selected values
%                   specified as a positive integer vector
% Arguments:
%       values      - values to select from
%                   must be a numeric vector
%       varargin    - 'NToSelect': number of values to select
%                   must be a positive integer scalar
%                   default == 5
%                   - 'EndPoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], vecs)
%                   - 'Direction': the selection direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'forward'   - select from the first indices
%                       'backward'  - select from the last indices
%                   default == 'forward'
%                   - 'MmaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 20%
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_element.m
%       cd/extract_subvector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_phase_average.m

% File History:
% 2019-05-12 Created by Adam Lu
% 

%% Hard-coded parameters
validDirections = {'forward', 'backward'};

%% Default values for optional arguments
nToSelectDefault = 5;           % select 5 values by default
endPointsDefault = [];          % set later
directionDefault = 'forward';   % select from the first indices by default
maxRange2MeanDefault = 20;      % range is not more than 20% of mean by default

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
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NToSelect', nToSelectDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, values, varargin{:});
nToSelect = iP.Results.NToSelect;
endPoints = iP.Results.EndPoints;
direction = validatestring(iP.Results.Direction, validDirections);
maxRange2Mean = iP.Results.MaxRange2Mean;

% Keep unmatched arguments for the TODO() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Save first index
idxFirst = extract_element(endPoints, 'first');

% Restrict to end points
valuesRestricted = extract_subvector(values, 'EndPoints', endPoints);

% Decide on the direction
switch direction
    case 'forward'
        % Don't flip the values
        valuesOfInterest = valuesRestricted;
    case 'backward'
        % Flip the values
        valuesOfInterest = flip(valuesRestricted);
    otherwise
        error('direction unrecognized!');
end

% Count the number of samples
nValues = count_samples(valuesOfInterest);

% Check if there are enough values to select from
if nToSelect > nValues
    valSelected = NaN;
    indSelected = [];
    return
end

%% Do the job
% Select the initial set of values
indSelected = 1:nToSelect;

% Get the corresponding set of values
valSelected = valuesOfInterest(indSelected);

% Compute the range to mean percentage ratio
meanSelected = mean(valSelected);
rangeSelected = range(valSelected);
range2mean = (rangeSelected / meanSelected) * 100;

% Perform similarity test 
while range2mean > maxRange2Mean
    % Take the most extreme index out
    % TODO

    % Add the next index
    % TODO

    % TODO indSelected = 
    valSelected = valuesOfInterest(indSelected);
    meanSelected = mean(valSelected);
    rangeSelected = range(valSelected);
    range2mean = (rangeSelected / meanSelected) * 100;
end

% Reconstruct original indices
switch direction
    case 'forward'
        indSelected = (idxFirst - 1) + indSelected;
    case 'backward'
        indSelected = (idxFirst - 1) + (nValues + 1 - indSelected);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%