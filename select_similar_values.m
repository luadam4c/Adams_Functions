function [valSelected, indSelected] = select_similar_values (values, varargin)
%% Selects values that are within a certain percentage range of the mean
% Usage: [valSelected, indSelected] = select_similar_values (values, varargin)
% Explanation:
%       TODO
% Example(s):
%       values = [10, NaN, 100, 12, 8, 90, 11, 9, 10]
%       [v, i] = select_similar_values(values)
%       [v, i] = select_similar_values(values, 'SelectionMethod', 'notNaN')
%       [v, i] = select_similar_values(values, 'NToSelect', 3)
%       [v, i] = select_similar_values(values, 'Direction', 'backward')
%       [v, i] = select_similar_values(values, 'NToSelect', 3, 'MaxRange2Mean', 20)
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
%                   - 'Indices': indices for the subvectors to extract 
%                       Note: if provided, would override 'EndPoints'
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == set in extract_subvectors.m
%                   - 'EndPoints': endpoints for the subvectors to extract 
%                   must be a numeric vector with 2 elements
%                       or a numeric array with 2 rows
%                       or a cell array of numeric vectors with 2 elements
%                   default == find_window_endpoints([], values)
%                   - 'SelectionMethod': the selection method
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'notNaN'        - select any non-NaN value
%                       'maxRange2Mean' - select vales so that the maximum 
%                                           range is within a percentage 
%                                           of the mean
%                   default == 'maxRange2Mean'
%                   - 'Direction': the selection direction
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'forward'   - select from the first indices
%                       'backward'  - select from the last indices
%                   default == 'forward'
%                   - 'MaxRange2Mean': maximum percentage of range versus mean
%                   must be a nonnegative scalar
%                   default == 40%
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/extract_elements.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_phase_average.m

% File History:
% 2019-05-12 Created by Adam Lu
% 2019-05-15 Add 'SelectionMethod' as an optional argument
%               with possible values 'notNaN', 'maxRange2Mean'
% 

%% Hard-coded parameters
validSelectionMethods = {'notNaN', 'maxRange2Mean'};
validDirections = {'forward', 'backward'};

%% Default values for optional arguments
nToSelectDefault = 5;           % select 5 values by default
indicesDefault = [];            % set in extract_subvectors.m
endPointsDefault = [];          % set later
selectionMethodDefault = 'maxRange2Mean';   
                                % select using maxRange2Mean by default
directionDefault = 'forward';   % select from the first indices by default
maxRange2MeanDefault = 40;      % range is not more than 40% of mean by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'values', ...
    @(x) validateattributes(x, {'numeric', 'cell'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NToSelect', nToSelectDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer', 'scalar'}));
addParameter(iP, 'Indices', indicesDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['Indices must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'EndPoints', endPointsDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['EndPoints must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'SelectionMethod', selectionMethodDefault, ...
    @(x) any(validatestring(x, validSelectionMethods)));
addParameter(iP, 'Direction', directionDefault, ...
    @(x) any(validatestring(x, validDirections)));
addParameter(iP, 'MaxRange2Mean', maxRange2MeanDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'scalar'}));

% Read from the Input Parser
parse(iP, values, varargin{:});
nToSelect = iP.Results.NToSelect;
indices = iP.Results.Indices;
endPoints = iP.Results.EndPoints;
selectionMethod = validatestring(iP.Results.SelectionMethod, ...
                                    validSelectionMethods);
direction = validatestring(iP.Results.Direction, validDirections);
maxRange2Mean = iP.Results.MaxRange2Mean;

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Restrict to end points
valuesRestricted = extract_subvectors(values, 'EndPoints', endPoints, ...
                                        'Indices', indices);
% Save first index
if ~isempty(indices)
    idxFirst = indices(1);
elseif ~isempty(endPoints)
    endPoints = force_column_vector(endPoints);
    idxFirst = extract_elements(endPoints, 'first');
else
    idxFirst = 1;
end

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

% Count the number of values
nValues = count_samples(valuesOfInterest);

% Check if there are enough values to select from
if nToSelect > nValues
    valSelected = nan(nToSelect, 1);
    indSelected = nan(nToSelect, 1);
    return
end

%% Do the job
% Select the initial set of values
isSelectedOfInterest = false(size(valuesOfInterest));
isSelectedOfInterest(1:nToSelect) = true;

% Get the selected values and the corresponding indices in valuesOfInterest
valSelected = valuesOfInterest(isSelectedOfInterest);
indSelectedOfInterest = find(isSelectedOfInterest);

% Store the next index to use
idxNext = nToSelect + 1;

% Determine whether valid values were selected
success = test_success(selectionMethod, valuesOfInterest, ...
                        isSelectedOfInterest, maxRange2Mean);

% Perform similarity test 
while ~success && idxNext <= nValues
    %   TODO: Make a function find_most_extreme.m
    % Find an NaN if any
    iMostExtreme = find(isnan(valSelected), 1);

    % Find the most extreme index based on the distance to the mean
    if isempty(iMostExtreme)
        [~, iMostExtreme] = max(abs(valSelected - mean(valSelected)));
    end

    % Take the most extreme index out
    isSelectedOfInterest(indSelectedOfInterest(iMostExtreme)) = false;

    % Add the next index
    isSelectedOfInterest(idxNext) = true;

    % Update the next index to use
    idxNext = idxNext + 1;

    % Update the selected values
    valSelected = valuesOfInterest(isSelectedOfInterest);
    indSelectedOfInterest = find(isSelectedOfInterest);

    % Update whether valid values were selected
    success = test_success(selectionMethod, valuesOfInterest, ...
                            isSelectedOfInterest, maxRange2Mean);
end

% If not found, return NaNs
if ~success
    valSelected = nan(nToSelect, 1);
    indSelected = nan(nToSelect, 1);
    return    
end

% Reconstruct original indices
switch direction
    case 'forward'
        indSelected = (idxFirst - 1) + indSelectedOfInterest;
    case 'backward'
        indSelected = (idxFirst - 1) + ...
                        flip(nValues + 1 - indSelectedOfInterest);
end

% Output original values
valSelected = values(indSelected);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function success = test_success (selectionMethod, valuesOfInterest, ...
                                isSelectedOfInterest, maxRange2Mean)
% Test whether the selected values satify the selection criteria

switch selectionMethod
    case 'maxRange2Mean'
        % Compute the initial percentage range relative to the mean
        range2mean = compute_range2mean(valuesOfInterest, isSelectedOfInterest);

        if range2mean > maxRange2Mean
            success = false;
        else
            success = true;
        end
    case 'notNaN'
        success = ~any(isnan(valuesOfInterest(isSelectedOfInterest)));
    otherwise
        error('selectionMethod unrecognized!')
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function range2mean = compute_range2mean (values, isSelected)

% Get the corresponding set of values
valSelected = values(isSelected);

% Compute the range to mean percentage ratio
range2mean = compute_stats(valSelected, 'range2mean');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

meanSelected = mean(valSelected);
rangeSelected = range(valSelected);
range2mean = (rangeSelected / meanSelected) * 100;

indSelected = transpose(1:nToSelect);
range2mean = compute_range2mean(valuesOfInterest, indSelected);
% Take the most extreme index out
indSelected(iMostExtreme) = [];
% Add the next index
indSelected = [indSelected; idxNext];
valSelected = valuesOfInterest(isSelected);

% Find default end points if not provided
if isempty(endPoints)
    endPoints = find_window_endpoints([], values);
else
    endPoints = force_column_vector(endPoints);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%