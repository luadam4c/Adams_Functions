function varargout = array_fun (myFunc, varargin)
%% Applies cellfun or arrayfun based on the input type, or use parfor if not already in a parallel loop
% Usage: varargout = array_fun (myFunc, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [a, b] = array_fun(@(x, y) min([x, y]), [3; 1], [4; 5])
%       [a, b] = array_fun(@(x, y) min([x, y]), [3, 1], [4, 5])
%       [a, b] = array_fun(@(x, y) min([x, y]), {3; 1}, {4; 5})
%       [a, b] = array_fun(@(x, y) min([x, y]), {3, 1}, {4, 5})
%       [a, b] = array_fun(@(x, y) min([x, y]), {3, 1; 2, 4}, {4, 5; 1, 2})
%
% Outputs:
%       varargout   - outputs of cellfun or arrayfun
%
% Arguments:
%       myFunc      - function to apply over cells or arrays
%                   must be a function handle
%       varargin    - 2nd to last arguments to cellfun or arrayfun
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/force_column_cell.m
%       cd/is_in_parallel.m
%
% Used by:
%       cd/argfun.m
%       cd/extract_parameter_value_pairs.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/check_fullpath.m
%       cd/compute_combined_trace.m
%       cd/compute_rms_error.m
%       cd/compute_weighted_average.m
%       cd/count_samples.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/extract_fields.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_in_list.m
%       cd/find_window_endpoints.m
%       cd/force_string_end.m
%       cd/ismember_custom.m
%       cd/load_neuron_outputs.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_network_analyze_spikes.m
%       cd/vecfun.m

% File History:
% 2020-01-01 Created by Adam Lu
% 2020-01-02 Fixed to work with 2D arrays
% TODO: Renew parpool if memory usage is too high
% TODO: Convert all arguments to a cell array (with num2cell) 
%       if any argument is a cell array

%% Hard-coded parameters
minItemsForParfor = 12;

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
addRequired(iP, 'myFunc', ...
    @(x) validateattributes(x, {'function_handle'}, {'scalar'}));

% Read from the Input Parser
parse(iP, myFunc);

%% Do the job
if is_in_parallel || numel(varargin{1}) < minItemsForParfor
    % Use cellfun or arrayfun
    if iscell(varargin{1})
        [varargout{1:nargout}] = cellfun(myFunc, varargin{:});
    else
        [varargout{1:nargout}] = arrayfun(myFunc, varargin{:});
    end
else
    % Use parfor
    % Count the number of output arguments
    nArgOut = nargout;

    % Count the number of items
    nItems = numel(varargin{1});

    % Save old dimensions
    oldDimensions = size(varargin{1});

    % Remove any parameter-value pairs
    [params, inputList] = extract_parameter_value_pairs(varargin);

    % Force all arguments as column cell vectors of dimension nItems x 1
    inputColumn = cellfun(@force_column_cell_individually, inputList, ...
                            'UniformOutput', false);

    % Concatenate horizontally to make a cell matrix (nItems x nargin)
    inputMatrix = horzcat(inputColumn{:});

    % Initialize a cell matrix for all outputs requested (nItems x nArgOut)
    outputMatrix = cell(nItems, nArgOut);

    % Run through all items
    parfor iItem = 1:nItems
        % Get all the arguments for this item
        inputsThis = inputMatrix(iItem, :);

        % Apply myFunc to these arguments
        outputsThis = apply_func(myFunc, inputsThis, nArgOut);

        % Save in output
        outputMatrix(iItem, :) = outputsThis;
    end

    % Reorganize outputs
    outputCells = arrayfun(@(x) outputMatrix(:, x), 1:nArgOut, ...
                            'UniformOutput', false);

    % TODO: Deal with other parameter-value pairs in cellfun or arrayfun

    % Apply old dimensions
    if oldDimensions(2) ~= 1 || ...
            numel(oldDimensions) > 2 && oldDimensions(3) ~= 1
        outputCells = cellfun(@(x) reshape(x, oldDimensions), ...
                                    outputCells, 'UniformOutput', false);
    end

    % Concatenate outputs unless requested not to
    if isfield(params, 'UniformOutput') && ~params.UniformOutput || ...
            isfield(params, 'uniformOutput') && ~params.uniformOutput
        varargout = outputCells;
    else
        % Try to concatenate outputs as non-cell arrays
        varargout = cellfun(@cell2num, outputCells, 'UniformOutput', false);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function array = force_column_cell_individually (array)
%% Forces an array as a column cell vector
%   Note: this behaves differently from force_column_cell.m

% Force as a cell array using num2cell()
if ~iscell(array)
    array = num2cell(array);
end

% Force as a column vector
array = force_column_cell(array, 'ToLinearize', true);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function outputList = apply_func(myFunc, argList, nArgOut)
%% Applies a function to inputs from a cell array and returns outputs in a cell array

[outputList{1:nArgOut}] = myFunc(argList{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of arguments
nArgs = numel(varargin);
% Create variables for each argument
arrayfun(@(x) eval('arg%d = varargin{%d}'), 1:nArgs, 1:nArgs);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%