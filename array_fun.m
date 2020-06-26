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
%                   - 'RenewParpool': whether to renew parallel parpool 
%                                       every batch to release memory
%                   must be logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/cleanup_parcluster.m
%       cd/create_error_for_nargin.m
%       cd/extract_parameter_value_pairs.m
%       cd/force_column_cell.m
%       cd/is_in_parallel.m
%       cd/rmfield_custom.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/align_vectors_by_index.m
%       cd/argfun.m
%       cd/extract_vars.m
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/check_fullpath.m
%       cd/combine_param_tables.m
%       cd/compute_combined_trace.m
%       cd/compute_confidence_ellipse.m
%       cd/compute_rms_error.m
%       cd/compute_weighted_average.m
%       cd/count_samples.m
%       cd/decide_on_colormap.m
%       cd/extract_columns.m
%       cd/extract_elements.m
%       cd/extract_fields.m
%       cd/extract_subvectors.m
%       cd/find_closest.m
%       cd/find_in_list.m
%       cd/find_window_endpoints.m
%       cd/force_string_end.m
%       cd/ismatch.m
%       cd/ismember_custom.m
%       cd/read_neuron_outputs.m
%       cd/m3ha_compute_statistics.m
%       cd/m3ha_import_raw_traces.m
%       cd/m3ha_plot_figure03.m
%       cd/m3ha_plot_figure07.m
%       cd/parse_atf_swd.m
%       cd/sscanf_full.m
%       cd/vecfun.m

% File History:
% 2020-01-01 Created by Adam Lu
% 2020-01-02 Fixed to work with 2D arrays
% 2020-03-09 Added RenewParpool as an optional argument
% 2020-04-26 Fixed bug when nArgOut is 0
% 2020-05-16 Fixed concatenation using parfor when uniform output
% TODO: Convert all arguments to a cell array (with num2cell) 
%       if any argument is a cell array

%% Hard-coded parameters
minItemsForParfor = 12;
maxNumWorkers = 12;

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

%% Preparation
% Count the number of output arguments
nArgOut = nargout;

% Remove any parameter-value pairs
[params, inputList] = extract_parameter_value_pairs(varargin);

% Deal with parameter-value pairs for this function
% TODO: Use parse_and_remove_from_struct.m
if isfield(params, 'RenewParpool')
    renewParpool = params.RenewParpool;
    params = rmfield_custom(params, 'RenewParpool');
else
    renewParpool = false;
end

% Count the number of items
nItems = numel(inputList{1});

% Save old dimensions
oldDimensions = size(inputList{1});

%% Do the job
if is_in_parallel || nItems < minItemsForParfor
    % Place parameters in a list
    paramList = struct2arglist(params);

    %% Use cellfun or arrayfun
    if iscell(varargin{1})
        [varargout{1:nargout}] = cellfun(myFunc, inputList{:}, paramList{:});
    else
        [varargout{1:nargout}] = arrayfun(myFunc, inputList{:}, paramList{:});
    end
else
    %% Use parfor

    % Force all arguments as column cell vectors of dimension nItems x 1
    inputColumn = cellfun(@force_column_cell_individually, inputList, ...
                            'UniformOutput', false);

    % Concatenate horizontally to make a cell matrix (nItems x nargin)
    inputMatrix = horzcat(inputColumn{:});

    % Initialize a cell matrix for all outputs requested (nItems x nArgOut)
    outputMatrix = cell(nItems, nArgOut);

    % Cleanup previous paralle pools and start a new one
    if renewParpool
        cleanup_parcluster;

        % Get current parallel pool object or create a new one
        poolObj = gcp;

        % number of workers in the current or default parallel pool object
        oldNumWorkers = poolObj.NumWorkers;

        % Maximum number of workers to use
        numWorkers = min(oldNumWorkers, maxNumWorkers); 
    end

    % Initialize count of number of items completed
    ct = 0;

    % Run through all items
    while ct < nItems
        % Locate the first item in this batch
        first = ct + 1;

        % Locate the last item in this batch
        if renewParpool && ct + numWorkers <= nItems
            % If memory is to be released, limit the batch to numWorkers
            last = ct + numWorkers;
        else
            last = nItems;
        end

        parfor iItem = first:last
            % Get all the arguments for this item
            inputsThis = inputMatrix(iItem, :);

            % Apply myFunc to these arguments
            outputsThis = apply_func(myFunc, inputsThis, nArgOut);

            % Save in output
            if ~isempty(outputsThis)
                outputMatrix(iItem, :) = outputsThis;
            end
        end

        % Renew parallel pool object to clear memory
        if renewParpool
            % Delete the parallel pool object to release memory
            delete(poolObj);

            % Recreate a parallel pool object with desired number of workers
            poolObj = parpool('local', numWorkers);
        end

        % Update the number of trials completed
        ct = last;
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
        varargout = cellfun(@cell2num_custom, outputCells, 'UniformOutput', false);
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

if nArgOut >= 1
    [outputList{1:nArgOut}] = myFunc(argList{:});
else
    myFunc(argList{:});
    outputList = {};
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function array = cell2num_custom (cellArray)
%% Converts a cell array to a non-cell array by concatenation
% TODO: Pull out as its own function

% List all contents and concatenate
array = [cellArray{:}];

% Reshape as original dimensions
array = reshape(array, size(cellArray));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of arguments
nArgs = numel(varargin);
% Create variables for each argument
arrayfun(@(x) eval('arg%d = varargin{%d}'), 1:nArgs, 1:nArgs);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
