function vecs = force_matrix (vecs, varargin)
%% Forces vectors into a non-cell array matrix
% Usage: vecs = force_matrix (vecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       force_matrix({1:5, 1:3, 1:4})
%       force_matrix({1:5, 1:3, 1:4}, 'AlignMethod', 'leftadjust')
%       force_matrix({1:5, 1:3, 1:4}, 'AlignMethod', 'none')
%       force_matrix({1:5, magic(3)})
%       force_matrix({{1:3, 1:3}, {1:3, 1:3}})
%       force_matrix({1:5, 1:3, []})
%       force_matrix({{'a', 'b'}; {'b'; 'c'; 'd'}})
%       force_matrix({{'a', 'b'}; {'b'; 'c'; 'd'}}, 'AlignMethod', 'leftadjust')
%       load_examples;
%       force_matrix(myCellCellNumeric, 'TreatCellAsArray', true);
%       force_matrix({{1:5; 1:3}, {1:2; 1:6}}, 'TreatCellAsArray', true)
%       force_matrix({{1:5}, {1:2; 1:6}}, 'TreatCellNumAsArray', true)
%
% Outputs:
%       vecs        - vectors as a matrix
%                   specified as a matrix
%
% Arguments:
%       vecs        - vectors
%                   must be an array
%       varargin    - 'AlignMethod': method for truncation or padding
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'leftAdjust'  - align to the left and truncate
%                       'rightAdjust' - align to the right and truncate
%                       'leftAdjustPad'  - align to the left and pad
%                       'rightAdjustPad' - align to the right and pad
%                       'none'        - no alignment/truncation
%                   default == 'leftAdjustPad'
%                   - 'TreatCellAsArray': whether to treat a cell array
%                                           as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellNumAsArray': whether to treat a cell array
%                                       of numeric arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'TreatCellStrAsArray': whether to treat a cell array
%                                       of character arrays as a single array
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'Verbose' - whether to print to standard output
%                                   regardless of message mode
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - Any other parameter-value pair for force_column_vector()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%       cd/force_column_vector.m
%       cd/iscellnumeric.m
%       cd/isnum.m
%
% Used by:
%       cd/alternate_elements.m
%       cd/combine_abf_data.m
%       cd/compute_combined_data.m
%       cd/compute_combined_trace.m
%       cd/compute_relative_event_times.m
%       cd/create_indices.m
%       cd/create_synced_movie_trace_plot_movie.m
%       cd/detect_spikes_current_clamp.m
%       cd/find_window_endpoints.m
%       cd/force_column_vector.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/parse_multiunit.m
%       cd/plot_chevron.m
%       cd/plot_measures.m
%       cd/plot_swd_histogram.m
%       cd/plot_traces_spike2_mat.m
%       cd/plot_tuning_curve.m
%       cd/spike2Mat2Text.m
%       cd/transform_vectors.m
%       cd/vecfun.m

% File History:
% 2019-01-03 Created by Adam Lu
% 2019-01-22 Added a quick return for performance
% 2019-04-26 Fixed bug for 'AlignMethod' == 'none'
% 2019-09-07 Added 'Verbose' as an optional argument
% 2019-10-02 Now returns if already a matrix
% 2019-10-03 Added 'TreatCellNumAsArray' as an optional argument
% TODO: Restrict the number of samples if provided
% 

%% Quick return for performance
% Do nothing if already a matrix
if ~iscell(vecs) && ismatrix(vecs)
    return
end

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad', 'none'};

%% Default values for optional arguments
alignMethodDefault  = 'leftAdjustPad';   % pad on the right by default
treatCellAsArrayDefault = false;% treat cell arrays as many arrays by default
treatCellNumAsArrayDefault = false; % treat cell arrays of numeric arrays
                                    %   as many arrays by default
treatCellStrAsArrayDefault = true;  % treat cell arrays of character arrays
                                    %   as an array by default
verboseDefault = true;              % print warning by default

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
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));
addParameter(iP, 'TreatCellAsArray', treatCellAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellNumAsArray', treatCellNumAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'TreatCellStrAsArray', treatCellStrAsArrayDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs, varargin{:});
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);
treatCellAsArray = iP.Results.TreatCellAsArray;
treatCellNumAsArray = iP.Results.TreatCellNumAsArray;
treatCellStrAsArray = iP.Results.TreatCellStrAsArray;
verbose = iP.Results.Verbose;

% Keep unmatched arguments for the force_column_vector() function
otherArguments = iP.Unmatched;

%% Do the job
% Return if already a matrix
if isnum(vecs) || ~isvector(vecs) && ...
        treatCellAsArray && iscell(vecs) || ...
        treatCellNumAsArray && iscellnumeric(vecs) || ...
        treatCellStrAsArray && iscellstr(vecs)
    return
end

% Extract vectors padded on the right
%   Note: don't do this if alignMethod is set to 'none'
%           Otherwise, extract_subvectors.m uses create_indices.m,
%         	which uses force_matrix.m and will enter infinite loop
if ~strcmpi(alignMethod, 'none')
    vecs = extract_subvectors(vecs, 'AlignMethod', alignMethod, ...
                            'TreatCellAsArray', treatCellAsArray, ...
                            'TreatCellNumAsArray', treatCellNumAsArray, ...
                            'TreatCellStrAsArray', treatCellStrAsArray);
else
    vecs = force_column_vector(vecs, 'TreatCellAsArray', treatCellAsArray, ...
                    'TreatCellNumAsArray', treatCellNumAsArray, ...
                    'TreatCellStrAsArray', treatCellStrAsArray, otherArguments);
end

% Count the number of samples for each vector
nSamples = cellfun(@numel, vecs);

% Find the maximum and minimum number of samples
maxNSamples = max(nSamples);
minNSamples = min(nSamples);

% Put together as an array
if maxNSamples == minNSamples
    if maxNSamples == 0
        % Don't force as matrix
        if verbose
            disp(['Warning: Vector lengths are all zero, ', ...
                    'concatenation aborted!']);
        end
    else
        vecs = horzcat(vecs{:});
    end
else
    if verbose
        disp(['Warning: Vector lengths are not consistent, ', ...
                'concatenation aborted!']);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

nUniqueNSamples = numel(unique(cellfun(@numel, vecs)));
if nUniqueNSamples == 1

% Count the number of vectors
nVectors = numel(vecs);

% Create an empty matrix that matches the number of vectors
vecs = double.empty(0, nVectors);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
