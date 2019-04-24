function vecs = force_matrix (vecs, varargin)
%% Forces vectors into a non-cell array matrix
% Usage: vecs = force_matrix (vecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       force_matrix({1:5, 1:3, 1:4})
%       force_matrix({1:5, magic(3)})
% Outputs:
%       vecs        - vectors as a matrix
%                   specified as a matrix
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
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%
% Used by:
%       cd/compute_combined_data.m
%       cd/compute_combined_trace.m
%       cd/create_indices.m
%       cd/find_window_endpoints.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_swd_histogram.m

% File History:
% 2019-01-03 Created by Adam Lu
% 2019-01-22 Added a quick return for performance
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AlignMethod', alignMethodDefault, ...
    @(x) any(validatestring(x, validAlignMethods)));

% Read from the Input Parser
parse(iP, vecs, varargin{:});
alignMethod = validatestring(iP.Results.AlignMethod, validAlignMethods);

%% Do the job
% Extract vectors padded on the right
%   Note: don't do this if alignMethod is set to none
%           Otherwise, extract_subvectors.m uses create_indices.m,
%         	which uses force_matrix.m and will enter infinite loop
if ~strcmpi(alignMethod, 'none')
    vecs = extract_subvectors(vecs, 'AlignMethod', alignMethod);
end

% Put together as an array
vecs = horzcat(vecs{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
