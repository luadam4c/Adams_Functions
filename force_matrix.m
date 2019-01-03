function vecs = force_matrix (vecs, varargin)
%% Forces vectors into a non-cell array matrix
% Usage: vecs = force_matrix (vecs, varargin)
% Explanation:
%       TODO
% Example(s):
%       force_matrix({1:5, 1:3, 1:4})
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
%                   default == 'leftAdjust'
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_subvectors.m
%
% Used by:
%       cd/compute_combined_trace.m
%       cd/find_window_endpoints.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/plot_swd_histogram.m

% File History:
% 2019-01-03 Created by Adam Lu
% 

%% Hard-coded parameters
validAlignMethods = {'leftAdjust', 'rightAdjust', ...
                    'leftAdjustPad', 'rightAdjustPad'};

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
% Do nothing if already a matrix
if ~iscell(vecs) && ismatrix(vecs)
    return
end

% Extract vectors padded on the right
vecs = extract_subvectors(vecs, 'AlignMethod', alignMethod);

% Put together as an array
vecs = horzcat(vecs{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
