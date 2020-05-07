function vecsNew = align_vectors_by_index (vecsOrig, idxToAlign, varargin)
%% Aligns vectors by an index from each vector
% Usage: vecsNew = align_vectors_by_index (vecsOrig, idxToAlign, varargin)
% Explanation:
%       TODO
%       cf. force_matrix.m
%
% Example(s):
%       vecsNew = align_vectors_by_index(magic(4), [3, 2, 4, 1])
%
% Outputs:
%       vecsNew     - TODO: Description of vecsNew
%                   specified as a TODO
%
% Arguments:
%       vecsOrig     - TODO: Description of vecsOrig
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/array_fun.m
%       cd/create_error_for_nargin.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/set_default_flag.m
%
% Used by:
%       cd/create_plot_movie.m

% File History:
% 2020-05-06 Created by Adam Lu
% TODO: Use in force_matrix.m?

%% Hard-coded parameters

% TODO: Make this an optional argument
idxAlignedNew = [];

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
addRequired(iP, 'vecsOrig');
addRequired(iP, 'idxToAlign', ...
    @(x) validateattributes(x, {'numeric'}, {'vector', 'positive', 'integer'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, vecsOrig, idxToAlign, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Determine whether vectors was a matrix
wasMatrix = set_default_flag([], ~iscell(vecsOrig));

% Force as a column cell array of column vectors
vecsOrig = force_column_vector(vecsOrig, 'ForceCellOutput', true);

% Force as a column vector
idxToAlign = idxToAlign(:);

% Count the number of samples in each vector
nSamplesEachVec = cellfun(@numel, vecsOrig);

% Find the maximum of idxToAlign
maxIdx = max(idxToAlign);

% Find the maximum of nSamples - idxToAlign
maxNMinusIdx = max(nSamplesEachVec - idxToAlign);

% Set the new index for the aligned sample points as the maximum of idxToAlign
if isempty(idxAlignedNew)
    idxAlignedNew = maxIdx;
end

% Compute the new number of samples
nSamplesNew = idxAlignedNew + maxNMinusIdx;

% Compute the number to NaNs to pad before
nToPadPre = idxAlignedNew - idxToAlign;

% Compute the number to NaNs to pad after
nToPadPost = (nSamplesNew - idxAlignedNew) - (nSamplesEachVec - idxToAlign);

% Pad with NaNs
vecsNew = array_fun(@(a, b, c) [nan(a, 1); b; nan(c, 1)], ...
                    num2cell(nToPadPre), vecsOrig, num2cell(nToPadPost), ...
                    'UniformOutput', false);

% Convert back to numeric array
if wasMatrix
    vecsNew = force_matrix(vecsNew);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
