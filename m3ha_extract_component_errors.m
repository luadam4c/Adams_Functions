function [componentErrors, componentLabels] = ...
                m3ha_extract_component_errors (errorTable, varargin)
%% Extracts component errors from an error table
% Usage: [componentErrors, componentLabels] = ...
%               m3ha_extract_component_errors (errorTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       componentErrors - TODO: Description of output1
%                       specified as a TODO
%       componentLabels - TODO: Description of output1
%                       specified as a TODO
%
% Arguments:
%       errorTable  - TODO: Description of errorTable
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/argfun.m
%       cd/create_error_for_nargin.m
%       cd/create_labels_from_numbers.m
%       cd/force_matrix.m
%       cd/is_field.m
%
% Used by:
%       cd/m3ha_rank_neurons.m

% File History:
% 2020-01-03 Moved from m3ha_rank_neurons.m
% 

%% Hard-coded parameters
%   Note: The following must be consistent with compute_single_neuron_errors.m
idxSweep = 1;
idxMatch = 2;
idxAmp = 3;
idxTime = 4;
idxSlope = 5;
avgSwpErrorStr = 'avgSwpError';
ltsMatchErrorStr = 'ltsMatchError';
avgLtsAmpErrorStr = 'avgLtsAmpError';
avgLtsDelayErrorStr = 'avgLtsDelayError';
avgLtsSlopeErrorStr = 'avgLtsSlopeError';
errorWeightsStr = 'errorWeights';

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

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
addRequired(iP, 'errorTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
parse(iP, errorTable, varargin{:});
% param1 = iP.Results.param1;

%% Do the job
% Extract fields
[avgSwpErrors, ltsMatchErrors, avgLtsAmpErrors, ...
        avgLtsDelayErrors, avgLtsSlopeErrors] = ...
    argfun(@(x) errorTable.(x), ...
            avgSwpErrorStr, ltsMatchErrorStr, avgLtsAmpErrorStr, ...
            avgLtsDelayErrorStr, avgLtsSlopeErrorStr);


% Extract errorWeights
if is_field(errorTable, errorWeightsStr)
    % Errors weights is a cell array of numeric vectors
    errorWeights = errorTable.(errorWeightsStr);

    % Force as a matrix with each column corresponding to a type of error
    errorWeights = transpose(force_matrix(errorWeights));
else
    % Create error weight labels for separated columns
    errorWeightsStrs = create_labels_from_numbers(1:5, ...
                                            'Prefix', [errorWeightsStr, '_']);

    % Exract each column
    errorWeights = cellfun(@(x) errorTable.(x), errorWeightsStrs, ...
                            'UniformOutput', false);

    % Force as a matrix with each column corresponding to a type of error
    errorWeights = force_matrix(errorWeights);
end

% Make sure error weights are normalized
errorWeights = errorWeights ./ repmat(sum(errorWeights, 2), 1, ...
                                        size(errorWeights, 2));

% Compute components of total error
%   Note: must match componentLabels
componentErrors = [ltsMatchErrors, avgSwpErrors, avgLtsAmpErrors, ...
                    avgLtsDelayErrors, avgLtsSlopeErrors] .* ...
        errorWeights(:, [idxMatch, idxSweep, idxAmp, idxTime, idxSlope]);

% Decide on component labels
%   Note: must match componentErrors
componentLabels = {'LTS Match Error', 'Sweep Error', 'LTS Amp Error', ...
                    'LTS Time Error', 'LTS Slope Error'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%