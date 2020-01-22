function gGababNs = compute_gabab_conductance (tVec, timeStart, amplitudeNs, ...
                            tauRise, tauFallFast, tauFallSlow, weight, varargin)
%% Computes a the conductance over time for a GABAB-IPSC
% Usage: gGababNs = compute_gabab_conductance (tVec, timeStart, amplitudeNs, ...
%                           tauRise, tauFallFast, tauFallSlow, weight, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       gGababNs    - GABAB-IPSC conductance (nS)
%                   specified as a numeric vector
%
% Arguments:
%       tVec        - time vector (ms)
%                   must be a numeric vector
%       timeStart   - time of IPSC start (ms)
%                   must be a numeric scalar
%       amplitudeNs - amplitude (nS, but saved in spreadsheet as uS)
%                   must be a numeric scalar
%       tauRise     - time constant for rising phase (ms)
%                   must be a numeric scalar
%       tauFallFast - time constant for fast falling phase (ms)
%                   must be a numeric scalar
%       tauFallSlow - time constant for slow falling phase (ms)
%                   must be a numeric scalar
%       weight      - weight of fast falling phase in the falling phase
%                   must be a numeric scalar
%       varargin    - 'SheetName': spreadsheet name for 
%                               preparing for m3ha_neuron_run_and_analyze.m
%                   must be a string scalar or a character vector
%                   default == '' (no output)
%                   - Any other parameter-value pair for TODO()
%
% Requires:
%       cd/find_closest.m
%       cd/force_column_vector.m
%       cd/force_matrix.m
%       cd/force_row_vector.m
%       cd/match_format_vectors.m
%
% Used by:
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_network_compare_ipsc.m
%       cd/m3ha_resave_sweeps.m
%       cd/m3ha_trace_comparison.m

% File History:
% 2016-11-07 Moved from /media/adamX/m3ha/data_dclamp/trace_comparison.m
% 2020-01-03 Now accepts multiple time vectors and/or parameters
%               at the same time
% 2020-01-04 Added 'SaveParams' as an optional argument

%% Hard-coded constants
NS_PER_US = 1000;

%% Hard-coded parameters
% Note: Must be consistent with m3ha_neuron_run_and_analyze.m
sheetParamNames = {'IpscTimeOrig', 'GababAmpIpscr', 'GababTriseIpscr', ...
        'GababTfallFastIpscr', 'GababTfallSlowIpscr', 'GababWeightIpscr'};

%% Default values for optional arguments
sheetNameDefault = '';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 7
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
% addRequired(iP, 'reqarg1', ...                  % TODO: Description of reqarg1
%     % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SheetName', sheetNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
% parse(iP, reqarg1, varargin{:});
parse(iP, varargin{:});
sheetName = iP.Results.SheetName;

%% Preparation
% Make sure time vectors are column vectors
tVec = force_column_vector(tVec);

% Force as a matrix
tVec = force_matrix(tVec);

% Make sure parameters are row vectors
[timeStart, amplitudeNs, tauRise, tauFallFast, tauFallSlow, weight] = ...
    argfun(@force_row_vector, ...
            timeStart, amplitudeNs, tauRise, tauFallFast, tauFallSlow, weight);

% Compute the shifted time vectors
timeShifted = tVec - timeStart;

% Compute the exponential components
[expRise, expFallFast, expFallSlow] = ...
    argfun(@(x) exp(-timeShifted ./ x), tauRise, tauFallFast, tauFallSlow);

% Compute the GABA-B conductance
gGababNs = amplitudeNs .* (1 - expRise) .^ 8 .* ...
                (expFallFast .* weight + expFallSlow .* (1 - weight));  

% Compute the index of the starting time
idxTimeStart = find_closest(tVec, timeStart);

% For each column, set everything before timeStart to be zero
if numel(idxTimeStart) > 1
    for iColumn = 1:size(gGababNs, 2)
        gGababNs(1:idxTimeStart(iColumn), iColumn) = 0;
    end
else
    gGababNs(1:idxTimeStart, :) = 0;
end

% Back up parameters in a spreadsheet if requested
if ~isempty(sheetName)
    % Convert amplitude for nS to uS
    amplitudeUs = amplitudeNs ./ NS_PER_US;

    % Make sure all are column vectors with the same number of rows
    [timeStart, amplitudeUs, tauRise, tauFallFast, tauFallSlow, weight] = ...
        match_format_vectors(timeStart, amplitudeUs, tauRise, ...
                    tauFallFast, tauFallSlow, weight, 'RowInstead', false);

    % Create the table
    paramTable = table(timeStart, amplitudeUs, tauRise, tauFallFast, ...
                    tauFallSlow, weight, 'VariableNames', sheetParamNames);

    % Save the table
    writetable(paramTable, sheetName);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
