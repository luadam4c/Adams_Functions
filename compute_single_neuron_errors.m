function errors = compute_single_neuron_errors (vSim, vReal, varargin)
%% Computes all errors for single neuron data
% Usage: errors = compute_single_neuron_errors (vSim, vReal, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:    
%       vSim        - simulated voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       vReal       - recorded voltage traces
%                   must be a numeric vector or a cell array of numeric vectors
%       varargin    - 'TimeVecs': TODO: Description of TimeVecs
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsSim': TODO: Description of IvecsSim
%                   must be a TODO
%                   default == TODO
%                   - 'IvecsReal': TODO: Description of IvecsReal
%                   must be a TODO
%                   default == TODO
%                   - 'SimMode': TODO: Description of SimMode
%                   must be a TODO
%                   default == TODO
%                   - 'FitWindow': TODO: Description of FitWindow
%                   must be a TODO
%                   default == TODO
%                   - 'BaseWindow': TODO: Description of BaseWindow
%                   must be a TODO
%                   default == TODO
%                   - 'BaseNoise': TODO: Description of BaseNoise
%                   must be a TODO
%                   default == TODO
%                   - 'SweepWeights': TODO: Description of SweepWeights
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       cd/find_window_endpoints.m
%       cd/force_row_numeric.m
%       cd/iscellnumeric.m
%       cd/match_row_count.m
%
% Used by:    
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

% File History:
% 2018-10-24 Adapted from code in run_neuron_once_4compgabab.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeVecsDefault = ;         % default TODO: Description of param1
ivecsSimDefault = ;         % default TODO: Description of param1
ivecsRealDefault = ;        % default TODO: Description of param1
simModeDefault = ;          % default TODO: Description of param1
fitWindowDefault = ;        % default TODO: Description of param1
baseWindowDefault = ;       % default TODO: Description of param1
baseNoiseDefault = ;        % default TODO: Description of param1
sweepWeightsDefault = ;     % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'vSim', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vSim must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'vReal', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vReal must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeVecs', timeVecsDefault, ...
    % TODO: validation function %);
addParameter(iP, 'IvecsSim', ivecsSimDefault, ...
    % TODO: validation function %);
addParameter(iP, 'IvecsReal', ivecsRealDefault, ...
    % TODO: validation function %);
addParameter(iP, 'SimMode', simModeDefault, ...
    % TODO: validation function %);
addParameter(iP, 'FitWindow', fitWindowDefault, ...
    % TODO: validation function %);
addParameter(iP, 'BaseWindow', baseWindowDefault, ...
    % TODO: validation function %);
addParameter(iP, 'BaseNoise', baseNoiseDefault, ...
    % TODO: validation function %);
addParameter(iP, 'SweepWeights', sweepWeightsDefault, ...
    % TODO: validation function %);

% Read from the Input Parser
parse(iP, vSim, vReal, varargin{:});
tBoth = iP.Results.TimeVecs;
iSim = iP.Results.IvecsSim;
iReal = iP.Results.IvecsReal;
simMode = iP.Results.SimMode;
fitWindow = iP.Results.FitWindow;
baseWindow = iP.Results.BaseWindow;
baseNoise = iP.Results.BaseNoise;
sweepWeights = iP.Results.SweepWeights;

% Make sure all windows, if a vector and not a matrix, are row vectors
if isvector(fitWindow)
    fitWindow = force_row_numeric(fitWindow);
end

%% Preparation
% Count the number of sweeps
nSweeps = numel(vSim);

% Match row counts for sweep-dependent variables with the number of sweeps
fitWindow = match_row_count(fitWindow, nSweeps);

% Simulation mode specific operations
switch simMode
    case 'passive'
    case 'active'
    otherwise
        error('simMode unrecognized!');
end

% Get the lower and upper bounds of the fit windows
fitWindowLowerBounds = fitWindow(:, 1);
fitWindowUpperBounds = fitWindow(:, 2);

% Extract the regions to fit
tBoth = cellfun();
vSim = cellfun();
vReal = cellfun();
iSim = cellfun();
iReal = cellfun();

% compute root-mean-squared error
parfor iSwp = 1:nSweeps
    % Get the time vector for this sweep
    timeVec = realData{iSwp}(:, 1);

    % Find the indices of the time vector for fitting
    indFitWin = find(timeVec >= fitWindowLowerBound(iSwp) & ...
                     timeVec <= fitWindowUpperBound(iSwp));

    % Extract the vectors for fitting
    tFit{iSwp} = timeVec(indFitWin);
    vreal{iSwp} = realData{iSwp}(indFitWin, 2);
    vsim{iSwp} = simData{iSwp}(indFitWin, 2);
    ireal{iSwp} = realData{iSwp}(indFitWin, 3);
    isim{iSwp} = simData{iSwp}(indFitWin, 9);

    % Compute root-mean-squared error (mV) over the fit window
    rmse(iSwp) = compute_rms_error(vreal{iSwp}, vsim{iSwp});
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

fitWindowLowerBounds = match_vector_counts(fitWindowLowerBounds, [nSweeps, 1]);
fitWindowUpperBounds = match_vector_counts(fitWindowUpperBounds, [nSweeps, 1]);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%