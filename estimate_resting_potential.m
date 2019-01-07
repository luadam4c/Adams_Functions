function [epas, Rin] = estimate_resting_potential (holdPotential, holdCurrentPa, varargin)
%% Estimates the resting membrane potential (mV) and the input resistance (MOhm) from holding potentials and holding currents
% Usage: [epas, Rin] = estimate_resting_potential (holdPotential, holdCurrentPa, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       epas        - resting membrane potential (mV)
%                   specified as a positive scalar
%       Rin         - input resistance (MOhm)
%                   specified as a positive scalar
% Arguments:
%       holdPotential   - holding potentials in mV
%                       must be a numeric vector
%       holdCurrentPa   - holding currents in pA
%                       must be a numeric vector
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/find_passive_params.m

% File History:
% 2018-11-13 Adapted from ~/m3ha/optimizer4gabab/import_rawtraces.m
% 

%% Hard-coded constants
PA_PER_NA = 1000;

%% Default values for optional arguments

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
addRequired(iP, 'holdPotential', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'holdCurrentPa', ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Read from the Input Parser
parse(iP, holdPotential, holdCurrentPa, varargin{:});

%% Preparation
% Count the number of holding potentials and holding currents
nPotentials = length(holdPotential);
nCurrents = length(holdCurrentPa);

% Make sure they are the same
if nPotentials ~= nCurrents
    fprintf(['The number of holding potentials provided do not ', ...
                'match the number of holding currents provided!']);
    Rin = [];
    epas = [];
    return
end

% Force all vectors to be column vectors
holdPotential = holdPotential(:);
holdCurrentPa = holdCurrentPa(:);

%% Do the job
% Convert the holding current from pA to nA
holdCurrentNa = holdCurrentPa / PA_PER_NA;

% Estimate epas and Rin with linear least squares
% Define the matrices
% Ohm's Law: I * R + epas = V
%     units: [nA] * [MOhm] + [mV] = [mV]
% X * w = V
X = ones(nPotentials, 2);
X(:, 1) = holdCurrentNa;
V = holdPotential;

% Compute the estimates
w = pinv(X) * V;

%% Output results
% Extract the input resistance (MOhm)
Rin = w(1);         

% Extract the resting membrane potential (mV)
epas = w(2);

%% Plot I-V curve
%{
% Check if the values make sense
if Rin <= 0
    colorRin = 'r';
else
    colorRin = 'k';
end

% Construct a vector of holding currents
holdCurrentToPlot = linspace(min(holdCurrentNa), max(holdCurrentNa), 1000);

% Compute predicted values
holdPotentialPredicted = holdCurrentToPlot * Rin + epas;

% Plot holding potential versus holding current
h = figure('Visible', 'off');
clf(h);
hold on;
plot(holdCurrentNa, holdPotential, 'o', 'LineWidth', 2);
plot(holdCurrentToPlot, holdPotentialPredicted, 'r')
text(0.1, 0.9, ['Rin = ', num2str(Rin), ' MOhm'], ...
    'Units', 'normalized', 'Color', colorRin);
text(0.1, 0.85, ['epas = ', num2str(epas), ' mV'], ...
    'Units', 'normalized', 'Color', 'k');
% text(0.1, 0.8, ['Slope = ', num2str(RinRegression), ' MOhm'], ...
%     'Units', 'normalized', 'Color', 'r');
title(['Voltage-Current relationship for ', outparams.cellname]);
ylabel('Holding potential (mV)');
xlabel('Holding current (nA)');
figName = fullfile(outparams.outFolderName, ...
                    [outparams.prefix, ...
                    '_voltage-current-relationship.png']);
saveas(h, figName);
close(h)
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the regression coefficient
RinRegression = (holdCurrent - mean(holdCurrent)) \ ...
                    (holdPotential - mean(holdPotential));

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%