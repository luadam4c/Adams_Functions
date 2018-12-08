function [params, algorithmInfo] = ...
                estimate_passive_params (fitCoeffs, pulseAmplitude, varargin)
%% Estimates passive parameters from fitted coefficients, current pulse amplitude and some constants
% Usage: [params, algorithmInfo] = ...
%               estimate_passive_params (fitCoeffs, pulseAmplitude, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       params      - estimated parameters inferred from fits, 
%                       estimated from the extrapolated coefficients of a long 
%                       pulse response, with fields:
%                       C0      - amplitude of slow exponential (mV)
%                       tau0    - time constant of slow exponential (ms)
%                       C1      - amplitude of fast exponential (mV)
%                       tau1    - time constant of fast exponential (ms)
%                       pulseAmplitude    
%                               - current pulse amplitude (pA)
%                       responseAmplitude
%                               - response amplitude at steady state (mV)
%                       Rinput  - input resistance (MOhm)
%                       alpha1  - alpha1 = sqrt(tau0/tau1 - 1)
%                       LInit   - initial guess for electrotonic length
%                       L       - electrotonic length
%                       rho     - dendritic to somatic conductance ratio
%                       Rs      - series (pipette) resistance (MOhm)
%                       Rmemb   - input resistance of cell (MOhm)
%                       Rsoma   - somatic resistance (MOhm)
%                       Rdend   - dendritic resistance (MOhm)
%                       Gmemb   - input conductance of cell (uS)
%                       Gsoma   - somatic conductance (uS)
%                       Gdend   - dendritic conductance (uS)
%                       taum    - membrane time constant (ms)
%                       Cm      - specific membrane capacitance (uF/cm^2)
%                       Rm      - specific membrane resistivity (Ohm-cm^2)
%                       radiusSoma       - radius of soma (um)
%                       Ra               - axial resistivity (Ohm-cm)
%                       diameterDendrite - diameter of dendrite (um)
%                       radiusDendrite   - radius of dendrite (um)
%                       lambda           - space constant (um)
%                       lengthDendrite   - length of dendrite (um)
%                   specified as a scalar structure
% Arguments:    
%       fitCoeffs   - fitted coefficients
%                   must be a structure including the fields:
%                       tauSlow (in ms)
%                       tauFast (in ms)
%                       ampSlowLpr (in mV)
%                       ampFastLpr (in mV)
%       pulseAmplitude  - current pulse amplitude (pA)
%                   must be a numeric scalar
%       varargin    - 'MembraneCapacitance': specific membrane capacitance 
%                                                   (uF/cm^2)
%                   must be a nonnegative scalar
%                   default == 0.88 uF/cm^2
%                   - 'AxialResistivity': axial resistivity (Ohm-cm)
%                   must be a nonnegative scalar
%                   default == 173 Ohm-cm
%                   - 'SeriesResistance': series (pipette) resistance (MOhm)
%                   must be a nonnegative scalar
%                   default == 10 MOhm
%
% Used by:    
%       cd/fit_and_estimate_passive_params.m

% File History:
% 2018-10-11 Moved code from find_passive_params.m
% 2018-12-07 Changed the ranges to look for (version 15-2)
% 

%% Hard-coded constants
OHM_PER_MOHM = 1e6;
V_PER_MV = 1e-3;
A_PER_PA = 1e-12;
S_PER_MS = 1e-3;
F_PER_UF = 1e-6;
UM_PER_CM = 1e4;

%% Hard-coded parameters
minLs = [1; 0.5; 1.5; 0; 2; 2.5; 3; 3.5; 4; 4.5; 5];
maxLs = [1.5; 1; 2; 0.5; 2.5; 3; 3.5; 4; 4.5; 5; Inf];
minRho = 0;

%% Default values for optional arguments
membraneCapacitanceDefault = 0.88;  % specific membrane capacitance (uF/cm^2)
axialResistivityDefault = 173;      % axial resistivity (Ohm-cm)
seriesResistanceDefault = 10;       % series (pipette) resistance (MOhm)

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
addRequired(iP, 'fitCoeffs', ...
    @(x) validateattributes(x, {'struct'}, {'2d'}));
addRequired(iP, 'pulseAmplitude', ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'MembraneCapacitance', membraneCapacitanceDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'AxialResistivity', axialResistivityDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'SeriesResistance', seriesResistanceDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, fitCoeffs, pulseAmplitude, varargin{:});
Cm = iP.Results.MembraneCapacitance;
Ra = iP.Results.AxialResistivity;
Rs = iP.Results.SeriesResistance;

%% Preparation
% Check if fitCoeffs has all the required fields
if ~isfield(fitCoeffs, 'tauSlow') || ~isfield(fitCoeffs, 'tauFast') || ...
    ~isfield(fitCoeffs, 'ampSlowLpr') || ~isfield(fitCoeffs, 'ampFastLpr')
    fprintf(['First argument must have these 4 fields: ', ...
                'tauSlow, tauFast, ampSlowLpr, ampFastLpr!!\n']);
    params = [];
    return
end

% Extract parameters and rename to be consistent with Johnston & Wu, pp. 94~96
tau0 = fitCoeffs.tauSlow;
tau1 = fitCoeffs.tauFast;
amp0 = fitCoeffs.ampSlowLpr;
amp1 = fitCoeffs.ampFastLpr;

%% Find Rinput based on amp0, amp1, pulseAmplitude
% Compute the change in membrane potential at steady state (mV)
responseAmplitude = amp0 + amp1;

% Compute input resistance (MOhm)
Rinput = (responseAmplitude * V_PER_MV) / (pulseAmplitude * A_PER_PA) / ...
            OHM_PER_MOHM;

%% Convert amp0, amp1 to C0, C1 (always positive)
C0 = abs(amp0);
C1 = abs(amp1);

%% Compute alpha1, L, rho based on C0, C1, tau0, tau1
% Compute alpha1
%   Note: This should be real because tau0 should be the slower time constant
alpha1 = sqrt(tau0/tau1 - 1);
if ~isreal(alpha1)
    error('code is wrong!');
end
fprintf('alpha1 == %g\n', alpha1);

% Compute the electrotonic length (L) 
%   and the dendritic-to-somatic conductance ratio (rho)
[L, rho, LInit, h, finalMaxL, f, x, maxLOfInterest] = ...
    compute_L_and_rho(C0, C1, tau0, tau1, alpha1, minLs, maxLs, minRho);

%% Compute geometric parameters based on Rinput, rho, L, Cm, Ra, Rs
% Compute the input resistance of the cell membrane (MOhm)
Rmemb = Rinput - Rs;

% Compute the input conductance of the cell membrane (uS)
Gmemb = 1/Rmemb;

% Compute the input conductance of the soma (uS)
Gsoma = Gmemb/(1 + rho);

% Compute the input conductance of the dendrite (uS)
Gdend = Gsoma * rho;

% Compute the input resistance of the soma and dendrite (MOhm)
Rsoma = 1/Gsoma;
Rdend = 1/Gdend;

% Estimate the membrane time constant (ms) by using the time constant of 
%   the slow component
taum = tau0;

% Compute the specific membrane resistivity (Ohm-cm^2)
%   Note: 1 F = 1 s / 1 Ohm
Rm = (taum * S_PER_MS) / (Cm * F_PER_UF);

% Assume the soma is a sphere or a cylinder with length == diameter 
%   and no ends), find the radius of the soma (um): 
%       Rsoma = Rm/(4*pi*)
%       => r = sqrt(Rm/(4*pi*Rsoma))
radiusSoma = UM_PER_CM * sqrt(Rm / (4 * pi * (Rsoma * OHM_PER_MOHM)));

% Modeling the dendrite as a cylinder, compute the diameter of the dendrite (um)
%       Rdend = (2/pi) * (Rm*Ra)^(1/2) * d^(-3/2) * coth(L)
%       => d = ((2/pi) * (Rm*Ra)^(1/2) * Rdend^(-1) * coth(L) )^(2/3)
if isinf(Rdend)
    diameterDendrite = 0;
else
    diameterDendrite = ...
        UM_PER_CM * ( (2/pi) * (Rm*Ra)^(1/2) * ...
                      (Rdend * OHM_PER_MOHM)^(-1) * coth(L) )^(2/3);
end

% Compute the radius of the dendrite (um)
radiusDendrite = diameterDendrite / 2;

% Compute the space constant lambda (um)
lambda = sqrt((radiusDendrite * Rm) / (2 * Ra));

% Compute the length of the dendrite (um)
lengthDendrite = L * lambda;

%% Store outputs
% Output algorithm info
algorithmInfo.minLs = minLs;
algorithmInfo.maxLs = maxLs;
algorithmInfo.minRho = minRho;
algorithmInfo.f = f;
algorithmInfo.x = x;
algorithmInfo.maxLOfInterest = maxLOfInterest;
algorithmInfo.finalMaxL = finalMaxL;

% Output parameters into a params structure
params.C0 = C0;
params.tau0 = tau0;
params.C1 = C1;
params.tau1 = tau1;
params.pulseAmplitude = pulseAmplitude;
params.responseAmplitude = responseAmplitude;
params.Rinput = Rinput;
params.alpha1 = alpha1;
params.LInit = LInit;
params.L = L;
params.rho = rho;
params.Rs = Rs;
params.Rmemb = Rmemb;
params.Gmemb = Gmemb;
params.Gsoma = Gsoma;
params.Gdend = Gdend;
params.Rsoma = Rsoma;
params.Rdend = Rdend;
params.taum = taum;
params.Cm = Cm;
params.Rm = Rm;
params.radiusSoma = radiusSoma;
params.Ra = Ra;
params.diameterDendrite = diameterDendrite;
params.radiusDendrite = radiusDendrite;
params.lambda = lambda;
params.lengthDendrite = lengthDendrite;

% Close the L equation figure
close(h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [L, rho, LInit, h, finalMaxL, f, x, maxLOfInterest] = ...
                compute_L_and_rho (C0, C1, tau0, tau1, alpha1, ...
                                    minLs, maxLs, minRho)
%% Computes the electrotonic length (L) and the dendritic-to-somatic conductance ratio (rho)

% Hard-coded parameters
nTrialsPerRange = 5;

% Count the number of ranges to try
nMinLs = length(minLs);
nMaxLs = length(maxLs);
if nMinLs ~= nMaxLs
    error('Number of minimums and maximums do not match!');
else
    nRanges = nMaxLs;
end

% Count the total number of trials to try
nTrials = nRanges * nTrialsPerRange;

% All minLs and maxLs to try
minLsRepeated = reshape(transpose(repmat(minLs, 1, nTrialsPerRange)), ...
                        nTrials, 1);
maxLsRepeated = reshape(transpose(repmat(maxLs, 1, nTrialsPerRange)), ...
                        nTrials, 1);

if alpha1 == 0                    % tau0 == tau1, no dendrite
    % If alpha1 is zero, that means tau0 == tau1, so there is no dendrite
    %   Thus, the electrotonic length is zero
    LInit = 0;
    L = 0;
    rho = 0;
else
    % Otherwise, one must solve the equation in Johnston & Wu, p. 94
    %   which has an infinite number of solutions
    % We will use the first solution that converges in the range 
    %   L = minL ~ maxL

    % Make x a symbolic variable and f a symbolic function
    syms x f(x)

    % Define the equation to solve for L
    f(x) = abs(C1 / (2*C0*tau1/tau0 - C1)) - ...
                cot(alpha1*x) * (cot(alpha1*x) - 1/(alpha1*x));

    % Plot it
    h = figure; clf(h); hold on
    minLOfInterest = min(minLs(~isinf(minLs)));
    maxLOfInterest = max(maxLs(~isinf(maxLs)));
    fplot(f, [minLOfInterest, maxLOfInterest]);
    line([minLOfInterest, maxLOfInterest], [0, 0], ...
            'Color', 'r', 'LineStyle', '--');
    xlabel('x');
    ylabel('f(x)');
    title('Electrotonic length f(L) == 0')
    drawnow
            
    % Initialize a trial count
    ctTrial = 1;

    % Make an initial guess for L (cf. Johnston & Wu, p. 94)
    LInit = pi / alpha1;

    % Solve for L and rho once
    L = double(vpasolve(f(x) == 0, x, LInit));
    rho = compute_rho(alpha1, L);
    fprintf('L == %g\n', L);
    fprintf('rho == %g\n', rho);

    % Start with the first range
    minL = minLsRepeated(ctTrial);
    maxL = maxLsRepeated(ctTrial);

    % Find another solution if L or rho is out of range
    while ctTrial <= nTrials && ...
            (isempty(L) || L < minL || L > maxL || rho < minRho)  
        % If no more trials left, break and return NaNs
        if ctTrial == nTrials
            L = NaN;
            rho = NaN;
            break
        end
        
        % Increment the trial count
        ctTrial = ctTrial + 1;

        % Start with a random number in the range
        LInit = minL + (maxL - minL) * rand(1);

        % Solve for L and rho again
        L = double(vpasolve(f(x) == 0, x, LInit));
        rho = compute_rho(alpha1, L);

        % Get the next minL and maxL
        minL = minLsRepeated(ctTrial);
        maxL = maxLsRepeated(ctTrial);

        fprintf('L == %g\n', L);
        fprintf('rho == %g\n', rho);
    end

    % Update the final minL and maxL
    finalMinL = minL;
    finalMaxL = maxL;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function rho = compute_rho (alpha1, L)

rho = -alpha1 * cot(alpha1 * L) / coth(L);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%{
OLD CODE:

responseAmplitude = -(C0 + C1);

while L < 0 || rho < 0

% Start with a random number between 0 and 1
LInit = rand(1);

eqnForL = abs(C1/(2*C0*tau1/tau0-C1)) == ...
            cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x));

% eqnForL = generate_equation_for_L(C0, C1, tau0, tau1, alpha1, x);
% L = double(vpasolve(eqnForL, x, LInit));

function eqnForL = generate_equation_for_L (C0, C1, tau0, tau1, alpha1, x)

eqnForL = abs(C1/(2*C0*tau1/tau0-C1)) - ...
            cot(alpha1*x)*(cot(alpha1*x)-1/(alpha1*x));
L = double(vpasolve(eqnForL, x, LInit)); 

minL = 0;

[L, rho, LInit, h, finalMaxL, f, x, maxLOfInterest] = ...
    compute_L_and_rho(C0, C1, tau0, tau1, alpha1, minL, maxLs, minRho);

algorithmInfo.minL = minL;

%    fplot(f, [minL, maxL]);

% Version 15-1
minL = 0;
maxLs = [0.5; 1; 1.5; 2; 2.5; 3; 3.5; 4; 4.5; 5; Inf];

minLs = [0; 0.5; 1; 1.5; 2; 2.5; 3; 3.5; 4; 4.5; 5];
maxLs = [0.5; 1; 1.5; 2; 2.5; 3; 3.5; 4; 4.5; 5; Inf];

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%