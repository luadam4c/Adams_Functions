function xolotlObject = xolotl_set_simparams (xolotlObject, varargin)
%% Runs a simulation on a xolotl object
% Usage: xolotlObject = xolotl_set_simparams (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
% Note: Defaults are for xolotl version 12-Dec-2018
%       xolotlObject    - a created neuron
%                       must be a xolotl object
%       varargin    - 'TimeStep': output time step in ms
%                   must be a positive scalar
%                   default == 0.05 ms
%                   - 'TimeEnd': time end in ms
%                   must be a nonnegative scalar
%                   default == 5000 ms
%                   - 'SimTimeStep': simulation integration time step in ms
%                   must be a positive scalar
%                   default == same as 'TimeStep'
%                   - 'Temperature': temperature of the experiment 
%                                   in degrees Celsius
%                   must be a numeric scalar
%                   default == 11
%                   - 'InitialVoltage': initial voltage values (mV)
%                                       for each compartment
%                   must be a numeric vector
%                   default == -60 mV for all compartments
%                   - 'ClosedLoop': whether the final condition is used as the 
%                                   initial condition for the next simulation
%                   must be a logical or numeric binary scalar
%                   default == true
%                   - 'SolverOrder': order of numerical integration method
%                   must be one of:
%                       0 - standard solvers are used 
%                           (the Exponential Euler method 
%                               for single-compartment models 
%                               and the implicit Crank-Nicholson method
%                               for multi-compartment models)
%                       4 - Runge-Kutta method
%                   default == 0
%
% Requires:
%       cd/force_column_numeric.m
%       cd/match_row_count.m
%       cd/parse_xolotl_object.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Created by Adam Lu
% 2018-12-13 Now sets the temperature
% 2018-12-13 Added solverOrder
% 2018-12-13 Updated to use xolotl defaults
% 2018-12-13 Now sets initial voltage value
% TODO: approx_channels: whether approximations to computing gating functions should be used
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeStepDefault = [];       % set later
timeEndDefault = [];        % set later
simTimeStepDefault = [];    % set later
temperatureDefault = [];    % set later
initialVoltageDefault = []; % set later
closedLoopDefault = [];     % set later
solverOrderDefault = [];    % set later

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'xolotlObject');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TimeStep', timeStepDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'TimeEnd', timeEndDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));
addParameter(iP, 'SimTimeStep', simTimeStepDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'Temperature', temperatureDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));
addParameter(iP, 'InitialVoltage', initialVoltageDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addParameter(iP, 'ClosedLoop', closedLoopDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'scalar'}));
addParameter(iP, 'SolverOrder', solverOrderDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative', 'integer'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
timeStep = iP.Results.TimeStep;
timeEnd = iP.Results.TimeEnd;
simTimeStep = iP.Results.SimTimeStep;
temperature = iP.Results.Temperature;
initialVoltage = iP.Results.InitialVoltage;
closedLoop = iP.Results.ClosedLoop;
solverOrder = iP.Results.SolverOrder;

%% Preparation
% Parse the xolotl object
parsedParams = parse_xolotl_object(xolotlObject);

% Extract the number of compartments
compartments = parsedParams.compartments;
nCompartments = parsedParams.nCompartments;

%% Set simulation parameters
% Set the time step in ms
if ~isempty(timeStep)
    xolotlObject.dt = timeStep;
end

% Set the end time in ms
if ~isempty(timeEnd)
    xolotlObject.t_end = timeEnd;
end

% Set the simulation time step in ms
%   also if timeStep is newly set
if ~isempty(simTimeStep)
    xolotlObject.sim_dt = simTimeStep;
elseif ~isempty(timeStep)
    xolotlObject.sim_dt = timeStep;
end

% Set the temperature
if ~isempty(temperature)
    xolotlObject.temperature = temperature;
end

% Set whether to simulate in a closed loop
if ~isempty(closedLoop)
    xolotlObject.closed_loop = closedLoop;
end

% Set the method for numerical integration
if ~isempty(solverOrder)
    xolotlObject.solver_order = solverOrder;
end

%% Set initial voltages
if ~isempty(initialVoltage)
    % Force as a numeric column
    initialVoltage = force_column_numeric(initialVoltage);

    % Match the number of voltages to the number of compartments
    initialVoltage = match_row_count(initialVoltage, nCompartments);

    % Set initial voltages
    for iCompartment = 1:nCompartments
        % Get the current compartment name
        compartmentName = compartments{iCompartment};

        % Set the initial voltage for the current compartment
        xolotlObject.(compartmentName).V = initialVoltage(iCompartment);
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

timeStep = 0.05;
timeEnd = 5000;
simTimeStep = timeStep;
temperature = 11;
closedLoop = true;
solverOrder = 0;

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%