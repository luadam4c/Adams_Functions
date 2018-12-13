function xolotlObject = xolotl_set_simparams (xolotlObject, varargin)
%% Runs a simulation on a xolotl object
% Usage: xolotlObject = xolotl_set_simparams (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       xolotlObject    - the created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - the created neuron
%                       must be a xolotl object
%       varargin    - 'TimeStep': output integration time step in ms
%                   must be a TODO
%                   default == TODO
%                   - 'TimeEnd': time end in ms
%                   must be a TODO
%                   default == TODO
%                   - 'SimTimeStep': simulation integration time step in ms
%                   must be a TODO
%                   default == same as 'TimeStep'
%                   - 'ClosedLoop': whether the final condition is used as the 
%                                   initial condition for the next simulation
%                   must be a TODO
%                   default == false
%
% Requires:
%       TODO
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Created by Adam Lu
% TODO: Set temperatures
% TODO: Set initial voltage value
% 

%% Hard-coded parameters

%% Default values for optional arguments
timeStepDefault = 0.1;      % sample at 10 kHz by default
timeEndDefault = 1350;      % simulate for 1350 ms by default
simTimeStepDefault = [];    % same as timeStep by default
closedLoopDefault = false;  % don't use the final condition as the initial
                            %   condition for the next simulation by default

%{
celsius = 33            // Temperature: Christine probably incubated her slices 
                        //    at 33 degrees Celsius
                        // Note: Desteshe used 36 degrees Celsius
ipscTimeOrig = 1000     // Time of IPSC application in the experiments (ms), 
                        //    original value
cpStart = 100           // Current pulse start time (ms)
                        // This is only approximate in the experiments, 
                        //  but the data trace is realigned 
                        //  in singleneuronfitting.m
timeToStabilize = 2000  // Time to ensure that state variables 
                        //  for all channels are stabilized
                        //  2011-07-20 to allow sim to settle
ipscDur = 7000          // Duration of IPSC application (ms), for fitting
                        // This is probably shorter what was applied 
                        //  in some of the sweeps, but IPSC duration in 
                        //  all sweeps should be longer than this
vcRs = 0.01             // Voltage clamp series resistance in MOhm
secondorder = 0         // Backward Euler (default)
%}

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
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'SimTimeStep', simTimeStepDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'ClosedLoop', closedLoopDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
timeStep = iP.Results.TimeStep;
timeEnd = iP.Results.TimeEnd;
simTimeStep = iP.Results.SimTimeStep;
closedLoop = iP.Results.ClosedLoop;

% Check relationships between arguments
% TODO

%% Preparation
% Set default simulation time step
if isempty(simTimeStep)
    simTimeStep = timeStep;
end

%% Set simulation parameters
% Define the time step in ms
xolotlObject.dt = timeStep;

% Define the end time in ms
xolotlObject.t_end = timeEnd;

% Define the simulation time step in ms
xolotlObject.sim_dt = simTimeStep;

% Define whether to simulate in a closed loop
xolotlObject.closed_loop = closedLoop;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%