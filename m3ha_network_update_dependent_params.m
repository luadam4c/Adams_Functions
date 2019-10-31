function paramsTable = m3ha_network_update_dependent_params (paramsTable, varargin)
%% Update dependent parameters for particular experiments
% Usage: paramsTable = m3ha_network_update_dependent_params (paramsTable, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       paramsTable     - updated parameters table
%                       specified as a table
%
% Arguments:
%       paramsTable     - a parameters table with the row names being the 
%                           parameter names and a 'Value' column
%                       must be a table
%       varargin    - 'ExperimentName': name of the experiment of interest
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following: 
%                       'RTCl'  - chloride accumulation/extrusion model for 
%                                 Peter; updates 'cp_amp', 'cp_per', 'cp_num'
%                       'm3ha'  - GABA-B receptor network
%                       'noexp' - no experiment provided; do nothing
%                   default == 'noexp'
%
% Requires:    
%
% Used by:    
%        cd/m3ha_network_change_params.m
%

% File History:
% 2017-03-30 Created by Adam Lu
% 2017-05-03 Moved to /media/adamX/RTCl/
% 2017-05-03 Moved parts to m3ha_network_change_params.m
% 2017-05-03 Added inputParser scheme
% 2017-05-03 Added 'ExperimentName' as a parameter
% 2017-11-07 Fixed TCgabaaGmax for bicuculline mode
% 2018-05-08 Changed tabs to spaces and limited width to 80
% 2018-05-09 Fixed cpamp calculation for 'RTCl'
% 2019-10-31 Now uses tables

%% Hard-coded parameters
validExperiments = {'RTCl', 'm3ha', 'noexp'};

%% m3ha: Template parameters used for GABA-B IPSC conductance waveforms
% Note: amp in Christine's thesis was actually for 200 % G incr
%       This must be consistent with /media/adamX/m3ha/data_dclamp/CountSweeps.m
gabaaGmaxTemplate = [4.48; 6.72; 9.76; 4.48]./1000;    % maximal conductance (uS)
gababAmpTemplate = [16.00; 24.00; 8.88; 6.32]./1000;   % conductance amplitude (uS)
gababTriseTemplate = [52.00; 52.00; 38.63; 39.88]; % rising phase time constant (ms)
gababTfallFastTemplate = [90.10; 90.10; 273.40; 65.80];% fast decay time constant (ms)
gababTfallSlowTemplate = [1073.20; 1073.20; 1022.00; 2600.00];
                                                % slow decay time constant (ms)
gababWeightDefault = [0.952; 0.952; 0.775; 0.629];         % weight of the fast decay

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

% Add required inputs to an input Parser
addRequired(iP, 'paramsTable', ...
    @(x) validateattributes(x, {'table'}, {'2d'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'ExperimentName', 'noexp', ...
    @(x) any(validatestring(x, validExperiments)));

% Read from the input Parser
parse(iP, paramsTable, varargin{:});
experimentName = validatestring(iP.Results.ExperimentName, validExperiments);

%% Update dependent parameters, part I
switch experimentName
case {'RTCl', 'm3ha'}
    % Note: Must be consistent with m3ha_network_launch.m
    % Update current pulse period (ms)
    paramsTable('cpPer', 'Value') = ...
        floor(1000 / paramsTable('stimFreq', 'Value'));

    % Update number of current pulses
    paramsTable('cpNum', 'Value') = ...
        ceil(paramsTable('stimDur', 'Value')/paramsTable('cpPer', 'Value'));
case 'noexp'
    % Do nothing
end

% Update current pulse amplitude (nA)
switch experimentName
case 'RTCl'
    % Note: Must be consistent with m3ha_network_launch.m
    paramsTable('cpAmp', 'Value') = ...
        4 * (paramsTable('REdiam', 'Value') / 10) ^ 2;
case 'm3ha'
    paramsTable('cpAmp', 'Value') = ...
        0.2 * (paramsTable('REdiam', 'Value') / 10) ^ 2;
case 'noexp'
    % Do nothing
end

%% Update GABAB-receptor-dependent parameters
switch experimentName
case 'm3ha'
    % Note: Must be consistent with m3ha_network_launch.m 
    % Update parameters that are influenced by other parameters
    pCond = paramsTable('pCond', 'Value');
    gIncr = paramsTable('gIncr', 'Value');

    % Update maximal GABA-A conductance (uS)
    %   2 times that of GABA-B 
    %   based on Huguenard & Prince, 1994
    if paramsTable('TCgabaaGmax', 'Value') ~= 0      % Not in bicuculline mode
        paramsTable('TCgabaaGmax', 'Value') = ...
            gabaaGmaxTemplate(pCond) * 2 * gIncr/100;  
    end

    % Update GABA-B conductance amplitude (uS)
    paramsTable('TCgababAmp', 'Value') = gababAmpTemplate(pCond) * gIncr/100;

    % Update rising phase time constant (ms)
    paramsTable('TCgababTrise', 'Value') = gababTriseTemplate(pCond);

    % Update fast decay time constant (ms)
    paramsTable('TCgababTfallFast', 'Value') = gababTfallFastTemplate(pCond);

    % Update slow decay time constant (ms)
    paramsTable('TCgababTfallSlow', 'Value') = gababTfallSlowTemplate(pCond);

    % Update weight of the fast decay
    paramsTable('TCgababW', 'Value') = gababWeightDefault(pCond);
case {'RTCl', 'noexp'}
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
