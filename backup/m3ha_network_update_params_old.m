function paramValues = m3ha_network_update_dependent_params (paramNames, paramValues, varargin)
%% Update dependent parameters for particular experiments
% Usage: paramValues = m3ha_network_update_dependent_params (paramNames, paramValues, varargin)
% Outputs:
%       paramValues    - updated paramValues
% Arguments:
%       paramNames  - a cell array of all parameter names
%                   must be a cell array of strings or character arrays
%       paramValues   - a numeric array of all parameter values
%                   must be a numeric array of same length as paramNames
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
%        cd/find_in_strings.m
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
addRequired(iP, 'paramNames', ...
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'paramValues', ...
    @(x) validateattributes(x, {'numeric', 'logical'}, {'vector'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'ExperimentName', 'noexp', ...
    @(x) any(validatestring(x, validExperiments)));

% Read from the input Parser
parse(iP, paramNames, paramValues, varargin{:});
experimentname = validatestring(iP.Results.ExperimentName, validExperiments);

% Check relationships between arguments
if length(paramNames) ~= length(paramValues)
    error('length of paramNames and paramValues are unequal!!');
end

%% Update dependent parameters, part I
switch experimentname
case {'RTCl', 'm3ha'}
    % Note: Must be consistent with m3ha_network_launch.m
    % Find indices of parameters that influence other parameters
    indREdia = find_in_strings('REdiam', paramNames);
    indstimf = find_in_strings('stimFreq', paramNames);
    indstimd = find_in_strings('stimDur', paramNames);

    % Find indices of parameters that are influenced by other parameters
    indcpamp = find_in_strings('cpAmp', paramNames);
    indcpper = find_in_strings('cpPer', paramNames);
    indcpnum = find_in_strings('cpNum', paramNames);

    % Update current pulse period (ms)
    paramValues(indcpper) = floor(1000/paramValues(indstimf));

    % Update number of current pulses                                        
    paramValues(indcpnum) = ceil(paramValues(indstimd)/paramValues(indcpper));
case 'noexp'
    % Do nothing
end

% Update current pulse amplitude (nA)
switch experimentname
case 'RTCl'
    paramValues(indcpamp) = 4*(paramValues(indREdia)/10)^2;
case 'm3ha'
    paramValues(indcpamp) = 0.2*(paramValues(indREdia)/10)^2;
case 'noexp'
    % Do nothing
end

%% Update GABAB-receptor-dependent parameters
switch experimentname
case 'm3ha'
    % Find indices of parameters that influence other parameters
    indpCond = find_in_strings('pCond', paramNames);
    indgIncr = find_in_strings('gIncr', paramNames);

    % Find indices of parameters that are influenced by other parameters
    indgabaaGmax = find_in_strings('TCgabaaGmax', paramNames);
    indgababAmp = find_in_strings('TCgababAmp', paramNames);
    indgababTrise = find_in_strings('TCgababTrise', paramNames);
    indgababTfallFast = find_in_strings('TCgababTfallFast', paramNames);
    indgababTfallSlow = find_in_strings('TCgababTfallSlow', paramNames);
    indgababW = find_in_strings('TCgababW', paramNames);
    
    % Update parameters that are influenced by other parameters
    pCond = paramValues(indpCond);
    gIncr = paramValues(indgIncr);
    if paramValues(indgabaaGmax) ~= 0             % Not in bicuculline mode
        paramValues(indgabaaGmax) = gabaaGmaxTemplate(pCond) * 2 * gIncr/100;  
                                    % update maximal GABA-A conductance (uS)
                                    %   2 times that of GABA-B 
                                    %   based on Huguenard & Prince, 1994
    end
    paramValues(indgababAmp) = gababAmpTemplate(pCond) * gIncr/100;
                                    % update GABA-B conductance amplitude (uS)
    paramValues(indgababTrise) = gababTriseTemplate(pCond);
                                    % update rising phase time constant (ms)
    paramValues(indgababTfallFast) = gababTfallFastTemplate(pCond);
                                    % update fast decay time constant (ms)
    paramValues(indgababTfallSlow) = gababTfallSlowTemplate(pCond);
                                    % update slow decay time constant (ms)
    paramValues(indgababW) = gababWeightDefault(pCond);
                                    % update weight of the fast decay

case {'RTCl', 'noexp'}
    % Do nothing
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
