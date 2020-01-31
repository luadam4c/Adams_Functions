function [amp, tauRise, tauFallFast, tauFallSlow, weight] = ...
                m3ha_load_gabab_ipsc_params (varargin)
%% Loads GABA-B IPSC parameters
% Usage: [amp, tauRise, tauFallFast, tauFallSlow, weight] = ...
%               m3ha_load_gabab_ipsc_params (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [a, r, f1, f2, w] = m3ha_load_gabab_ipsc_params
%       [a, r, f1, f2, w] = m3ha_load_gabab_ipsc_params('AmpScaleFactor', 200)
%       [a, r, f1, f2, w] = m3ha_load_gabab_ipsc_params('AmpUnits', 'uS')
%
% Outputs:
%       TODO
%
% Arguments:
%       varargin    - 'AmpScaleFactor': amplitude scaling factor
%                   must be a numeric scalar
%                   default == 100%
%                   - 'AmpUnits': amplitude units
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'nS'    - nanoSiemens
%                       'uS'    - microSiemens
%                   default == 'nS'
%
% Requires:
%
% Used by:
%       cd/m3ha_compute_gabab_ipsc.m
%       cd/m3ha_network_plot_gabab.m
%       cd/m3ha_network_update_dependent_params.m
%       cd/m3ha_neuron_create_simulation_params.m
%       cd/m3ha_parse_sweep_settings.m

% File History:
% 2020-01-22 Created by Adam Lu
% 

%% Hard-coded constants
NS_PER_US = 1000;

%% Hard-coded parameters
validAmpUnits = {'nS', 'uS'};
ampOrig = [16.00; 24.00; 8.88; 6.32];                   % (nS)
tauRise = [52.00; 52.00; 38.63; 39.88];                 % (ms)
tauFallFast = [90.10; 90.10; 273.40; 65.80];            % (ms)
tauFallSlow = [1073.20; 1073.20; 1022.00; 2600.00];     % (ms)
weight = [0.952; 0.952; 0.775; 0.629];

%% Default values for optional arguments
ampScaleFactorDefault = 100;
ampUnitsDefault = 'nS';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AmpScaleFactor', ampScaleFactorDefault, ...
    @(x) assert(isempty(x) || isnumeric(x) && isscalar(x), ...
                ['AmpScaleFactor must be either empty ', ...
                    'or a numeric scalar!']));
addParameter(iP, 'ampUnits', ampUnitsDefault, ...
    @(x) any(validatestring(x, validAmpUnits)));

% Read from the Input Parser
parse(iP, varargin{:});
ampScaleFactor = iP.Results.AmpScaleFactor;
ampUnits = validatestring(iP.Results.ampUnits, validAmpUnits);

%% Do the job
% Scale the amplitudes as requested
amp = ampOrig .* (ampScaleFactor / 100);

% Convert amplitudes to desired units
switch ampUnits
    case 'nS'
        % Do nothiing
    case 'uS'
        % Convert to uS
        amp = amp ./ NS_PER_US;
    otherwise
        error('ampUnits unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
