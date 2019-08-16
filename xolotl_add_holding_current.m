function xolotlObject = xolotl_add_holding_current (xolotlObject, varargin)
%% Adds a holding current to a xolotl object
% Usage: xolotlObject = xolotl_add_holding_current (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       x = xolotl_create_model_soplata;
%       x = xolotl_add_holding_current(x, 'Amp', 10);
%       x.plot
%
% Outputs:
%       xolotlObject    - a created neuron with simulation parameters
%                       specified as a xolotl object
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters
%                       must be a xolotl object
%       varargin    - 'Amplitude': amplitude in nA 
%                   must be a numeric scalar
%                   default == 0 nA
%                   - Any other parameter-value pair for 
%                       xolotl_add_current_injection()
%
% Requires:
%       cd/xolotl_add_current_injection.m
%       cd/parse_xolotl_object.m
%
% Used by:
%       cd/m3ha_xolotl_test.m

% File History:
% 2018-12-12 Modified from xolotl_add_current_pulse.m
% 2019-08-15 Now uses xolotl_add_current_injections.m
% 

%% Hard-coded parameters

%% Default values for optional arguments
amplitudeDefault = 0;       % default amplitude in nA

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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'xolotlObject');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Amplitude', amplitudeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
amplitude = iP.Results.Amplitude;

% Keep unmatched arguments for the xolotl_add_current_injection() function
otherArguments = struct2arglist(iP.Unmatched);

%% Add the holding current vectors to the previously set current injections
xolotlObject = xolotl_add_current_injection(xolotlObject, ...
                    'CurrentVector', amplitude, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%