function varargout = xolotl_simulate (xolotlObject, varargin)
%% Simulates a voltage clamp experiment
% Usage: varargout = xolotl_simulate (xolotlObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       xolotlObject    - a created neuron with simulation parameters set
%                       must be a xolotl object
%       varargin    - 'OutputType': the output type
%                   must be a one of the following:
%                       0 - separate matrices
%                       1 - organized in a structure
%                       2 - organized in a structure and 
%                               enable spike-detection in C++ code (2). 
%                       Note: The 0 option is useful when you only want a 
%                       few outputs or don't care about lots of variable names. 
%                       The latter options are useful when it's important to 
%                       keep all the output data organized. In addition, 
%                       option 2 saves memory at the expense of detail.
%                   default == 0
%
% Used by:
%       cd/xolotl_estimate_holding_current.m

% File History:
% 2018-12-13 Created by Adam Lu

%% Hard-coded parameters

%% Default values for optional arguments
outputTypeDefault = [];     % set later

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
addParameter(iP, 'OutputType', outputTypeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'nonnegative'}));

% Read from the Input Parser
parse(iP, xolotlObject, varargin{:});
outputType = iP.Results.OutputType;

%% Preparation
% Prevent a MATLAB warning on gcc version (works for R2018b and beyond)
warning('off', 'MATLAB:mex:GccVersion');

%% Do the job
% Decide on the output type
if ~isempty(outputType)
    % Set the output type
    xolotlObject.output_type = outputType;
else
    % Read the current output type
    outputType = xolotlObject.output_type;
end

% Simulate and return output based on output type
switch outputType
    case 0
        % Arguments are in separate matrices
        [varargout{1}, varargout{2}] = xolotlObject.integrate;
    case {1, 2}
        % Arguments are in one single structure
        varargout{1} = xolotlObject.integrate;
    otherwise
        error('The output type %s is unrecognized!', outputType);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%