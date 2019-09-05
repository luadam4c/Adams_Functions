function phaseNumbers = combine_phase_numbers (allPhaseNumbers, varargin)
%% Combines (possibly multiple) phase number vectors into a single vector
% Usage: phaseNumbers = combine_phase_numbers (allPhaseNumbers, varargin)
% Explanation:
%       The default method uses unique_custom, ignoring NaNs
%
% Example(s):
%       combine_phase_numbers([[NaN; 1; 2; NaN], [NaN; 1; NaN; 3]])
%       combine_phase_numbers({[NaN 1 2 NaN]; [NaN 1 NaN 3]})
%       combine_phase_numbers({[NaN 1 2 NaN]; [NaN 1 1 3]})
%
% Outputs:
%       phaseNumbers    - combined phase numbers
%                       specified as a numeric column vector
%
% Arguments:
%       allPhaseNumbers     - vector(s) of phase numbers
%                   must be an array
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for compute_combined_trace
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/compute_combined_trace.m
%       cd/struct2arglist.m
%
% Used by:
%       cd/compute_population_average.m

% File History:
% 2019-08-27 Created by Adam Lu
% 2019-08-27 Now uses compute_combined_trace with 'unique' option
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'allPhaseNumbers');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, allPhaseNumbers, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the compute_combined_trace() function
otherArguments = struct2arglist(iP.Unmatched);

%% Apply the compute_combined_trace() function with the 'unique' option
phaseNumbers = compute_combined_trace(allPhaseNumbers, 'unique', ...
                                        otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%