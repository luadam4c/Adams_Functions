function phaseNumbers = combine_phase_numbers (allPhaseNumbers, varargin)
%% Combines phase number vectors using unique_custom, ignoring NaNs
% Usage: phaseNumbers = combine_phase_numbers (allPhaseNumbers, varargin)
% Explanation:
%       TODO
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
%                   - Any other parameter-value pair for the TODO() function
%
% Requires:
%       cd/cell2num.m
%       cd/count_samples.m
%       cd/create_error_for_nargin.m
%       cd/force_matrix.m
%       cd/struct2arglist.m
%       cd/unique_custom.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-27 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
% param1Default = [];             % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1    % TODO: 1 might need to be changed
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

% Keep unmatched arguments for the TODO() function
% otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Force as a matrix
allPhaseNumbers = force_matrix(allPhaseNumbers);

%% Do the job
% Get all unique phase numbers for each row in a cell array
uniquePhaseNumbers = ...
    arrayfun(@(x) unique_custom(allPhaseNumbers(x, :), 'IgnoreNan', true), ...
        transpose(1:size(allPhaseNumbers, 1)), 'UniformOutput', false);

% Count the number of unique phase numbers for each row
nUniquePhaseNumbers = count_samples(uniquePhaseNumbers);

% If any number of unique phase numbers is greater than 1, there is a conflict
if any(nUniquePhaseNumbers > 1)
    disp("Phase numbers don't match!");
    phaseNumbers = [];
else
    phaseNumbers = cell2num(uniquePhaseNumbers);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%