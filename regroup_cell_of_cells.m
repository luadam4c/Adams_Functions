function vecsRegrouped = regroup_cell_of_cells (vecs, varargin)
%% Regroup a cell array of cell arrays of numeric vectors
% Usage: vecsRegrouped = regroup_cell_of_cells (vecs, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       TODO
%
% Outputs:
%       vecsRegrouped   - regrouped vectors
%                       specified as a cell array of 
%                           cell arrays of numeric vectors
%
% Arguments:
%       vecs        - original vectors
%                   must be a cell array of cell arrays of numeric vectors
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%                   - Any other parameter-value pair for extract_columns()
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/extract_columns.m
%
% Used by:
%       cd/m3ha_simulate_population.m

% File History:
% 2020-XX-XX Created by TODO or Adapted from TODO
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
iP.KeepUnmatched = true;                        % allow extraneous options

% Add required inputs to the Input Parser
addRequired(iP, 'vecs');

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default);

% Read from the Input Parser
parse(iP, vecs, varargin{:});
% param1 = iP.Results.param1;

% Keep unmatched arguments for the extract_columns() function
otherArguments = iP.Unmatched;

%% Do the job
vecsRegrouped = extract_columns(vecs, 'TreatCnvAsColumns', true, ...
                                'OutputMode', 'single', otherArguments);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
