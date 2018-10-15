function [vecsMatched, vecsToBeMatched] = ...
            match_vector_numbers(vecsToMatch, vecsToBeMatched, varargin)
%% Matches a set of vecsToBeMatched to another set of vecsToBeMatched so that they become cell arrays of equal length
% Usage: [vecsMatched, vecsToBeMatched] = ...
%           match_vector_numbers(vecsToMatch, vecsToBeMatched, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vecsMatched - vector(s) matched that may be transformed
%                   specified as a cell array of numeric column vecsToBeMatched
%       vecsToBeMatched - vector(s) to be matched that may be transformed
%                   specified as a cell array of numeric column vecsToBeMatched
% Arguments:    
%       vecsToMatch - vector(s) to match
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vecsToBeMatched
%       vecsToBeMatched - vector(s) to be matched
%                   Note: If a cell array, each element must be a vector
%                         If an array, each column is a vector
%                   must be a numeric array or a cell array of numeric vecsToBeMatched
%
% Requires:
%       cd/count_vectors.m
%       cd/force_column_cell.m
%
% Used by:    
%       cd/find_pulse_response_endpoints.m

% File History:
% 2018-10-11 Created by Adam Lu
% TODO: Do not force into cell arrays if not necessary
% 

%% Hard-coded parameters

%% Default values for optional arguments

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

% Add required inputs to the Input Parser
addRequired(iP, 'vecsToMatch', ...
    @(x) isnumeric(x) || iscell(x) && all(cellfun(@isnumeric, x)) );
addRequired(iP, 'vecsToBeMatched', ...
    @(x) isnumeric(x) || iscell(x) && all(cellfun(@isnumeric, x)) );

% Read from the Input Parser
parse(iP, vecsToMatch, vecsToBeMatched, varargin{:});

%% Do the job
% Force both vecsToBeMatched and vecsToMatch to be column cell arrays
%   Note: This allows the user to provide a cell array for one 
%           and a numeric array for the other
vecsToBeMatched = force_column_cell(vecsToBeMatched);
vecsToMatch = force_column_cell(vecsToMatch);

% Count the number of vecsToBeMatched
nToBeMatched = count_vectors(vecsToBeMatched);

% Count the number of pulse vecsToBeMatched
nOrig = count_vectors(vecsToMatch);

% Make sure nOrig matches up with nToBeMatched
if nOrig ~= nToBeMatched
    if nOrig == 0 
        % Generate a cell array of empty vecsToBeMatched
        vecsMatched = cell(nToBeMatched, 1);
    elseif nOrig == 1
        % Generate a cell array of the same vector
        vecsMatched = repmat(vecsToMatch, nToBeMatched, 1);
    else
        error(['# of vecsToMatch provided must be 0, 1 or ', ...
                'the same as vecsToBeMatched!!']);
    end
else
    % Numbers already matched
    vecsMatched = vecsToMatch;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%