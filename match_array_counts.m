function [arrays1, arrays2] = match_array_counts (arrays1, arrays2, varargin)
%% Matches a set of arrays to another set of arrays so that they have equal number of arrays
% Usage: [arrays1, arrays2] = match_array_counts (arrays1, arrays2, varargin)
% Explanation:
%       TODO
%       cf. match_format_vector_sets.m
% Example(s):
%       [a, b] = match_array_counts({1:5, 2:6}, 1:5)
%       [a, b] = match_array_counts({1:5, [2:6]'}, 1:5)
%       [a, b] = match_array_counts([[1:5]', [2:6]'], [1:5]')
% Outputs:
%       arrays1     - new first set of arrays
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
%       arrays2     - new second set of arrays
%                   specified as a numeric array 
%                       or a cell array of numeric arrays
% Arguments:    
%       arrays1     - first set of arrays
%                   must be a numeric array or a cell array of numeric arrays
%       arrays2     - second set of arrays
%                   must be a numeric array or a cell array of numeric arrays
%       varargin    - 'ForceCellOutputs': whether to force outputs as 
%                                           cell arrays
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MatchDimensions': whether to always match dimensions
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   
%
% Requires:
%       cd/iscellnumeric.m
%       cd/match_dimensions.m
%
% Used by:
%

% File History:
% 2018-10-11 Created by Adam Lu
% 2018-10-23 Now matches arrays both ways
% 2018-10-23 Added 'ForceCellOutputs' as an optional argument
% 2018-10-24 Added 'MatchDimensions' as an optional argument
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputsDefault = false;    % don't force as cell array by default
matchDimensionsDefault = false;     % don't force match dimensions by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'arrays1', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['arrays1 must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'arrays2', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['arrays2 must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutputs', forceCellOutputsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchDimensions', matchDimensionsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, arrays1, arrays2, varargin{:});
forceCellOutputs = iP.Results.ForceCellOutputs;
matchDimensions = iP.Results.MatchDimensions;

%% Do the job
% Make sure arrays1 and arrays2 are both cell arrays
%   if one of them is a cell array
if iscell(arrays1) && isnumeric(arrays2)
    arrays2 = {arrays2};
elseif isnumeric(arrays1) && iscell(arrays2)
    arrays1 = {arrays1};
end

% Make sure arrays1 and arrays2 have the same length
if iscell(arrays1) && iscell(arrays2)
    % Count the number of arrays
    nVecs1 = numel(arrays1);
    nVecs2 = numel(arrays2);

    % If unequal, match numbers
    if nVecs1 ~= nVecs2
        if nVecs1 ~= 1 && nVecs2 ~= 1        
            % If both are not one, return error
            error(['Either arrays1 and arrays2 have ', ...
                    'equal numbers of arrays or one of them ', ...
                    'must have just one array!']);
        elseif nVecs1 == 1
            arrays1 = repmat(arrays1, size(arrays2));
        elseif nVecs2 == 1
            arrays2 = repmat(arrays2, size(arrays1));
        else
            error('Error in code logic!');
        end
    end

    % If requested, match dimensions
    if matchDimensions
        arrays1 = match_dimensions(arrays1, size(arrays2));
    end
else
    % Force outputs to be cell arrays if requested
    if forceCellOutputs
        arrays1 = {arrays1};
        arrays2 = {arrays2};
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%       cd/count_vectors.m
%       cd/force_column_cell.m

% Force both arrays2 and arrays1 to be column cell arrays
%   Note: This allows the user to provide a cell array for one 
%           and a numeric array for the other
arrays2 = force_column_cell(arrays2);
arrays1 = force_column_cell(arrays1);

% Count the number of arrays2
nToBeMatched = count_vectors(arrays2);

% Count the number of arrays1
nOrig = count_vectors(arrays1);

% Make sure nOrig matches up with nToBeMatched
if nOrig ~= nToBeMatched
    if nOrig == 0 
        % Generate a cell array of empty arrays2
        arrays1 = cell(nToBeMatched, 1);
    elseif nOrig == 1
        % Generate a cell array of the same array
        arrays1 = repmat(arrays1, nToBeMatched, 1);
    else
        error(['# of arrays1 provided must be 0, 1 or ', ...
                'the same as arrays2!!']);
    end
else
    % Numbers already matched
    arrays1 = arrays1;
end

arrays1 = match_dimensions(arrays1, arrays2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
