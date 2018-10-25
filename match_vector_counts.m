function [vecs1, vecs2] = match_vector_counts(vecs1, vecs2, varargin)
%% Matches a set of vectors to another set of vectors so that they have equal number of vectors
% Usage: [vecs1, vecs2] = match_vector_counts(vecs1, vecs2, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       vecs1       - new first set of vectors
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
%       vecs2       - new second set of vectors
%                   specified as a numeric vector 
%                       or a cell array of numeric vectors
% Arguments:    
%       vecs1       - first set of vectors
%                   must be a numeric vector or a cell array of numeric vectors
%       vecs2       - second set of vectors
%                   must be a numeric vector or a cell array of numeric vectors
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
%       cd/count_vectors.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_dimensions.m
%
% Used by:    
%       cd/compute_residuals.m
%       cd/find_pulse_response_endpoints.m

% File History:
% 2018-10-11 Created by Adam Lu
% 2018-10-23 Now matches vectors both ways
% 2018-10-23 Added 'ForceCellOutputs' as an optional argument
% 2018-10-24 Added 'MatchDimensions' as an optional argument
% TODO: Deal with vectors as an array (many columns)
% 

%% Hard-coded parameters

%% Default values for optional arguments
forceCellOutputsDefault = false;    % don't force as cell array by default
matchDimensionsDefault = false;     % don't force match dimensions by default

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
addRequired(iP, 'vecs1', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs1 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));
addRequired(iP, 'vecs2', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['vecs2 must be either a numeric array', ...
                    'or a cell array of numeric arrays!']));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ForceCellOutputs', forceCellOutputsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MatchDimensions', matchDimensionsDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, vecs1, vecs2, varargin{:});
forceCellOutputs = iP.Results.ForceCellOutputs;
matchDimensions = iP.Results.MatchDimensions;

%% Do the job
% Make sure vecs1 and vecs2 are both cell arrays
%   if one of them is a cell array
if iscell(vecs1) && isnumeric(vecs2)
    vecs2 = {vecs2};
elseif isnumeric(vecs1) && iscell(vecs2)
    vecs1 = {vecs1};
end

% Make sure vecs1 and vecs2 have the same length
if iscell(vecs1) && iscell(vecs2)
    % Count the number of vectors
    nVecs1 = numel(vecs1);
    nVecs2 = numel(vecs2);

    % If unequal, match numbers
    if nVecs1 ~= nVecs2
        if nVecs1 ~= 1 && nVecs2 ~= 1        
            % If both are not one, return error
            error(['Either vecs1 and vecs2 have ', ...
                    'equal numbers of vectors or one of them ', ...
                    'must have just one vector!']);
        elseif nVecs1 == 1
            vecs1 = repmat(vecs1, size(vecs2));
        elseif nVecs2 == 1
            vecs2 = repmat(vecs2, size(vecs1));
        else
            error('Error in code logic!');
        end
    end

    % If requested, match dimensions
    if matchDimensions
        vecs1 = match_dimensions(vecs1, size(vecs2));
    end
else
    % Force outputs to be cell arrays if requested
    if forceCellOutputs
        vecs1 = {vecs1};
        vecs2 = {vecs2};
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Force both vecs2 and vecs1 to be column cell arrays
%   Note: This allows the user to provide a cell array for one 
%           and a numeric array for the other
vecs2 = force_column_cell(vecs2);
vecs1 = force_column_cell(vecs1);

% Count the number of vecs2
nToBeMatched = count_vectors(vecs2);

% Count the number of vecs1
nOrig = count_vectors(vecs1);

% Make sure nOrig matches up with nToBeMatched
if nOrig ~= nToBeMatched
    if nOrig == 0 
        % Generate a cell array of empty vecs2
        vecs1 = cell(nToBeMatched, 1);
    elseif nOrig == 1
        % Generate a cell array of the same vector
        vecs1 = repmat(vecs1, nToBeMatched, 1);
    else
        error(['# of vecs1 provided must be 0, 1 or ', ...
                'the same as vecs2!!']);
    end
else
    % Numbers already matched
    vecs1 = vecs1;
end

vecs1 = match_dimensions(vecs1, vecs2);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%