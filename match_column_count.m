function arrayNew = match_column_count (arrayOld, nColumnsNew, varargin)
%% Expands or truncates an array to match a given number of columns (dimension #2)
% Usage: arrayNew = match_column_count (arrayOld, nColumnsNew, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       match_column_count([1, 2, 3], 6)
%       match_column_count([1; 2; 3], 6)
%       match_column_count([1, 2; 3, 4], 7)
%       match_column_count([1, 2; 3, 4], 7, 'ExpansionMethod', 'repeat')
%       match_column_count([1, 2; 3, 4], 7, 'ExpansionMethod', 'patchNaNs')
%
% Outputs:
%       arrayNew    - array matched
%                   specified as a numeric, cell or struct array
%
% Arguments:    
%       arrayOld    - array to match
%                   must be a numeric, cell or struct array
%       nColumnsNew    - new number of rows
%                   must be a positive integer scalar
%       varargin    - see 
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/match_format_vector_sets.m
%       cd/match_row_count.m
%
% Used by:
%       cd/plot_tuning_curve.m

% File History:
% 2020-01-02 Modified from match_row_count.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(create_error_for_nargin(mfilename));
end

%% Do the job
% Use match_row_count.m
arrayNew = match_row_count(arrayOld, nColumnsNew, 'DimToMatch', 2, varargin{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
