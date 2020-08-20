function [structure] = restore_fields (structure, varargin)
%% Set each field specified in varargin to previous values from the field strcat(field, '_prev')
% Usage: [structure] = restore_fields (structure, varargin)
% Arguments:
%   structure    - must include field & strcat(field, '_prev') for each field specified by varargin
%   varargin    - must include at least one argument
%              must be character arrays corresponding to numerical field values in structure
%
% Used by:
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_optimizergui_4compgabab.m
%       cd/m3ha_optimizer_4compgabab.m
%       cd/minEASE.m

% File History:
% 2016-10-06 Created
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
    error('Not enough input arguments, type ''help restore_fields'' for usage');
elseif ~isstruct(structure)
    error('First argument must be a structure array!');
elseif min(cellfun(@ischar, varargin)) < 1
    error('Second argument and beyond must be character arrays!');
end

%% Perform task
for f = 1:numel(varargin)
    field = varargin{f};
    fieldPrev = strcat(field, '_prev');

    % Check if field and its previous value exist in structure
    if ~isfield(structure, field)
        error('The field ''%s'' doesn''t exist in the given structure!\n', field);
    elseif ~isfield(structure, fieldPrev)
        error('The field ''%s'' doesn''t exist in the given structure!\n', fieldPrev);
    end

    % Restore field to previous value
    structure.(field) = structure.(fieldPrev);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
