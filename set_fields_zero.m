function [structure] = set_fields_zero (structure, varargin)
%% Set each field specified in varargin to zero and store previous values in a new field strcat(field, '_prev')
% Usage: [structure] = set_fields_zero (structure, varargin)
%   structure   - must include each field specified by varargin
%   varargin    - must include at least one argument
%               must be character arrays corresponding to numerical field values in structure
%
% Used by:
%       /media/adamX/m3ha/optimizer4gabab/optimizergui_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/optimizer_4compgabab.m
%       /media/adamX/m3ha/optimizer4gabab/fminsearch3_4compgabab.m
%       /home/Matlab/minEASE.m
%
% File History:
% 2016-10-06 Created
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
    error('Not enough input arguments, type ''help set_fields_zero'' for usage');
elseif ~isstruct(structure)
    error('First argument must be a structure array!');
elseif min(cellfun(@ischar, varargin)) < 1
    error('Second argument and beyond must be character arrays!');
end

%% Perform task
for f = 1:numel(varargin)
    field = varargin{f};

    % Check if field exists in structure
    if ~isfield(structure, field)
        error('The field ''%s'' doesn''t exist in the given structure!\n', field);
    end

    % Store previous value, i.e.,
    %   create structure.(field_prev) if not already exists
    field_prev = strcat(field, '_prev');    
    structure.(field_prev) = structure.(field);

    % Set field to zero
    structure.(field) = 0;
end
