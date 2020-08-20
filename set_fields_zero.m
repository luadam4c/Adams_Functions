function structure = set_fields_zero (structure, varargin)
%% Set each field specified in varargin to zero and store previous values in a new field strcat(field, '_prev')
% Usage: structure = set_fields_zero (structure, varargin)
%   structure   - must include each field specified by varargin
%   varargin    - must include at least one argument
%               must be character arrays corresponding to numerical field values in structure
%
% Used by:
%       cd/m3ha_fminsearch3.m
%       cd/m3ha_neuron_choose_best_params.m
%       cd/m3ha_optimizergui_4compgabab.m
%       cd/m3ha_optimizer_4compgabab.m
%       cd/minEASE.m

% File History:
% 2016-10-06 Created
% 2019-12-03 Now ignores fields if they don't exist

% TODO: Make optional argument
verbose = false;

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
        if verbose
            fprintf(['The field ''%s'' doesn''t exist ', ...
                    'in the given structure!\n'], field);
        end
    else
        % Store previous value, i.e.,
        %   create structure.(fieldPrev) if not already exists
        fieldPrev = strcat(field, '_prev');    
        structure.(fieldPrev) = structure.(field);

        % Set field to zero
        structure.(field) = 0;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
