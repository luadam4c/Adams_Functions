function print_all_fields (structure)
%% Print everything in a structure
% Usage: print_all_fields (structure)
%	structure	- must be a structure
%
% Used by:
%		/home/Matlab/Adams_Functions/find_passive_params.m
%
% 2016-11-02 Created
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
	error(['Not enough input arguments,', ...
		' type ''help print_all_fields'' for usage']);
elseif ~isstruct(structure)
	error('First argument must be a structure array!');
end

%% Perform task
fprintf('Printing structure ''%s'' \n', inputname(1));
fields = fieldnames(structure);
nfields = numel(fields);
for f = 1:nfields
	field_value = structure.(fields{f});
	if isnumeric(field_value)
		fprintf('%s == %g\n', fields{f}, field_value);
	elseif ischar(field_value)
		fprintf('%s == %s\n', fields{f}, field_value);		
	else
		display(field_value);
	end
end
fprintf('\n');
