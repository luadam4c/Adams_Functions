function [vectors, vector_names] = structs2vecs(structs)
%% Converts a cell array of structs with equal numbers of fields to a column cell array of row vectors or cell arrays
% Usage: vectors = structs2vecs(structs)
% Output:
%	vectors		- a column cell array of row vectors/cell arrays corresponding to the fields of the structure
%	vector_names	- a column cell array of field names of the structure corresponding to the row vectors/cell arrays
% Arguments:
%	structs		- a cell array (either column or row) of structures that contain an equal number of fields
%
% Used by:	
%		cd/m3ha_plot_histograms_refine_threshold.m
%		cd/m3ha_optimizer_4compgabab.m
%
% 2016-12-05 Moved from PlotHistogramsRefineThreshold.m, originally Brian Truong's code
% 2017-01-25 Now stores field names in vector_names
% 2017-01-25 Now allows fields to have different data types
% TODO: Use struct2table, than convert
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1 
	error('A cell array of structs is required, type ''help structs2vecs'' for usage');
elseif isempty(structs) || ~iscell(structs)
	error('Argument must be a nonempty cell array of structs, type ''help PlotHistogramsRefineThreshold'' for usage');
end

%% Perform task
array_length = numel(structs);			               % number of structs in the cell array
vector_names = fieldnames(structs{1});		% the names of the output vectors corresponding to the fields in the structs
num_fields = numel(vector_names);		% number of fields in a struct
struct_infos = cell(1, array_length);		% stores column cell arrays containing info of each struct
matrix = zeros(num_fields, array_length);	% stores all info if numeric
cell_matrix = cell(num_fields, array_length);	% stores all info if some are not numeric
for k = 1:array_length
	if isstruct(structs{k})
		% Convert the structure into a column vector and store in a matrix
		struct_infos{k} = struct2cell(structs{k});	% returns structure values in the same order as fieldnames
		if isnumeric(struct_infos{k})
			matrix(:, k) = cell2mat(struct_infos{k});
		else
			for j = 1:num_fields
				cell_matrix{j, k} = struct_infos{k}{j};
			end
		end
	end
end
if isnumeric(struct_infos{1})
	vectors = mat2cell(matrix, ones(1, num_fields));
else
	vectors = cell(num_fields, 1);
	for j = 1:num_fields
		if isnumeric(cell_matrix{j, 1})
			vectors{j} = cell2mat(cell_matrix(j, :));
		else
			vectors{j} = cell_matrix(j, :);
		end
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
