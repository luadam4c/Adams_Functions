function intersection = intersectm (arrays)
%% Find the intersection of multiple arrays
% Usage: intersection = intersectm (arrays)
% Outputs:	intersection	- intersection of input arrays
% Arguments:	arrays		- a cell array of arrays that will be intersected; if just an array, return the array
%				must be a cell array of input arrays that can be recognized by the built-in intersect function
% Used by:	
%		/home/Matlab/Adams_Functions/find_ind_str_in_cell.m


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
	error('Not enough input arguments, type ''help intersectm'' for usage');1
elseif isempty(arrays)
	error('Argument cannot be empty!');
end

%% Return the array if not a cell array
if ~iscell(arrays)
	intersection = arrays;
end

%% Return the intersection of the contents of the cell array
numarrays = numel(arrays);		% number of arrays
intersection = arrays{1};
if numarrays > 1
	for k = 2:numarrays
		intersection = intersect(intersection, arrays{k});
	end
end
