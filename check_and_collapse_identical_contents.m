function variable = check_and_collapse_identical_contents (array, arrayname)
%% Checks if a cell array or array has identical contents and collapse it to one copy of the content
% USAGE: variable = check_and_collapse_identical_contents (array, arrayname)
%
% Used by:
%		/media/adamX/RTCl/tuning_curves.m
%
% 2017-04-17 Created
% TODO: Make arrayname a Parameter-Value pair
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
%%% TODO

if isnumeric(array)
	if numel(unique(array)) ~= 1		% number of unique elements is not 1
		error('Elements of %s are not identical!', arrayname);
	else
		variable = unique(array);
	end

elseif iscell(array)
	array_cp = cell(size(array));		% stores repeats of the first element of array
	for k = 1:numel(array_cp)
		array_cp{k} = array{1};		% set each element to be the first element of array
	end
	if numel(unique(cellfun(@isequal, array, array_cp))) ~= 1
		error('Elements of %s are not identical!', arrayname);
	else
		variable = array{1};
	end
end
