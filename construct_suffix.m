function [final_suffix] = construct_suffix (varargin)
%% Constructs final suffix based on optional suffices and/or Name-Value pairs
% Usage: [final_suffix] = construct_suffix (varargin)
% Outputs:	final_suffix 	- a string (may be empty) that is a final suffix
% Arguments:	varargin	- 'Suffices': suffix(ces) to add to filebase
%				must be a string/character array or a cell array of strings/character arrays
%				default == 'nosuffices'
%				- 'NameValuePairs': Name-Value pairs that are changed
%				must be a 2-element cell array whose first element is a string/char array or cell array 
%							and whose second element is a numeric array
%				default == {'nopairs', NaN}
%		
%
% Used by:
%		cd/construct_fullfilename.m
%		/media/adamX/RTCl/raster_plot.m
% 
% 2017-05-04 Moved from construct_fullfilename.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Add parameter-value pairs to the input Parser
iP = inputParser;
addParameter(iP, 'Suffices', '', ...			% suffix(ces) to add to filebase
	@(x) assert(ischar(x) || iscell(x) && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))) || isstring(x) , ...
		'Suffices must be either a string/character array or a cell array of strings/character arrays!'));
addParameter(iP, 'NameValuePairs', {'', NaN}, ...	% Name-Value pairs that are changed
	@(x) assert(iscell(x) && numel(x) == 2 ...
			&& (ischar(x{1}) || iscell(x{1}) || isstring(x{1})) ...
			&& isnumeric(x{2}), ...
		['NameValuePairs must be a 2-element cell array whose first element ', ...
		'is a string/char array or cell array and whose second element is a numeric array!']));

% Read from the input Parser
parse(iP, varargin{:});
suffices = iP.Results.Suffices;
namevaluepairs = iP.Results.NameValuePairs;

%% Find all suffices
if isempty(suffices) && isempty(namevaluepairs{1})
	allsuffices = '';
else
	% Initialize a cell array for all suffices
	numsuffices = 0;
	if iscell(suffices)
		numsuffices = numsuffices + numel(suffices);
	elseif ~isempty(suffices)
		numsuffices = numsuffices + 1;
	end
	if iscell(namevaluepairs{1})
		numsuffices = numsuffices + numel(namevaluepairs{1});
	elseif ~isempty(namevaluepairs{1})
		numsuffices = numsuffices + 1;
	end
	allsuffices = cell(1, numsuffices);				% stores all suffices
	ct = 0;
	
	% Add premade suffices if necessary
	if ~isempty(suffices)		% if premade suffices are provided
		if iscell(suffices)
			for s = 1:numel(suffices)
				ct = ct + 1;
				allsuffices{ct} = suffices{s};
			end
		else
			ct = ct + 1;
			allsuffices{ct} = suffices;
		end
	end

	% Add premade Name-Value pairs if necessary
	if ~isempty(namevaluepairs{1})	% if Name-Value pairs that are changed are provided
		if iscell(namevaluepairs{1})	% if there might be more than one Name-Value pairs provided
			for p = 1:numel(namevaluepairs{1})
				ct = ct + 1;
				allsuffices{ct} = [namevaluepairs{1}{p}, '_', num2str(namevaluepairs{2}(p))];
			end
		else				% if there is only one Name-Value pair provided
			ct = ct + 1;
			allsuffices{ct} = [namevaluepairs{1}, '_', num2str(namevaluepairs{2})];
		end
	end
end

%% Construct final suffix
if isempty(allsuffices)			% if nothing provided
	% Final suffix is empty too
	final_suffix = allsuffices;
elseif ~isempty(allsuffices)		% suffix(ces) is(are) provided
	% Construct path based on directory, filebase, suffix(ces) and extension
	if iscell(allsuffices)		% if there might be more than one suffices provided
		% Construct final suffix by concatenating all suffices together with '_'
		final_suffix = allsuffices{1};
		if numel(allsuffices) > 1
			for s = 2:numel(allsuffices)
				final_suffix = [final_suffix, '_', allsuffices{s}];
			end
		end
	else				% if there is only one suffix provided
		% Final suffix is this suffix
		final_suffix = allsuffices;
	end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:


%}
