function [paramvals] = change_params(tochange_names, tochange_vals, paramnames, paramvals, varargin)
%% Change parameter values
% Usage: [paramvals] = change_params(tochange_names, tochange_vals, paramnames, paramvals, varargin)
% Outputs:	paramvals	- updated paramvals
% Arguments: 	tochange_names	- the name(s) of the parameter(s) to change
%				must be a string/char vec or a cell array of strings/char vecs
%		tochange_vals	- the values(s) of the parameter(s) to change
%				must be a numeric array
%		paramnames	- a cell array of all parameter names
%				must be a cell array of strings or character arrays
%		paramvals	- a numeric array of all parameter values
%				must be a numeric array of same length as paramnames
%		varargin	- 'ExperimentName': name of the experiment of interest
%				must be an unambiguous, case-insensitive match to one of the following: 
%					'RTCl'	- chloride accumulation/extrusion model for Peter; updates 'cp_amp', 'cp_per', 'cp_num'
%					'm3ha'	- GABA-B receptor network
%					'noexp'	- no experiment provided; do nothing
%				default == 'noexp'
%
% Requires:	
%		/home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%		/home/Matlab/Adams_Functions/update_params.m
%
% Used by:	
%		/media/adamX/RTCl/neuronlaunch.m
%		/media/adamX/m3ha/m3ha_launch1.m
%
% 2017-05-03 Moved from update_params.m
% 2017-05-03 Moved trial number update back to neuronlaunch.m
% 2019-05-03 Changed so that it reads tochange_names and tochange_vals

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
	error('Not enough input arguments, type ''help change_params'' for usage');
end

% Add required inputs to an input Parser
iP = inputParser;
addRequired(iP, 'tochange_names', ...			% the name(s) of the parameter(s) to change
	@(x) validateattributes(x, {'char', 'string', 'cell'}, {'nonempty'}));
addRequired(iP, 'tochange_vals', ...			% the values(s) of the parameter(s) to change
	@(x) validateattributes(x, {'numeric'}, {'nonempty'}));
addRequired(iP, 'paramnames', ...			% a cell array of all parameter names
	@(x) assert(iscell(x) && (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
		'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'paramvals', ...			% a numeric array of all parameter values
	@(x) validateattributes(x, {'numeric', 'logical'}, {'vector'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'ExperimentName', 'noexp', ...	% name of the experiment of interest
	@(x) any(validatestring(x, {'RTCl', 'm3ha', 'noexp'})));

% Read from the input Parser
parse(iP, tochange_names, tochange_vals, paramnames, paramvals, varargin{:});
experimentname = validatestring(iP.Results.ExperimentName, {'RTCl', 'm3ha', 'noexp'});

% Check relationships between arguments
if length(paramnames) ~= length(paramvals)
	error('length of paramnames and paramvals are unequal!!');
end

%% Change parameter(s)
if iscell(tochange_names)
	% One or more parameters to change
	ntochange = numel(tochange_names);	% number of parameters to change
	for p = 1:ntochange
		indp = find_ind_str_in_cell(tochange_names{p}, paramnames);
		paramvals(indp) = tochange_vals(p);
	end
else
	% Only one parameter to change
	indp = find_ind_str_in_cell(tochange_names, paramnames);
	paramvals(indp) = tochange_vals;
end

%% Update dependent parameters for particular experiments
[paramvals] = update_params(paramnames, paramvals, 'ExperimentName', experimentname);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:                                          

function [paramvals, pchanged] = change_params(trialn, k, pname, paramnames, paramsinit, params_min, params_inc, params_islog)

indp = find_ind_str_in_cell(pname, paramnames);
paramvals = paramsinit;
pchanged = pname;
if params_islog(indp)
	paramvals(indp) = params_min(indp) * params_inc(indp)^(k-1);
else
	paramvals(indp) = params_min(indp) + (k-1) * params_inc(indp);
end

trialn		- current trial number
%				must be an integer between 1 and length(pchnames)
%% Update trial number
indtrialn = find_ind_str_in_cell('trialn', paramnames);
paramvals(indtrialn) = trialn;							% update trial number

%}
