function [output1, output2] = function_template (arg1, arg2, arg3, outfolder, filebase)
%% Description of function_name
% Usage: [output1, output2] = function_template (arg1, arg2, arg3, outfolder, filebase)
% Outputs:	output1		- Description of output1
%		output2		- Description of output2
% Arguments:	arg1		- Description of arg1
%				must be a %%%
%				default == %%%
%		arg2		- Description of arg2
%				must be a %%%
%				default == %%%
%		arg3		- Description of arg3
%				must be a %%%
%				default == %%%
%		outfolder 	- (opt) directory to place outputs, e.g. '/media/shareX/share/'
%				must be a directory
%				default == pwd
%		filebase 	- (opt) base of filename (without extension), e.g. 'A100110_0008_18'
%				must be a char array
%				default == 'unnamed'
%
% File History:
% 201X-XX-XX Created
% 

%% Parameters used in the body of the function
% (Any number that was arbitrarily defined in the code should be place here)

%% Directory (ies) for placing outputs
directories = {'/dir1/', '/dir2/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 3
	error('Not enough input arguments, type ''help function_name'' for usage');
elseif isempty(arg1) || isempty(arg2)
	error('First two inputs cannot be empty!');
% elseif %% anything else
	% error('');
elseif nargin >= 4 && ~ischar(outfolder)
	error('outfolder must be a character array!');
elseif nargin >= 5 && ~ischar(filebase)
	error('filebase must be a character array!');
end

%% Set defaults for optional arguments
if nargin < 4
	outfolder = pwd;
end
if nargin < 5
	filebase = 'unnamed';
end

% Check if needed directories exist
for k = 1:numel(directories)
	if exist(fullfile(outfolder, directories{k}), 'dir') ~= 7
		mkdir(fullfile(outfolder, directories{k}));
		fprintf('New directory made: %s\n\n', fullfile(outfolder, directories{k}));
	end
end

%% Perform task

%%%
