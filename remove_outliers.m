function [new_data, indices_left] = remove_outliers (old_data, wl2IQR, plotflag)
%% Remove outliers from data matrix and return new matrix
% Usage: [new_data, indices_left] = remove_outliers (old_data, wl2IQR, plotflag)
% Outputs:	new_data	- data matrix with outlying data points removed
%		indices_left	- row indices of original data matrix that were left in place
% Arguments:	old_data	- a data matrix with each column being a condition and each row being a data point
%				must be a numeric array
%		wl2IQR		- (opt) the ratio of whisker length to interquartile range
%				must be a single number
%				default == 1.5 (same as the Matlab function boxplot())
%		plotflag	- (opt) whether to plot box plots or not
%				must be 0 or 1
%				default == 0
%
% Used by:	
%		/media/adamX/m3ha/data_dclamp/compare_sse.m
%
% 2016-12-08 Created
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
	error('Not enough input arguments, type ''help remove_outliers'' for usage');
elseif ~isnumeric(old_data)
	error('First argument must be a numeric array!');
elseif nargin >= 2 && (~isnumeric(wl2IQR) || length(wl2IQR) ~= 1)
	error('wl2IQR must be a single number!');
elseif nargin >= 3 && ~(plotflag == 0 || plotflag == 1)
	error('plotflag must be either 0 or 1!');
end

%% Set defaults for optional arguments
if nargin < 2
	wl2IQR = 1.5;
end
if nargin < 3
	plotflag = 0;
end

%% Remove outliers
ndp = size(old_data, 1);			% total number of data points
Q = quantile(old_data, [0.25; 0.5; 0.75]);	% each row corresponds to a quartile; each column corresponds to a vector
q1 = Q(1, :);					% first quartile
q3 = Q(3, :);					% third quartile
IQR = q3 - q1;					% the interquartile range
highbar = q3 + wl2IQR * (q3 - q1);		% max of whisker
lowbar = q1 - wl2IQR * (q3 - q1);		% min of whisker
toleave = ones(1, ndp);				% whether to leave data point or not, initially all true
for k = 1:ndp
	if sum([sum(old_data(k, :) > highbar), sum(old_data(k, :) < lowbar)])
		toleave(k) = false;
	end
end
indices_left = find(toleave);
new_data = old_data(indices_left, :);

%% Plot boxplots for verification
if plotflag
	figure();
	boxplot(old_data);
	figure();
	boxplot(new_data);
end
