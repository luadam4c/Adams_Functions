function [istart, istart_ind] = find_istart (tvec0, ivec0s, istartwin, plotflag, vvec0s, outfolder, filebase)
%% Finds time of current application from a series of current vectors
% Usage: [istart, istart_ind] = find_istart (tvec0, ivec0s, istartwin, plotflag, vvec0s, outfolder, filebase)
% Arguments:	
%		tvec0		- original time vector, must be a column vector in ms
%		ivec0s		- original current vectors
%				must be a matrix where each column is an original current trace in pA
%		istartwin	- (opt) Window in which istart would lie, e.g. [1000 1100]
%				must be within range of tvec0
%				default == [tbase tvec0(end)]
%		plotflag	- (opt) whether to plot current traces or not
%				must be 0 or 1
%				default == 0
%		vvec0s		- (opt) original voltage vectors for plotting
%				must be a matrix where each column is an original voltage trace in mV
%		outfolder 	- (opt) directory to place outputs, e.g. '/media/adamX/m3ha/data_dclamp/take4/test_sweeps/'
%				must be a directory
%				default == pwd
%		filebase 	- (opt) base of filename (without extension), e.g. 'A100110_0008_18'
%				must be a char array
%				default == 'unnamed'
%
% Requires:	
%		cd/check_subdir.m
%
% Used by:
%		/media/shareX/share/Adam/Sample_files_from_Katie/test_sweeps.m
%
% 2016-09-22 Adapted from dclampDataExtractor.m
% 2016-10-05 Accounted for the case that tvec0 doesn't start from 0
% 2016-10-16 Added back close(h)
% 2016-10-27 Replace each directory with directories{k}
% 2016-11-02 Moved check directories to check_subdir.m
% 

%% Parameters used for data analysis
mafw = 0.3;	% width in ms for the moving-average filter
slth = 5;	% slope threshold for finding istart

%% Directory (ies) for placing figures
directories = {'/istart/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
	error('Not enough input arguments, type ''help find_istart'' for usage');
elseif isempty(tvec0) || isempty(ivec0s)
	error('First two inputs cannot be empty!');
elseif ~isnumeric(tvec0) || ~isnumeric(ivec0s)
	error('First two inputs must be numbers or numeric arrays!');
elseif ~isequal(size(tvec0, 1), size(ivec0s, 1))
	error('Time and current vectors do not have the same length!');
elseif size(tvec0, 2) > 1
	error('tvec0 must be a column vector!');
elseif nargin >= 3 && length(istartwin) < 2
	error('istartwin must have a start and an end!');
end

%% Extract info from data
sims = tvec0(2) - tvec0(1);					% sampling interval in ms
nswps = size(ivec0s, 2);
tbase = tvec0(1) - sims;

%% Check more arguments
if nargin >= 3 && (istartwin(1) < tbase || istartwin(2) > tvec0(end))
	error('istartwin out of range!');
elseif nargin >= 4 && (plotflag ~= 1 && plotflag ~= 0 && plotflag ~= false && plotflag ~= true)
	error('plotflag out of range!');
elseif nargin >= 5 && ~isequal(size(ivec0s), size(vvec0s))
	error('ivec0s & vvec0s do not have the same size!');
elseif nargin >= 6 && ~ischar(outfolder)
	error('outfolder must be a character array!');
elseif nargin >= 7 && ~ischar(filebase)
	error('filebase must be a character array!');
end

%% Set defaults for optional arguments
if nargin < 3
	istartwin = [tbase tvec0(end)];
end
if nargin < 4
	plotflag = 0;
end
if nargin < 5
	vvec0s = [];
end
if nargin < 6
	outfolder = pwd;
end
if nargin < 7
	filebase = 'unnamed';
end

%% Display standard output header
fprintf('FINDING time of current application for %s ...\n', filebase);
fprintf('Sampling interval == %g ms\n', sims);
fprintf('Number of sweeps == %d\n', nswps);

%% Find istart
ndp_mafw = round(mafw/sims);
istdvec = std(ivec0s, 0, 2);					% Compute the standard deviation over all sweeps
istdvec_filt = smooth(istdvec, ndp_mafw);			% Smooth with a moving average filter spanning 5 ms
istdvec_filt_slope = diff(istdvec_filt)/sims;			% Slope of smoothed standard deviation
ind_begin = find(tvec0 >= istartwin(1), 1);
ind_end = find(tvec0 <= istartwin(2), 1, 'last');
ind = ind_begin:ind_end;					% indices of interest
i_first_change_pt = find(istdvec_filt_slope(ind) > slth, 1);	% Slope of standard deviation is positive, 
								% 	using 0.1 as the slope threshold
if isempty(i_first_change_pt)
	istart_ind = ind(1);					% For the rest of the analysis to work smoothly
	istart = NaN;						% To detect
else
	istart_ind = i_first_change_pt + ind(1) - 1;		% index of current start
	istart = istart_ind * sims;				% time of current start in ms
end
fprintf('Index of current application == %d\n', istart_ind);
fprintf('Time of current application == %g ms\n', istart);

%% Plot current and voltage traces
if plotflag
	% Check if needed directories exist in outfolder
	check_subdir(outfolder, directories);

	% Plot current traces
	h = figure(3000);
	set(h, 'Visible', 'off');
	set(h, 'Name', 'Current analysis');
	clf(h);
	if ~isempty(vvec0s)
		subplot(2,1,1);
	end
	for swp = 1:nswps					% Plot each current trace
		plot(tvec0(ind), ivec0s(ind, swp)); hold on; 
	end
	plot(tvec0(ind), istdvec(ind), 'LineWidth', 2, 'Color', 'r');
	plot(tvec0(ind), istdvec_filt(ind), 'LineWidth', 2, 'Color', 'g');
	plot(tvec0(ind), istdvec_filt_slope(ind), 'LineWidth', 2, 'Color', 'b');
	if ~isempty(i_first_change_pt)
		for swp = 1:nswps
			plot(tvec0(istart_ind), ivec0s(istart_ind, swp), ...
			    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
		end
	end
	title(['Current analysis for ', strrep(filebase, '_', '\_')]);
	xlabel('Time (ms)')
	ylabel('Current (pA)')
	xlim(istartwin);
%	ylim([-30 30]);				
	if  ~isempty(vvec0s)
		subplot(2,1,2);
		for swp = 1:nswps		% Plot each voltage trace
			plot(tvec0(ind), vvec0s(ind, swp)); hold on; 
		end
		if ~isempty(i_first_change_pt)
			for swp = 1:nswps
				plot(tvec0(istart_ind), vvec0s(istart_ind, swp), ...
				    'LineWidth', 2, 'Color', 'k', 'Marker', 'x', 'MarkerSize', 12);
			end
		end
		xlabel('Time (ms)')
		ylabel('Voltage (mV)')
		xlim(istartwin);
%		ylim([-90 40]);				
	end
	figname = fullfile(outfolder, directories{1}, [filebase, '_istart', '.png']);
	saveas(h, figname);
	close(h);
end

%{
%% OLD CODE

ind = round(istartwin(1)/sims):round(istartwin(2)/sims);	% indices of interest % NOT ROBUST


%}
