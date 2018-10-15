function [ipeak_time, ipeak_ind, ipeak_amp, ipeak_delay] = find_IPSC_peak (tvec0, ivec0s, istart, ipeakwin, plotflag, outfolder, filebase)
%% Finds time of current peak from a an inhibitory current trace (must be negative current)
% Usage: [ipeak_time, ipeak_ind, ipeak_amp, ipeak_delay] = find_IPSC_peak (tvec0, ivec0s, istart, ipeakwin, plotflag, outfolder, filebase)
% Outputs:	
%		ipeak_time	- time of ipeak (tvec0(ipeak_ind)) in ms (a row vector)
% 		ipeak_ind	- index of ipeak relative to ivec0s (a row vector)
%		ipeak_amp	- amplitude of ipeak in pA (a row vector)
%		ipeak_delay	- delay of ipeak (ipeak_time - istart) in ms (a row vector)
% Arguments:
%		tvec0		- original time vector, must be a column vector in ms
%		ivec0s		- original current vector in pA
%				must be a numeric array with columns the same length as tvec0
%		istart		- (opt) time relative to which the peak delay is computed
%				usually time of IPSC application, e.g. 1000
%				must be within range of tvec0
%				default == tbase, i.e., tvec0(1) - (tvec0(2) - tvec0(1))
%		ipeakwin	- (opt) Window in which ipeak would lie, e.g. [1000 1300]
%				must be within range of tvec0
%				default == [tbase tvec0(end)]
%		plotflag	- (opt) whether to plot current trace or not
%				must be 0 or 1
%				default == 0
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
%		/media/adamX/m3ha/data_dclamp/dclampDataExtractor.m
%		/media/adamX/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m
%
%
% 2016-10-05 Adapted from find_istart.m and dclampDataExtractor.m
% 2016-10-13 Can now read more than one current vector at a time
% 2016-10-16 Added back close(h)
% 2016-10-27 Replace each directory with directories{k}
% 2016-11-02 Moved check directories to check_subdir.m
% 

%% Directories for placing figures
directories = {'/IPSCpeak/'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 2
	error('Not enough input arguments, type ''help find_IPSC_peak'' for usage');
elseif isempty(tvec0) || isempty(ivec0s)
	error('First two inputs cannot be empty!');
elseif ~isnumeric(tvec0) || ~isnumeric(ivec0s)
	error('First two inputs must be numbers or numeric arrays!');
elseif ~isequal(size(tvec0, 1), size(ivec0s, 1))
	error('Time and current vectors do not have the same length!');
elseif size(tvec0, 2) > 1
	error('tvec0 must be a column vector!');
elseif nargin >= 4 && length(ipeakwin) < 2
	error('ipeakwin must have a start and an end!');
end

%% Extract info from data
sims = tvec0(2) - tvec0(1);					% sampling interval in ms
nswps = size(ivec0s, 2);
tbase = tvec0(1) - sims;

%% Check more arguments
if nargin >= 4 && (ipeakwin(1) < tbase || ipeakwin(2) > tvec0(end))
	error('ipeakwin out of range!');
elseif nargin >= 5 && (plotflag ~= 1 && plotflag ~= 0 && plotflag ~= false && plotflag ~= true)
	error('plotflag out of range!');
elseif nargin >= 6 && ~ischar(outfolder)
	error('outfolder must be a character array!');
elseif nargin >= 7 && ~ischar(filebase)
	error('filebase must be a character array!');
end

%% Set defaults for optional arguments
if nargin < 3
	istart = tbase;
end
if nargin < 4
	ipeakwin = [tbase tvec0(end)];	
end
if nargin < 5
	plotflag = 0;
end
if nargin < 6
	outfolder = pwd;
end
if nargin < 7
	filebase = 'unnamed';
end

%% Display standard output header
% fprintf('FINDING time of inhibitory current peak for %s ...\n', filebase);
% fprintf('Sampling interval == %g ms\n', sims);

%% Find ipeak
ind_begin = find(tvec0 >= ipeakwin(1), 1);
ind_end = find(tvec0 <= ipeakwin(2), 1, 'last');
ind = ind_begin:ind_end;			% indices of interest
ivec0s_part = ivec0s(ind, :);			% Use original current trace(s)
ipeak_amp = zeros(1, nswps);			% amplitude of ipeak in pA
ipeak_ind = zeros(1, nswps);			% index of ipeak relative to ivec0s
parfor swp = 1:nswps				% FOR each sweep
	[vtemp, itemp] = min(ivec0s_part(:, swp));
	ipeak_amp(swp) = vtemp;				% ipeak_amp is the amplitude of ipeak
	ipeak_ind(swp) = (ind_begin - 1) + itemp;	% index of ipeak relative to ivec0s
	% fprintf('Current peak found at index %g with amplitude = %g pA\n\n', ipeak_ind(swp), ipeak_amp(swp));
end
ipeak_time = tvec0(ipeak_ind)';			% time of ipeak in ms
ipeak_delay = ipeak_time - istart;		% delay of ipeak (ipeak_time - istart) in ms

%% Plot current trace
if plotflag == 1
	% Check if needed directories exist in outfolder
	check_subdir(outfolder, directories);

	% Plot current trace
	h = figure(4000);
	set(h, 'Visible', 'off');
	set(h, 'Name', 'IPSC peak amplitude analysis');
	clf(h);
	for swp = 1:nswps			% Plot each current trace and mark peak amplitude
		plot(tvec0(ind), ivec0s(ind, swp)); hold on; 
		plot(tvec0(ipeak_ind(swp)), ipeak_amp(swp), 'Marker', 'x', 'MarkerSize', 12);
	end
	xlabel('Time (ms)')
	ylabel('Current (pA)')
	title(['IPSC peak amplitude analysis for ', strrep(filebase, '_', '\_')]);
	figname = fullfile(outfolder, directories{1}, [filebase, '_IPSCpeak', '.png']);
	saveas(h, figname);
	close(h);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

IPSC_ind = round(ipscpwin(1)/sims);     % Assume no IPSC offset % NOT ROBUST
ind = IPSC_ind:round(ipscpwin(2)/sims);     % indices of interest   % NOT ROBUST

    [ipeak_amp, itemp1] = min(ivec0s_part1);    % ipeak_amp is the amplitude of ipeak
    ipeak_ind = (ind(1) - 1) + itemp1;      % index of ipeak relative to ivec0s
    ipeak_time = tvec0(ipeak_ind);          % time of ipeak in ms

    close(h);


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%