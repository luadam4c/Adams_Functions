function combine_looped_params (infolder1, infolder2, outfolder)
%% TODO
% USAGE: combine_looped_params (infolder1, infolder2, outfolder)
%
% Requires:
%		infolder/*loopedparams.mat
%		cd/extract_looped_params.m
%
% 2017-04-17 Created
% 2017-04-17 Removed latency_cells_to_plot
% 2017-05-03 Moved to Adams_Functions
% 2017-05-03 Renamed function combine_loopparams -> combine_looped_params
% 2017-05-03 Changed loopedparamsfile 'loopvariables.mat' -> 'loopedparams.mat'

%% Must be consistent with create_looped_params.m & extract_looped_params.m
loopedparamsfile = 'loopedparams.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(create_error_for_nargin(mfilename));
end

% TODO: Input Parser scheme

%% Find the loopedparams.mat in the infolders and create one in the outfolder
m3 = matfile(fullfile(outfolder, loopedparamsfile), 'Writable', true);

%% Extract parameters
[nump1, pnames1, plabels1, pislog1, pvalues1, ~, ~, ~, ncells1, actmode1] = extract_looped_params (infolder1);
[nump2, pnames2, plabels2, pislog2, pvalues2, ~, ~, ~, ncells2, actmode2] = extract_looped_params (infolder2);

%% Check if invariants are the same
if nump1 ~= nump2
	error('nump are not the same!');
else
	nump = nump1;
end
if ~isequal(pnames1, pnames2)
	error('pnames are not the same!');
else
	pnames = pnames1;
end
if ~isequal(plabels1, plabels2)
	error('plabels are not the same!');
else
	m3.plabels = plabels1;
end
if ~isequal(pislog1, pislog2)
	error('pislogs are not the same!');
else
	m3.pislog = pislog1;
end
if ncells1 ~= ncells2
	error('ncells are not the same!');
else
	m3.ncells = ncells1;
end
if ~isequal(actmode1, actmode2)
	error('actmodes are not the same!');
else
	m3.actmode = actmode1;
end

%% Combine pvalues and compute nperp & ntrials
pvalues = cellfun(@sort, ...
		cellfun(@uniquetol, ...
		cellfun(@vertcat, pvalues1, pvalues2, 'UniformOutput', false), ...
			'UniformOutput', false), ...
			'UniformOutput', false);
nperp = cellfun(@length, pvalues);
ntrials = sum(nperp);

%% Check number of .spi files
files = dirr(outfolder, '.spi');
nfiles = length(files);
if nfiles ~= ntrials
	error('Number of .spi files and %s don''t match!', loopedparamsfile);
end

%% Create pchnames & pchvalues
pchnames = cell(1, nfiles);		% changed parameter name for each file
pchvalues = zeros(1, nfiles);		% changed parameter value for each file
ct = 0;					% counts number of files used
for p = 1:nump
	% Store parameter names and values for each file (to keep things in order)
	indices = (ct + 1):(ct + nperp(p));					% indices for this parameter
	pchnames(indices) = repmat(pnames(p), 1, length(indices));
	pchvalues(indices) = transpose(pvalues{p});

	% Update count of number of files used
	ct = ct + nperp(p);
end

%% Save nump, pnames, pvalues, nperp, ntrials, pchnames & pchvalues
m3.nump = nump;
m3.pnames = pnames;
m3.pvalues = pvalues;
m3.nperp = nperp;
m3.ntrials = ntrials;
m3.pchnames = pchnames;
m3.pchvalues = pchvalues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

loopfiles1 = dirr(infolder1, loopedparamsfile);
loopfiles2 = dirr(infolder2, loopedparamsfile);
if length(loopfiles1) > 1 || length(loopfiles2) > 1
    error('Too many versions of %s!', loopedparamsfile);
elseif length(loopfiles1) < 1 || length(loopfiles2) < 1
    error('%s is required in each infolder!', loopedparamsfile);
else
    m1 = matfile(fullfile(infolder1, loopfiles1(1).name));      % there should be only one loopedparams.mat
    m2 = matfile(fullfile(infolder2, loopfiles2(1).name));
end

if ~isequal(latency_cells_to_plot1, latency_cells_to_plot2)
    error('latency_cells_to_plot are not the same!');
else
    m3.latency_cells_to_plot = latency_cells_to_plot1;
end

loopedparamsfile_old = 'loopvariables.mat'; %%% TODO: Make this an option too for compatibility with old data

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
