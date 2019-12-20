function [oscillatoryPeriod, oscillatoryIndex] = m3ha_network_autocorrelogram (infolder, varargin)
%% Shows an m3ha_network_autocorrelogram for each set of neurons (each .spi file in the infolder)
% Usage: [oscillatoryPeriod, oscillatoryIndex] = m3ha_network_autocorrelogram (infolder, varargin)
%
% Arguments:
%    infolder   - the name of the directory containing the .syn files, e.g. '20170317T1127_Ggaba_0.01'
%               must be a directory
%    varargin   - 'FigTypes': figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%               could be anything recognised by the built-in saveas() function
%               (see isfigtype.m under Adams_Functions)
%               default == 'png'
%               - 'OutFolder': the name of the directory that the plots will be placed
%               must be a directory
%               default: same as infolder
%               - 'MaxNumWorkers': maximum number of workers for running NEURON 
%               must a positive integer
%               default: 20
%               - 'SingleTrialNum': number of single trial ran
%                                       or 0 if not just a single trial
%               must a positive integer
%               default: 0
%               - 'RenewParpool': whether to renew parpool every batch to release memory
%               must be logical 1 (true) or 0 (false)
%               default: true, but alway false if plotspikes == false
%               - 'PlotCorrelogram': whether to plot correlograms
%               must be logical 1 (true) or 0 (false)
%               default == true
%
% Requires:
%       cd/combine_strings.m
%       cd/tuning_curves.m
%       cd/tuning_maps.m
%       cd/m3ha_network_define_actmode.m
%       infolder/*.spi
%       infolder/*loopedparams.mat
%       infolder/['sim_params_', pstring, '.csv'] for all the possible parameter strings
%       /home/Matlab/Downloaded_Functions/dirr.m
%       /home/Matlab/Downloaded_Functions/subaxis.m
%       /home/Matlab/Adams_Functions/isfigtype.m
%       /home/Matlab/Adams_Functions/save_all_figtypes.m
%       /home/Matlab/Adams_Functions/plot_tuning_curve.m (through tuning_curves.m)
%       /home/Matlab/Adams_Functions/plot_tuning_map.m (through tuning_maps.m)
%       /home/Matlab/Adams_Functions/find_ind_str_in_cell.m
%       /home/Matlab/Adams_Functions/extract_looped_params.m
%
% Used by:
%       /home/Matlab/Adams_Functions/m3ha_network_plot_oscillations.m

% 2018-02-09 BT - Adapted from raster_plot.m
% 2018-04-03 BT - Correlogram made using weighted, binned spikes 
% 2018-04-04 BT - Oscillatory period using peak prom
% 2018-04-07 BT - Removed period with peak prom, now smooth data to find peak. Only use data after end of stimulation
% 2018-04-10 BT - Fixed axes
% 2018-04-24 BT - Smoothed with moving average filter, find first peak by filtering out small prominences

bins = 5;   % bin for average trace (ms)
smoothing = 50;    % data points to use for smoothing 
labels = {'Left Electrode', 'Middle Electrode', 'Right Electrode'};
elecPos = [100/6, 100/6*3, 100/6*5] + 0.5;
colors = {'b', 'r', 'g'};
colheader = {'Left Period', 'Middle Period', 'Right Period', 'Left Index', 'Middle Index', 'Right Index'}; % for csv

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath_custom(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                    % for dirr.m & subaxis.m
addpath_custom(fullfile(functionsdirectory, '/Adams_Functions/'));
                                    % for isfigtype.m, find_ind_str_in_cell.m 
                                    %   extract_looped_params.m

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Add required inputs to an Input Parser
iP = inputParser;
addRequired(iP, 'infolder', @isdir);            % the name of the directory containing the .syn files

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'FigTypes', 'png', ...         % figure type(s) for saving; e.g., 'png', 'fig', or {'png', 'fig'}, etc.
    @(x) min(isfigtype(x, 'ValidateMode', true)));
addParameter(iP, 'OutFolder', '@infolder', @isdir); % the name of the directory that the plots will be placed
addParameter(iP, 'MaxNumWorkers', 20, ...       % maximum number of workers for running NEURON
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'RenewParpool', true, ...      % whether to renew parpool every batch to release memory
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SingleTrialNum', 0, ...       % number of single trial ran
    @(x) validateattributes(x, {'numeric'}, {'nonnegative'}));
addParameter(iP, 'PlotCorrelogram', true, ...        % whether to plot raster plots
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, infolder, varargin{:});
[~, figtypes] = isfigtype(iP.Results.FigTypes, 'ValidateMode', true);
outfolder = iP.Results.OutFolder;
maxnumworkers = iP.Results.MaxNumWorkers;
renewparpool = iP.Results.RenewParpool;
singletrialnum = iP.Results.SingleTrialNum;
plotcorrelogram = iP.Results.PlotCorrelogram;

% Change default arguments if necessary
if strcmp(outfolder, '@infolder')
    outfolder = infolder;
end

%% Find all .spi files
REfiles = dir(fullfile(infolder, 'RE*.spi'));
nREfiles = length(REfiles);
TCfiles = dir(fullfile(infolder, 'TC*.spi'));
nTCfiles = length(TCfiles);

%% Get loop parameters
[nump, pnames, plabels, pislog, pvalues, nperp, pchnames, pchvalues, ncells, actmode, loopmode] = extract_looped_params(infolder);

%% Check number of .spi files
ntrials = numel(pchnames);
if nREfiles ~= ntrials
    fprintf('Warning: Number of RE*.spi files and number of trials don''t match!\n');
end
if nTCfiles ~= ntrials
    fprintf('Warning: Number of TC*.spi files and number of trials don''t match!\n');
end

%% Create raster plot
REspikes = cell(ntrials, 1);        % raw spike data for RE neurons
TCspikes = cell(ntrials, 1);        % raw spike data for TC neurons
numActiveTC = zeros(ntrials, 1);    % number of cells activated
numActiveRE = zeros(ntrials, 1);    % number of cells activated
oscillatoryIndex = cell(ntrials, 1);        % oscillation index
oscillatoryPeriod = cell(ntrials, 1);   % oscillation period, cell = {left middle right}

ct = 0;                             % counts number of trials completed
poolobj = gcp('nocreate');          % get current parallel pool object without creating a new one
if isempty(poolobj)
    poolobj = parpool;              % create a default parallel pool object
    oldnumworkers = poolobj.NumWorkers;    % number of workers in the default parallel pool object
else
    oldnumworkers = poolobj.NumWorkers;    % number of workers in the current parallel pool object
end
numworkers = min(oldnumworkers, maxnumworkers);    % number of workers to use for running NEURON
if renewparpool
    delete(poolobj);                % delete the parallel pool object to release memory
end
while ct < ntrials                  % while not trials are completed yet
    if singletrialnum    % if running only one trial
        first = singletrialnum; % run that trial
    else
        first = ct + 1;         % first trial in this batch
    end
    if singletrialnum           % if running only one trial
        last = singletrialnum;      % run only that trial
    elseif renewparpool && ...
        ct + numworkers <= ntrials  % if memory is to be released
        last = ct + numworkers;     % limit the batch to numworkers
    else
        last = ntrials;             % run all trials at once
    end
    if renewparpool
        % Recreate a parallel pool object using fewer workers to prevent running out of memory
        poolobj = parpool('local', numworkers);    
    end
    parfor i = first:last
    %for i = first:last
        % Construct current parameter string
        pstring = '';               % initialize for parfor
        if iscell(pchvalues)
            pstring = combine_strings('NameValuePairs', {pchnames{i}, pchvalues{i}});
        elseif isnumeric(pchvalues)
            pstring = combine_strings('NameValuePairs', {pchnames{i}, pchvalues(i)});
        end
        pspifile = [pstring, '.spi'];

        % Load spike data
        jnowRE = 0;                 % initialize for parfor
        jnowTC = 0;                 % initialize for parfor
        for j = 1:nREfiles
            if ~isempty(strfind(REfiles(j).name, pspifile))
                jnowRE = j;         % index in REfiles of current file
                REspikes{i} = load(fullfile(infolder, REfiles(jnowRE).name));
            end
        end
        for j = 1:nTCfiles
            if ~isempty(strfind(TCfiles(j).name, pspifile))
                jnowTC = j;         % index in TCfiles of current file
                TCspikes{i} = load(fullfile(infolder, TCfiles(jnowTC).name));
            end
        end

        % Don't plot if file not found or no spikes
        stimcellIDs = 0;            % initialize for parfor    
        if isempty(REspikes{i}) && isempty(TCspikes{i})
            fprintf('The file %s does not exist or has no spikes!\n', pspifile);
            numActiveTC(i) = 0;     % number of cells activated is zero
            numActiveRE(i) = 0;        % number of cells activated is zero
        else
            % Extract parameters from sim_params file
            fid = fopen(fullfile(infolder, ['sim_params_', pstring, '.csv']));
            simfilecontent = textscan(fid, '%s %f %s', 'Delimiter', ',');
            paramnames = simfilecontent{1};
            params_val = simfilecontent{2};
            tstart = params_val(find_ind_str_in_cell('tstart', paramnames, 'SearchMode', 'exact'));
            tstop = params_val(find_ind_str_in_cell('tstop', paramnames, 'SearchMode', 'exact'));
            stim_start = params_val(find_ind_str_in_cell('stim_start', paramnames, 'SearchMode', 'exact'));
            stim_dur = params_val(find_ind_str_in_cell('stim_dur', paramnames, 'SearchMode', 'exact'));
            stim_freq = params_val(find_ind_str_in_cell('stim_freq', paramnames, 'SearchMode', 'exact'));
            fclose(fid);

            % Extract info
            if ~isempty(REspikes{i})
                REspikecelln = REspikes{i}(:, 1);
                REspiketimes = REspikes{i}(:, 2);
            else
                REspikecelln = [];
                REspiketimes = [];
            end
            if ~isempty(TCspikes{i})
                TCspikecelln = TCspikes{i}(:, 1);
                TCspiketimes = TCspikes{i}(:, 2);
            else
                TCspikecelln = [];
                TCspiketimes = [];
            end
            xlim1 = tstart;
            xlim2 = tstop;
            stim_start_plot = stim_start;
            stim_dur_plot = stim_dur;
            
            % Bin number of spikes by cell number
            % Create edges for the bins, bins = length(edges)-1
            edges = (stim_start+stim_dur):bins:tstop;

            % Store number of spikes per bin
            nSpikes = zeros(200,length(edges)-1);

            % 100 TC cells & 100 RE cells 
            for celln = 1:100
                % Bin all spikes by a single cell
                celltimeTC = TCspiketimes(find(TCspikecelln == celln));
                celltimeRE = REspiketimes(find(REspikecelln == celln));
                for binn = 2:length(edges)
                    % Cell specific spikes in bin
                    bin_tindicesTC = find(celltimeTC > edges(binn-1) & celltimeTC <= edges(binn));
                    bin_tindicesRE = find(celltimeRE > edges(binn-1) & celltimeRE <= edges(binn));
                    
                    % Total number of spikes by cell in bin
                    nSpikes(celln,binn) = length(bin_tindicesTC);
                    nSpikes(celln+100,binn) = length(bin_tindicesRE);
                end
            end

            % Weighted average
            % nSpikes = [2*cell# (RE & TC) x bin], val = number of spikes
            w = cell(length(elecPos),1);
            elec = cell(length(elecPos),1);
            for j = 1:length(elecPos)
                w{j} = [1./abs((1:100)' - elecPos(j)); 1./abs((1:100)' - elecPos(j))];
                elec{j} = (w{j}/sum(w{j}))' * nSpikes;
            end
            
            % Calculate correlogram differences for TC
            correlogram = cellfun(@(x) xcorr(x), elec, 'UniformOutput', false);
            
            % Halve xcorr and normalize to 1
            correlogram = cellfun(@(x) x((length(x)+1)/2:end)/max(x), correlogram, 'UniformOutput', false); 

            % Smooth correlogram 
            correlogram = cellfun(@(x) smooth(x, smoothing), correlogram, 'UniformOutput', false);

            % Find oscillatory period by finding distance from stimulation end to first peak
            [pks,lcs,wid,prom] = cellfun(@(x) findpeaks(x),correlogram, 'UniformOutput', false);
            
            % Using prominence to filter out minor peaks
            prom_av = cellfun(@(x) mean(x), prom, 'UniformOutput', false);
           
            % Store temp variables for period & index
            oscillatoryPeriod_temp = cell(length(elecPos),1);
            oscillatoryIndex_temp = cell(length(elecPos),1);
            
            % For every electrode
            for j = 1:length(elecPos)
                % Only using peaks with a greater than average prominence   
                prom_ind = find(prom{j} >= prom_av{j});
                
                % Extract relevant peaks and indices
                pks{j} = pks{j}(prom_ind);
                lcs{j} = lcs{j}(prom_ind);

                % Heavier weighting towards initial peaks
                pos_weight = (pks{j}).*(length(edges)-lcs{j}).^1.5;

                % Get peak, supposed to be 'period'
                [~, pos_ind] = max(pos_weight);

                % Get max peak & index
                pks{j} = pks{j}(pos_ind);
                lcs{j} = lcs{j}(pos_ind);
                
                % If no peaks -> no period or index
                if length(pks{j}) == 0   
                    oscillatoryPeriod_temp{j} = 0;
                    pks{j} = NaN;
                    lcs{j} = NaN;
                    oscillatoryIndex_temp{j} = 0;
                else
                    % Period is location of first peak in time scale
                    oscillatoryPeriod_temp{j} = lcs{j}(1)*bins-bins;
                    
                    % Index given by (A(peak)-A(trough))/A(max) [A(max) is 1 after normalization]
                    [minpk, ~] = min(correlogram{j}(1:lcs{j}(1)));
                    oscillatoryIndex_temp{j} = (pks{j}(1)-minpk)/1;
                end
            end
            oscillatoryPeriod{i} = oscillatoryPeriod_temp;
            oscillatoryIndex{i} = oscillatoryIndex_temp;
            if plotcorrelogram == true
                % Plot average trace
                f = figure(40000);
                %f.Visible = 'off';
                clf(f);
                hold on;
                for j = 1:length(elec)
                    % Don't plot rows with no spike values
                    if all(elec{j} == 0)
                        elec{j}(elec{j} == 0) = NaN;
                    end
                    plot(edges, elec{j}, 'Color', colors{j});
                end
                % When no spikes at all, set axes to [0 1]
                if cellfun(@(x) all(isnan(x)), elec)
                    ylim([0 1]);
                end
                xlim([0 tstop]);
                xlabel('Spike time (ms)');
                ylabel(['Weighted spikes in bin = ' num2str(bins) ' ms']);
                legend(labels, 'Location', 'northeast');
                figbaseWext = strrep(REfiles(jnowRE).name, 'RE_', '');
                figbase = strrep(figbaseWext, '.spi', '');
                title([strrep(figbase, '_', '\_') ' Weighted Spikes']);
    %                suptitle(strrep(figbase, '_', '\_'));

                % Save figure
                figname = fullfile(outfolder, strrep(figbaseWext, '.spi', '_average_trace.png'));
                save_all_figtypes(f, figname, figtypes);
                close(f);

                % Plot correlogram
                h = figure(50000);
                %h.Visible = 'off';
                clf(h);
                hold on;
                p = zeros(length(elec),1);
                for j = 1:length(elec)
                    % x: end of stimulation to end of recording
                    p(j) = plot([0:bins:tstop-(stim_start+stim_dur)], correlogram{j}, 'Color', colors{j});
                    text(xlim1, 0.925-.05*(j-1), [labels{j} ' Oscillatory Period: ' num2str(oscillatoryPeriod{i}{j}) ' ms'], 'Color', colors{j});
                    plot(oscillatoryPeriod{i}{j},pks{j}(1),'*','Color',colors{j});
                    text(xlim1, 0.775-.05*(j-1), [labels{j} ' Oscillatory Index: ' num2str(oscillatoryIndex{i}{j})], 'Color', colors{j});
                end
                xlim([0 tstop]);
                xlabel(['Time lag (ms)']);
                ylim([0 1]);
                ylabel('');
                legend(p, labels, 'Location', 'northeast');
                figbaseWext = strrep(REfiles(jnowRE).name, 'RE_', '');
                figbase = strrep(figbaseWext, '.spi', '');
                title([strrep(figbase, '_', '\_') ' Autocorrelation']);
    %                suptitle(strrep(figbase, '_', '\_'));

                % Save figure
                figname = fullfile(outfolder, strrep(figbaseWext, '.spi', '_correlogram.png'));
                save_all_figtypes(h, figname, figtypes);
                close(h); 
            end  
        end
%        close all;
    end

    if renewparpool
        delete(poolobj);    % delete the parallel pool object to release memory
    end
    if singletrialnum
        ct = ntrials + 1;   % don't run again
    else
       ct = last;          % update number of trials completed
    end
end
if renewparpool
    % Recreate a parallel pool object using the previous number of workers
    poolobj = parpool('local', oldnumworkers);
end
% Converting cells to matrices
oscillatoryPeriodMat = zeros(ct,length(elecPos));
for i = 1:length(oscillatoryPeriod)
    for j = 1:length(elecPos)
        if isnan(oscillatoryPeriod{i}{j})
            oscillatoryPeriodMat(i,j) = NaN;
        else
            oscillatoryPeriodMat(i,j) = oscillatoryPeriod{i}{j};
        end
    end
end
% Converting cells to matrices
oscillatoryIndexMat = zeros(ct,length(elecPos));
for i = 1:length(oscillatoryIndex)
    for j = 1:length(elecPos)
        if isnan(oscillatoryIndex{i}{j})
            oscillatoryIndexMat(j,i) = NaN;
        else
            oscillatoryIndexMat(i,j) = oscillatoryIndex{i}{j};
        end
    end
end
rowheader = cell(1,length(pchnames));
for i = 1:length(pchnames)
    rowheader{i} = [pchnames{i} '_' num2str(pchvalues(i))];
end
dlmwrite_with_header(fullfile(infolder, 'oscillatoryPeriodIndex.csv'), [oscillatoryPeriodMat oscillatoryIndexMat], 'RowHeader', rowheader, 'ColumnHeader', colheader)
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

            % lcs = cellfun(@(x) x+spacingEst, lcs, 'UniformOutput', false);
spacingEst = 450 / bins;    % min peak time (ms) between peaks in correlogram. Divide by bins to get bin index spacing. Base = 450 / bins
sims = 0.1;     % sampling interval in ms, 1/sims must be integer
                % If adjusting sims, must also adjust xlabel of correlogram to proper units
                % RESOLUTION OF SPIKE TIMES
function [Rxx]=autom(x)
    % [Rxx]=autom(x)
    % This function Estimates the autocorrelation of the sequence of
    % random variables given in x as: Rxx(1), Rxx(2),…,Rxx(N), where N is
    % Number of samples in x.
    N = length(x);
    Rxx = zeros(1,N);
    for m = 1:N+1
        for n = 1:N-m+1
            Rxx(m) = Rxx(m) + x(n) * x(n + m - 1);
        end
    end
end
  [pks,lcs] = findpeaks(correlogramTC);
            if length(pks) < 2  % if not at least two peaks, no period
                oscillatoryPeriodTC(i) = 0;
            else
                [pks,lcs] = findpeaks(correlogramTC,1,'MinPeakDistance',spacingEst);    % find peaks in range of expected period
                if length(pks) < 2
                    oscillatoryPeriodTC(i) = 0;
                else
                    oscillatoryPeriodTC(i) = mean(diff(lcs)) * sims;    % average spacing between peaks is period
                end
            end
            [pks,lcs] = findpeaks(correlogramRE);
            if length(pks) < 2  % if not at least two peaks, no period
                oscillatoryPeriodRE(i) = 0;
            else
                [pks,lcs] = findpeaks(correlogramRE,1,'MinPeakDistance',spacingEst);    % find peaks in range of expected period
                if length(pks) < 2
                    oscillatoryPeriodRE(i) = 0;
                else
                    oscillatoryPeriodRE(i) = mean(diff(lcs)) * sims;    % average spacing between peaks is period
                end
            end 

            figbaseWext = strrep(TCfiles(jnowTC).name, 'TC_', '');
            correlogramTC = autom(TCcollapse);
            correlogramRE = autom(REcollapse);
            correlogramTC = correlogramTC / max(correlogramTC);    % normalize correlograms to 0-1
            correlogramRE = correlogramRE / max(correlogramRE);

            % Find ratio of currently activated cells to total cells
            TCcollapse = zeros(length(unique(TCspiketimes)),1);
            TCuniquetimes = unique(TCspiketimes)';
            for j = 1:length(TCcollapse)
                TCcollapse(j) = length(find(TCspiketimes == TCuniquetimes(j))) / numActiveTC(i);
            end
            TCcollapse = TCcollapse - 0.5;    % shift down to center y-range from -0.5 to 0.5

            REcollapse = zeros(length(unique(REspiketimes)),1);
            REuniquetimes = unique(REspiketimes)';
            for j = 1:length(REcollapse)
                REcollapse(j) = length(find(REspiketimes == REuniquetimes(j))) / numActiveRE(i);
            end
            REcollapse = REcollapse - 0.5;    % shift down to center y-range from -0.5 to 0.5
                        TCspiketimeshold = TCspiketimes;    % copy spike times for backup
            TCspiketimes(TCspiketimes < stim_start) = [];    % delete times before stimulation start
            TCcollapse = zeros((tstop - stim_start)/sims, 1);    % template for ratio vector
            for j = stim_start:sims:tstop-1    % from time stim start to end of recording
                TCcollapse(j/sims-stim_start+1) = length(find(TCspiketimes == j)) / numActiveTC(i);    % ratio of spikes at time to total num spikes
            end
            REspiketimeshold = REspiketimes;
            REspiketimes(REspiketimes < stim_start) = [];
            REcollapse = zeros((tstop - stim_start)/sims, 1);
            for j = stim_start:sims:tstop-1
                REcollapse(j/sims-stim_start+1) = length(find(REspiketimes == j)) / numActiveRE(i);
            end
            %{
            if numActiveTC(i) == 0
                oscillatoryPeriodTC(i) = 0;
                oscillatoryIndexTC(i) = 0;
                correlogramTC = NaN([samplingUB+1 1]);
            else
                correlogramTC = zeros(samplingUB/sims+1,1);
                for j = 0:sims:samplingUB
                    for k = 1:length(TCcollapse)-j
                        correlogramTC(j+1) = correlogramTC(j+1) + (TCcollapse(k) * TCcollapse(k + j));
                    end
                end
                correlogramTC = correlogramTC / max(correlogramTC);    % normalize
                [corrMaxPeaks, corrMaxLocs] = findpeaks(correlogramTC);
                [corrMinPeaks, corrMinLocs] = findpeaks(-correlogramTC);
                if isempty(corrMaxPeaks)
                    oscillatoryPeriodTC(i) = 0;
                    oscillatoryIndexTC(i) = 0;
                elseif isempty(corrMinPeaks)
                    oscillatoryPeriodTC(i) = corrMaxLocs(1);
                    oscillatoryIndexTC(i) = 0;
                else
                    oscillatoryPeriodTC(i) = corrMaxLocs(1);
                    oscillatoryIndexTC(i) = (corrMaxPeaks(1) - corrMinPeaks(1)) / correlogramTC(1);
                end
            end
            % Calculate correlogram differences for RE
            if numActiveRE(i) == 0
                oscillatoryPeriodRE(i) = 0;
                oscillatoryIndexRE(i) = 0;
                correlogramRE = NaN([samplingUB+1 1]);
            else
                correlogramRE = zeros(samplingUB/sims+1,1);
                for j = 0:sims:samplingUB
                    for k = 1:length(REcollapse)-j
                        correlogramRE(j+1) = correlogramRE(j+1) + (REcollapse(k) * REcollapse(k + j));
                    end
                end
                correlogramRE = correlogramRE / max(correlogramRE);    % normalize
                [corrMaxPeaks, corrMaxLocs] = findpeaks(correlogramRE);
                [corrMinPeaks, corrMinLocs] = findpeaks(-correlogramRE);
                if isempty(corrMaxPeaks)
                    oscillatoryPeriodRE(i) = 0;
                    oscillatoryIndexRE(i) = 0;
                elseif isempty(corrMinPeaks)
                    oscillatoryPeriodRE(i) = corrMaxLocs(1);
                    oscillatoryIndexRE(i) = 0;
                else
                    oscillatoryPeriodRE(i) = corrMaxLocs(1);
                    oscillatoryIndexRE(i) = (corrMaxPeaks(1) - corrMinPeaks(1)) / correlogramRE(1);
                end
            end

                        RERErad = params_val(find_ind_str_in_cell('RERErad', paramnames, 'SearchMode', 'exact'));
            RETCrad = params_val(find_ind_str_in_cell('RETCrad', paramnames, 'SearchMode', 'exact'));
            actcellID = params_val(find_ind_str_in_cell('actcellID', paramnames, 'SearchMode', 'exact'));
            figbaseWext = strrep(TCfiles(jnowTC).name, 'TC_', '');
           %pos_weight = cellfun(@(x,y) (x.^.7).*((length(edges)-y).^2.5), pks, lcs, 'UniformOutput', false);
                %[~, prom_ind] = max(pos_weight{j});

%}
%}


