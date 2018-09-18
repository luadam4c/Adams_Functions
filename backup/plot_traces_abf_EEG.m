function [alldata, si] = plot_traces_abf_EEG (abffilename, left, right)
%% (May be obsolete) Plots traces all in one place for EEG data
% Usage: [alldata, si] = plot_traces_abf_EEG (abffilename, left, right)
%        abffilename either is the full address or must be in current directory
%        left is the start of the time interval of interest
%        right is the end of the time interval of interest

% File history
% 20160615 - created
% 20160922 - renamed

if isempty(strfind(abffilename, '/'))
	current_folder = pwd;
	filename_full = fullfile(current_folder, abffilename);
else
	filename_full = abffilename;
end

% Parameters to set

% Load abf file, si is in us
[alldata, si] = abf2load(filename_full);

% Find data parameters
nsweeps = size(alldata, 2);     % Number of sweeps
ntps = size(alldata, 1);        % Number of time points

% Set up vector for timepoints in seconds
tps = ( si*1e-6 : si*1e-6 : ntps * si*1e-6 )';

% Plot raw data (each sweep individually)
for i = 1:nsweeps
    cdata = alldata(:, i);
    figure(i)
    plot(tps, cdata, 'k')
    hold on
    % axis([0 ntps*si/1000 min(alldata(:,i))*1.2 max(alldata(:,i)*1.2)])
    axis([left right min(alldata(:,i))*1.2 max(alldata(:,i)*1.2)])
    title(sprintf('Data between %.1f s and %.1f s', left, right));
    xlabel('Time (s)')
	ylabel('EEG amplitude (uV)')
    saveas(gcf, strcat(cd, sprintf('/%.1f_%.1f_Sweep#%d', left, right, i), '.png'), 'png')
    hold off
end
