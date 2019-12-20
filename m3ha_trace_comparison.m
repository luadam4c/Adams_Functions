% m3ha_trace_comparison.m
%% Compares current and conductance traces with theoretical traces
%
% Requires:    
%       cd/compute_gabab_conductance.m
%       cd/rescale_vec.m
%       cd/compute_elcurr.m

% File History
% 2016-09-04 - Created
% 2016-09-07 - Changed method of saving ioffset_tentative
% 2016-11-06 - Added I_theo2_corr
% 2016-11-07 - Moved code to compute_gabab_conductance.m, compute_elcurr.m
% 2016-11-07 - Added scaling factor, G_data_rescaled, Rs, V_corr1, V_corr2, I_theo1_corr & I_theo2_corr
% 2016-11-07 - I_theo1 now comes from G_data_rescaled instead of G_data
% 2016-11-09 - Modified Rs so that it is fitted
% 2016-11-09 - Added text labels
% 2016-11-10 - Created logheader & logvariables; made variables row vectors instead of column vectors
% 2016-12-08 - Changed the series resistance compensation to a fixed voltage offset

%% Flags
debugflag = 0;

%% Parameters
Voff0 = 1;        % Initial value for fixed voltage offset (mV)

%% Debug files
debug_files = {'G091810_0003_14', 'B092710_0009_25', 'G091810_0002_19'};
debug_files = {'A092810_0006', 'E091710_0000'};

%% Fixed parameters used in the experiments
IPSC_start = 1000;    % Supposed time of IPSC application (ms)
E_rev = -105;        % Reversal potential for GABA_B IPSC (mV)
E_rev_ljp = -115;    % Possible reversal potential for GABA_B IPSC used by Christine (mV)

%% Parameters used in analyses
precision = 0.5;    % precision of the scaling factor for correcting conductance traces

%% Headers and variables for plotting
logheader = {'Data filename', 'Conductance scaling factor', ...
        'Current scaling factor', ...
        'Second conductance scaling factor', ... 
        'Tentative current offset (ms)', ...
        'Fitted voltage offset using rescaled G data (mV)', ...
        'Fitted voltage offset using G theoretical (mV)', ...
        'Sum of squares error using rescaled G data', ...
        'Sum of squares error using G theoretical'};

logvariables = {'fnrow', 'condscale', 'currscale', ...
        'condscale2', 'ioffset_tentative', ...
        'Voff1', 'Voff2', 'sse1', 'sse2'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
addpath_custom(fullfile(functionsdirectory, '/Adams_Functions/'));        % for m3ha_find_ipsc_peak.m, m3ha_find_lts.m, etc.

%% Set folders for reading and saving files
if exist('/media/adamX/m3ha/', 'dir') == 7
    homedirectory = '/media/adamX/m3ha/';
elseif exist('/scratch/al4ng/m3ha/', 'dir') == 7
    homedirectory = '/scratch/al4ng/m3ha/';
else
    error('Valid homedirectory does not exist!');
end
if debugflag
%    infolder = fullfile(homedirectory, '/data_dclamp/take4/');
    infolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
    outmatfolder = fullfile(homedirectory, '/data_dclamp/take4/debug/');
    outfigfolder = fullfile(homedirectory, '/data_dclamp/take4/debug/trace_comparison/');
else
    infolder = fullfile(homedirectory, '/data_dclamp/take4/');
    outmatfolder = fullfile(homedirectory, '/data_dclamp/take4/');
    outfigfolder = fullfile(homedirectory, '/data_dclamp/take4/trace_comparison/');
end
if exist(outmatfolder, 'dir') ~= 7
    mkdir(outmatfolder);
end
if exist(outfigfolder, 'dir') ~= 7
    mkdir(outfigfolder);
end

%% Set a random figure number for the plots 
% (This prevents the generation of many figures, which takes up too much memory)
fign = 10^6 + floor(rand(1)*10^6);

%% Import data
m1 = matfile(fullfile(infolder, 'dclampdatalog_take4.mat'), 'Writable', true);
m2 = matfile(fullfile(infolder, 'dclampdatalog_take4_resave.mat'), 'Writable', true);
fnrow = m1.fnrow;
gabab_amp = m1.gabab_amp;
gabab_Trise = m1.gabab_Trise;
gabab_TfallFast = m1.gabab_TfallFast;
gabab_TfallSlow = m1.gabab_TfallSlow;
gabab_w = m1.gabab_w;
condscale = m2.condscale;
currscale = m2.currscale;
ntraces = numel(fnrow);

%% Compare each current and conductance trace with the corresponding theoretical trace
ioffset_tentative = zeros(1, ntraces);
condscale2 = zeros(1, ntraces);
Voff1 = zeros(1, ntraces);    % fitted voltage offset (mV) using G_data_rescaled
Voff2 = zeros(1, ntraces);    % fitted voltage offset (mV) using G_theo 
sse1 = zeros(1, ntraces);    % sum of squares error of current trace using G_data_rescaled
sse2 = zeros(1, ntraces);    % sum of squares error of current trace using G_theo
parfor k = 1:ntraces
%for k = 1:ntraces
    % Check filename if in debug mode
    filebase = strrep(fnrow{k}, '.mat', '');    % file base name
    setbase = filebase(1:12);            % set base name
    if debugflag && ~ismember(setbase, debug_files)
        continue;
    end

    % Import traces
    m3 = matfile(fullfile(infolder, '/matfiles/', fnrow{k}));
    t_data = m3.d_orig(:, 1);
    G_data = m3.d_orig(:, 2);    % Conductance (nS)
    I_data = m3.d_orig(:, 3);    % Current (pA)
    V_data = m3.d_orig(:, 4);    % Voltage (mV)
    sims = t_data(2) - t_data(1);
    ndps = length(t_data);        % number of data points

    % Calculate the theoretical conductance trace and rescale the conductance data accordingly
    G_theo = compute_gabab_conductance(t_data, IPSC_start, gabab_amp(k), gabab_Trise(k), ...
                gabab_TfallFast(k), gabab_TfallSlow(k), gabab_w(k));
    G_data_pruned =    G_data;            
    G_data_pruned(G_data_pruned < 0) = 0;    % remove negative values from G_data first
    [G_data_rescaled, condscale2(k)] = rescale_vec(G_data_pruned, G_theo, precision);

    % Calculate current traces that might have been applied by dynamic clamp
    I_theo1 = compute_elcurr(G_data_rescaled, V_data, E_rev);
    I_theo1_ljp = compute_elcurr(G_data_rescaled, V_data, E_rev_ljp);
    I_theo2 = compute_elcurr(G_theo, V_data, E_rev);
    I_theo2_ljp = compute_elcurr(G_theo, V_data, E_rev_ljp);

    % Assuming that the dynamic clamp had a fixed voltage offset
    % estimate it by fitting to the current trace
    errorfun1 = @(x) sum((I_data - (G_data_rescaled .* (E_rev_ljp - (V_data + x) ) ) ).^2);
    [Voff1(k), sse1(k)] = fminsearch(errorfun1, Voff0);
    errorfun2 = @(x) sum((I_data - (G_theo .* (E_rev_ljp - (V_data + x) ) ) ).^2);
    [Voff2(k), sse2(k)] = fminsearch(errorfun2, Voff0);

    % Correct voltage trace by accounting for voltage offset
    V_corr1 = V_data + Voff1(k);
    V_corr2 = V_data + Voff2(k);

    % Calculate new current traces that might have been applied by dynamic clamp
    I_theo1_corr = compute_elcurr(G_data_rescaled, V_corr1, E_rev_ljp);
    I_theo2_corr = compute_elcurr(G_theo, V_corr2, E_rev_ljp);

    % Calculate tentative IPSC offset
    [gmax_data, idata] = max(G_data);
    [gmax_theo, itheo] = max(G_theo);
    ioffset_tentative(k) = (idata - itheo) * sims;    % IPSC offset (ms)

    % Plot and compare
    fprintf('Plotting %s ...\n', filebase);
    h = figure(fign);
    set(h, 'Visible', 'off')
    clf(h);

    subplot(3, 1, 1);
    plot(t_data, G_data, 'k', ...
        'DisplayName', ['rescaled data (scaling factor = ', num2str(condscale(k)), ')']); hold on;
    plot(t_data, G_theo, 'r', 'DisplayName', 'theoretical');
    plot(t_data, G_data_rescaled, 'g:', ...
        'DisplayName', ['data further rescaled (scaling factor = ', num2str(condscale2(k)), ')']);
    gub = max([gmax_data, gmax_theo]);            % rough upper bound for all conductance traces
    text(0.01, 0.9, ['Tentative IPSC offset = ', num2str(ioffset_tentative(k)), ' ms'], ...
        'Units', 'normalized', 'FontSize', 9);
    ylim([0 1.2*gub]);
    xlabel('Time (ms)')
    ylabel('Conductance (nS)')
    legend('Location', 'northeast');
    title(['Comparison of traces for ', strrep(filebase, '_', '\_')])

    subplot(3, 1, 2);
    plot(t_data, V_data, 'k', 'DisplayName', 'actual data (Vp)'); hold on;
    plot(t_data, V_corr1, 'c:', 'DisplayName', ['Vm = Vp + ', num2str(Voff1(k), 3), ' mV']);
    plot(t_data, V_corr2, 'b:', 'DisplayName', ['Vm = Vp + ', num2str(Voff2(k), 3), ' mV']);
    line([t_data(1), t_data(end)], [E_rev, E_rev], ...
        'LineStyle', '-', 'Color', 'm', 'DisplayName', 'Theoretical E\_rev');
    line([t_data(1), t_data(end)], [E_rev_ljp, E_rev_ljp], ...
        'LineStyle', '--', 'Color', 'm', 'DisplayName', 'Effective E\_rev if not ljp-corrected');
    ylim([-120 -60]);
    xlabel('Time (ms)')
    ylabel('Voltage (mV)')
    legend('Location', 'northeast');

    subplot(3, 1, 3);
    plot(t_data, I_data, 'k', ...
        'DisplayName', ['rescaled data (scaling factor = ', num2str(currscale(k)), ')']); hold on;
    plot(t_data, I_theo1, 'g', 'DisplayName', 'theoretical based on rescaled G data');
    plot(t_data, I_theo1_ljp, 'g:', 'DisplayName', 'above if not ljp-corrected');
    plot(t_data, I_theo1_corr, 'c:', 'DisplayName', 'above using Vm instead of Vp');
    plot(t_data, I_theo2, 'r', 'DisplayName', 'theoretical based on G theoretical');
    plot(t_data, I_theo2_ljp, 'r:', 'DisplayName', 'above if not ljp-corrected');
    plot(t_data, I_theo2_corr, 'b:', 'DisplayName', 'above using Vm instead of Vp');
    text(0.01, 1.1, 'I = -G*(V - E\_rev)', ...
        'Units', 'normalized');
    ylim([1.2*min([I_data; I_theo1; I_theo1_ljp; I_theo2; I_theo2_ljp]) 0]);
    xlabel('Time (ms)')
    ylabel('Current (pA)')
    legend('Location', 'northeast');

    figname = fullfile(outfigfolder, [filebase, '_trace_comp.png']);
    saveas(h, figname);

end

%% Save all variables
command = ['save(fullfile(outmatfolder, ''trace_comparison.mat''), ', ...
        '''logheader'', ''logvariables'', '];
for f = 1:numel(logvariables)
    command = [command, sprintf('''%s'', ', logvariables{f})];
end
command = [command, '''-v7.3'');'];
eval(command);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{

logheader = {'Data filename', 'Conductance scaling factor', ...
        'Current scaling factor', ...
        'Second conductance scaling factor', ... 
        'Tentative current offset (ms)', ...
        'Fitted Rs using rescaled G data (Mohm)', ...
        'Fitted Rs using G theoretical (Mohm)', ...
        'Sum of squares error using rescaled G data', ...
        'Sum of squares error using G theoretical'};

logvariables = {'fnrow', 'condscale', 'currscale', ...
        'condscale2', 'ioffset_tentative', ...
        'Rs1', 'Rs2', 'sse1', 'sse2'};
Rs1 = zeros(1, ntraces);    % fitted series resistance (MOhm) using G_data_rescaled
Rs2 = zeros(1, ntraces);    % fitted series resistance (MOhm) using G_theo 

    Rs0 = 10;        % Initial value for series resistance (MOhm)
    % Assuming the the dynamic clamp protocol corrected for series resistance, 
    % estimate the series resistance by fitting to the current trace
    errorfun1 = @(x) sum((I_data - (G_data_rescaled .* (E_rev_ljp - (V_data - I_data .* (10^-3 * x)) ) ) ).^2);
    [Rs1(k), sse1(k)] = fminsearch(errorfun1, Rs0);
    errorfun2 = @(x) sum((I_data - (G_theo .* (E_rev_ljp - (V_data - I_data .* (10^-3 * x)) ) ) ).^2);
    [Rs2(k), sse2(k)] = fminsearch(errorfun2, Rs0);

    % Correct voltage trace by accounting for series resistance
    V_corr1 = V_data - I_data * (10^-3 * Rs1(k));
    V_corr2 = V_data - I_data * (10^-3 * Rs2(k));

    plot(t_data, V_corr1, 'c:', 'DisplayName', ['Vm = Vp - I\_data * Rs (', num2str(Rs1(k), 3), ' MOhm)']);
    plot(t_data, V_corr2, 'b:', 'DisplayName', ['Vm = Vp - I\_data * Rs (', num2str(Rs2(k), 3), ' MOhm)']);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
function trace_comparison(varargin)

    legend('real data', 'theoretical')
    legend('real data', 'theoretical based on G data', ...
        'above if not ljp-corrected', ...
        'theoretical based on G theoretical', ...
        'above if not ljp-corrected')

%    h = figure('Visible', 'off');
%    close(h);

    text(0.01, 0.9, ['Second scaling factor = ', num2str(condscale2(k, 1))], ...
        'Units', 'normalized');

save(fullfile(infolder, 'trace_comparison.mat'), ...
    'condscale', 'currscale', 'ioffset_tentative', ...
    'condscale2', 'Rs1', 'Rs2', '-v7.3');

%for k = 1:ntraces            % For debug


%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
