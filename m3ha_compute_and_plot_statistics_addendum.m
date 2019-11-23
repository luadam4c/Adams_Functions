% m3ha_compute_and_plot_statistics_addendum.m
%% Analyzes additional aspects of dclampdatalog_take4.mat
%
% 2016-11-10 Created
% 2018-02-04 Added calculation of spikelatency
% 2019-11-22 Renamed from analyze_dclampdatalog_take4.m

%% Specify which matfile to use; assumed to be in infolder
filetouse = 'dclampdatalog_take4.mat';

%% New variables to add
newheader = {'Spike latency (ms)'};
newvariables = {'spikelatency'};

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%infolder = '/media/adamX/m3ha/data_dclamp/take4/debug/';
infolder = '/media/adamX/m3ha/data_dclamp/take4/';

%% Find path of matfile to use
fullmatfilepath = fullfile(infolder, filetouse);

%% Add variables to matfile
% Open matfile for writing
m = matfile(fullmatfilepath, 'Writable', true);

% Compute “spike latency”: time from maxslopetime to firstspiketime
%   and store in matfile
maxslopetime = m.maxslopetime;
firstspiketime = m.firstspiketime;
spikelatency = firstspiketime - maxslopetime;
m.spikelatency = spikelatency;

% If new variables not already saved, save it
logheader = m.logheader;
logvariables = m.logvariables;
if ~any(strcmp(logvariables, 'spikelatency'))
    m.logheader = [logheader, newheader];
    m.logvariables = [logvariables, newvariables];
end

%% Perform principal component analysis
% Open matfile for reading
m = matfile(fullmatfilepath);

% Extract variables from matfile
fnrow = m.fnrow;                    % data file name
ntraces = numel(fnrow);             % number of traces

var_labels = m.logheader(1, 2:end);         % labels are for all sweep info except data filename
var_names = m.logvariables(1, 2:end);       % variable names corresponding to each label
nvars = numel(var_names);                   % number of variables
bursttime = m.bursttime;
goodb_ind = find(bursttime > 0);
ltspeaktime = m.ltspeaktime;
goodlts_ind = find(ltspeaktime > 0);
goodb_ind = find(bursttime > 0);
ngoodlts = length(goodlts_ind);
ngoodb = length(goodb_ind);
all_data = zeros(ntraces, nvars);           % Puts all sweep info in a single matrix
all_data_goodlts_norm = zeros(ngoodlts, nvars - 4);    % Puts all sweep info with LTSs in a single matrix
all_data_goodb_norm = zeros(ngoodb, nvars);    % Puts all sweep info with bursts in a single matrix
ct = 0;
for k = 1:nvars
    var_vector = m.(var_names{k});    % Read in the sweep information from matfile
    all_data(:, k) = var_vector;
    if k ~= 13 && k ~= 18 && k ~= 28 && k ~= 29
        ct = ct + 1;
        all_data_goodlts_norm(:, ct) = var_vector(goodlts_ind)/norm(var_vector(goodlts_ind));
        lts_var_names{ct} = var_names{k};
        lts_var_labels{ct} = var_labels{k};
    end
    all_data_goodb_norm(:, k) = var_vector(goodb_ind)/norm(var_vector(goodb_ind));
end
zsc_lts = zeros(size(all_data_goodlts_norm, 1), size(all_data_goodlts_norm, 2));
zsc_b = zeros(size(all_data_goodb_norm, 1), size(all_data_goodb_norm, 2));
for g = 1:size(zsc_b, 2)
    zsc_b(:,g) = zscore(all_data_goodb_norm(:,g));
    if g <= size(zsc_b, 2)    
        zsc_lts(:,g) = zscore(all_data_goodlts_norm(:,g));
    end
end

% Perform principal component analysis
[coeffs_lts, scores_lts, latent_lts, tsquared_lts, explained_lts, mu_lts] = pca(zsc_lts);
fprintf('Each principal component explains this much percent of variation in LTS data:\n');
disp(explained_lts);
for n = 1:4
    [B, I1] = sort(abs(coeffs_lts(:, n)), 'descend');
    fprintf(['Important variables in component #', num2str(n), ' are: \n']);
    disp(var_names(I1(B > 0.1)));
end
%biplot(coeffs_lts(:,1:3),'scores',scores_lts(:,1:3),'varlabels',lts_var_names);

[coeffs_b, scores_b, latent_b, tsquared_b, explained_b, mu_b] = pca(zsc_b);
fprintf('Each principal component explains this much percent of variation in burst data:\n');
disp(explained_b);
for n = 1:4
    [B, I1] = sort(abs(coeffs_b(:, n)), 'descend');
    fprintf(['Important variables in component #', num2str(n), ' are: \n']);
    disp(var_names(I1(B > 0.1)));
end
biplot(coeffs_b(:,1:3),'scores',scores_b(:,1:3),'varlabels',var_labels);

%% Find possible burst thresholds
narrowpeak2ndder = m.narrowpeak2ndder;
burst_thr = max(narrowpeak2ndder(goodb_ind));
fprintf('Possible 2nd derivative threshold for bursts is %g\n', burst_thr);

%{
%% OLD CODE

%% Load statistics
load('take4/dclampdatalog_take4.mat')

%}
