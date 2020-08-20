function m3ha_network_plot_oscillations (infolder, oscFile)
%% Shows a swarm plot for each set of neurons (each .spi file in the infolder)
% USAGE: function m3ha_network_plot_oscillations (infolder, varargin)
% Arguments:
%    infolder   - the name of the directory containing the .syn files, 
%                   e.g. '20170317T1127_Ggaba_0.01'
%               must be a directory
%    oscFile    - the name of the file containing oscillatory periods
%                   and indices created by m3ha_network_autocorrelogram.m
%               must be a csv
%
% Requires:
%       cd/m3ha_network_autocorrelogram.m
%       cd/save_all_figtypes.m
%       /home/Matlab/Downloaded_Functions/plotSpread/plotSpread.m
% 
% Used by:
%       
%
% 2018-04-20 BT - Created
% 2019-01-24 BT - Moved to Adams_Functions from network_model

labels = {'Left Electrode', 'Middle Electrode', 'Right Electrode'}; % Table labels
colors = {'b', 'r', 'g'}; % Colors per electrode
outfolder = infolder; % Output plot in infolder

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% If csv file doesn't exist, run autocorrelogram to create it
if ~exist(fullfile(infolder, oscFile))  
    [op, oi] = autocorrelogram(infolder);
end

% Read in full table from csv
tab = readtable(fullfile(infolder, oscFile));   

% Record electrode order
headers = tab{:,1};

% Record condition order
values = tab{1:end,tab.Properties.VariableNames(2:end-1)};

% First half is period
oscillatoryPeriodMat = values(:,1:size(values,2)/2)';

% Second half is index
oscillatoryIndexMat = values(:,size(values,2)/2+1:end)';

% Remove gIncr
oscillatoryPeriodMat = oscillatoryPeriodMat(1:4,:);
oscillatoryIndexMat = oscillatoryIndexMat(1:4,:);

% Group same electrode values together
categ = ones(size(oscillatoryPeriodMat));   
for i = 2:size(oscillatoryPeriodMat,1)
    categ(i,:) = i;
end

% Plot oscillatory period spread
h = figure(60000);
h.Visible = 'off';
clf(h);
rowheader = strrep(headers, '_', '\_');
plotSpread(oscillatoryPeriodMat,'categoryIdx', categ, 'categoryColors', colors, 'xNames',rowheader,'ylabel','Oscillatory Period');
title('Oscillatory Period Spread');
legend(labels, 'Location', 'northeast');
% Save figure
figname = fullfile(outfolder, 'period_spread.png');
save_all_figtypes(h, figname, 'png');
close(h);

% Plot oscillatory index spread
h = figure(70000);
h.Visible = 'off';
clf(h);
plotSpread(oscillatoryIndexMat,'categoryIdx', categ, 'categoryColors', colors, 'xNames',rowheader,'ylabel','Oscillatory Index');
title('Oscillatory Index Spread');
legend(labels, 'Location', 'northeast');
% Save figure
figname = fullfile(outfolder, 'index_spread.png');
save_all_figtypes(h, figname, 'png');
close(h);
end
