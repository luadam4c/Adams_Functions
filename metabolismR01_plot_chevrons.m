% metabolismR01_plot_chevrons.m
%% Plot Chevron plots and computes p values for the metabolism R01 grant renewal 20200225
%
% Requires:
%       cd/all_files.m
%       cd/plot_chevron.m

% Paths
grantDir = '/media/shareX/2020marchR01';
chevronDir = fullfile(grantDir, 'chevronStats');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Find all Chevron spreadsheet files
[~, chevronPaths] = all_files('Directory', chevronDir, ...
                                'Keyword', 'chevron', 'Extension', 'csv');

% Plot all Chevron plots
cellfun(@(file) plot_chevron(file), chevronPaths);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%