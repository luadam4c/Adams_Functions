function m3ha_plot_correlations (fitmode, infolder, outfolder)
%% Plot Correlation diagrams for data that will be used for fitting
% Usage: m3ha_plot_correlations (fitmode, infolder, outfolder)
% Arguments: 
%        fitmode        - 0 - all data
%                    - 1 - all of g incr = 100%, 200%, 400%
%                    - 2 - all of g incr = 100%, 200%, 400% 
%                    but exclude cell-pharm-g_incr sets containing problematic sweeps
%        infolder    - (opt) the directory that contains the matfile to read
%                    must be a directory
%                    default == //media/adamX/m3ha/data_dclamp/take4/
%        outfolder    - (opt) the directory to output correlation diagrams
%                    (different subdirectories will be created for each fitmode)
%                    must be a directory
%                    default == //media/adamX/m3ha/data_dclamp/take4/
% 
% Requires:    
%       "infolder"/dclampdatalog_take4.mat
%       cd/find_in_strings.m
%       cd/m3ha_find_ind_to_fit.m
%       cd/m3ha_specs_for_fitmode.m
%
% Used by:    
%       cd/m3ha_parse_dclamp_data.m
%

% File History:
% 2016-10-13 - Created by BT, adapted from PlotHistogramsRefineThreshold.m
% 2016-10-17 - BT - Changed scatterplot to color points by peak class and added legend.
% 2016-10-18 - BT - Changed color scheme of existing figures. 
% 2016-10-20 - BT - Changed importing of vectors to all vectors within mat file
% 2016-10-20 - BT - Made correlations combinatorial and to generate in parfor
% 2016-10-20 - BT - Added optional parameters for infolder and outfolder
% 2016-10-20 - AL - Made variables consistent with PlotHistogramsRefineThreshold.m
% 2016-10-20 - AL - Fixed error when ioffset_old was not skipped during plotting
% 2016-10-20 - AL - Preallocate memory
% 2016-10-21 - AL - Added classExists to account for the case when peakclass does not include all classes
% 2016-10-31 - AL - Placed suffix & titleMod into specs_for_fitmode.m
% 2016-11-03 - BT - Added label for correlation coefficient
% 2016-11-08 - BT - Copy correlation plots with |correlation coefficient| 
%                       above threshold (corr_thr) in a subdirectory called “interesting”
% 2016-11-09 - AL - Now saves correlation matrix to a mat file in outfolder
% 2017-11-09 - AL - Updated peakclassLabels
% 2018-02-11 - AL - Changed variable names to camel form
% 2018-02-11 - AL - Streamlined some parts of the code
% 2018-02-11 - AL - Annotated code
% 2018-02-12 - AL - Changed parfor to iterate over iAll
% 2018-02-12 - AL - Changed default for correlation coefficients from 0 to 1
% 2018-02-12 - AL - Now checks whether the vector lengths are identical 
%                       (nTraces & nEntries)
% 2018-02-12 - AL - Now checks whether peak classes exist 
%                       after restricting vectors

%% Flags
ltsOnly = 0;        % whether to include only the traces with LTSs

%% Parameters
corr_thr = 0.5;        % threshold for absolute value of correlation coefficient 
                    %   that is deemed "interesting"

%% Peak Classification (corresponds to peakclass == 1 ~ 9):
peakclassLabels = {'Not LTS by prominence', ...
                    'Not LTS by narrowness', ...
                    'Not LTS by shape', ...
                    'Not LTS by overrule', ...
                    'LTS with no burst by overrule', ...
                    'LTS with no burst; contentious', ...
                    'LTS with burst; contentious', ...
                    'LTS with no burst; definite', ...
                    'LTS with burst; definite'};

%% Peak colors
peakclassColor = colormap(parula(numel(peakclassLabels)));

%% Specify which matfile to use; assumed to be in infolder
filetouse = 'dclampdatalog_take4.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check arguments
if nargin < 1
    error('A fitmode is required, type ''help m3ha_plot_correlations'' for usage');
elseif isempty(fitmode) || ~isnumeric(fitmode) || ...
        ~(fitmode == 0 || fitmode == 1 || fitmode == 2)
    error('fitmode out of range!, type ''help m3ha_plot_correlations'' for usage');
elseif nargin >= 2 && ~isdir(infolder)
    error('infolder must be a directory!');
elseif nargin >= 3 && ~isdir(outfolder)
    error('outfolder must be a directory!');
end

%% Set defaults for optional arguments
if nargin < 2
    infolder = '/media/adamX/m3ha/data_dclamp/take4/';
end
if nargin < 3
    outfolder = '/media/adamX/m3ha/data_dclamp/take4/';
end

%% Find path of matfile to use
fullmatfilepath = fullfile(infolder, filetouse);
m = matfile(fullmatfilepath);

%% Print user specifications
fprintf('Using fit mode == %d ... \n', fitmode);
fprintf('Using matfile == %s ... \n', fullmatfilepath);

%% Set suffix according to fitmode
suffix = m3ha_specs_for_fitmode(fitmode);

%% Create output folders for saving files
outfolderCorr = fullfile(outfolder, ['correlations', suffix]);
if exist(outfolderCorr, 'dir') ~= 7
    mkdir(outfolderCorr);
    fprintf('Made directory %s \n', outfolderCorr);
end
subfolderCorr = fullfile(outfolderCorr, 'interesting');
if exist(subfolderCorr, 'dir') ~= 7
    mkdir(subfolderCorr);
    fprintf('Made directory %s \n', subfolderCorr);
end

%% Extract information for the variables (sweep information) to correlate
varLabels = m.logheader(1, 2:end);  % variable labels excludes filenames
varNames = m.logvariables(1, 2:end);% variable names corresponding to each label
nVars = numel(varNames);            % number of variables to correlate
allVecs = cell(1, nVars);           % stores a vector for each variable
nTraces = [];                       % stores the number of traces used
for iVar = 1:nVars                     % for each variable
    % Read in the vector for this variable from matfile
    allVecs{iVar} = m.(varNames{iVar});

    % Check if the length is identical to previous vectors
    nEntries= length(allVecs{iVar});
    if isempty(nTraces)
        nTraces = nEntries;
    elseif nEntries ~= nTraces
        error(['The variable %s has a different number of entries (%d)', ...
                'than the previous vectors (%d)!'], ...
                varNames{iVar}, nEntries, nTraces);
    end
end


%% Restrict data to fit if fitmode > 0
if fitmode > 0
    % Find indices of fnrow in dclampdatalog_take4.mat that will be used for fitting
    fnrowOld = m.fnrow;             % file names for each sweep
    cellidrowOld = m.cellidrow;     % Cell ID # for each sweep
    prowOld = m.prow;               % Pharm condition for each sweep
    growOld = m.grow;               % G incr for each sweep
    indtofit = m3ha_find_ind_to_fit(fnrowOld, cellidrowOld, prowOld, ...
                                growOld, fitmode, infolder);

    % Temporarily take out the traces with spikelatency < -10 ms
    spikelatency = m.spikelatency;
    indNotWeird = find(isnan(spikelatency) | spikelatency >= -10);

    % Restrict to those traces of interest for fitting
    indToRestrict = intersect(indtofit, indNotWeird);
    for k = 1:nVars
        allVecs{k} = allVecs{k}(indToRestrict);
    end
end

%% Restrict vectors if limiting to traces with an LTS
if ltsOnly
    % TODO: Always plot these but with a different filename
    % Find indices of those traces with an LTS
    ltspeaktime = allVecs{find_in_strings('ltspeaktime', varNames)};
    indGoodLts = find(~isnan(ltspeaktime));

    % Restrict to those traces with an LTS
    for k = 1:nVars
        allVecs{k} = allVecs{k}(indGoodLts);
    end
end

%{
%% Plot correlations
% Extract needed vectors and information
peakclass = allVecs{find_in_strings('peakclass', varNames)};
nPeakclass = length(peakclassLabels);

% Compute and plot correlations
corrMat = ones(nVars, nVars);           % stores correlation coefficients
parfor iAll = 1:nVars^2                 % for each pair of variables
    % Get the indices of each variable
    [iX, iY] = ind2sub([nVars, nVars], iAll);

    % Don't autocorrelate and skip ioffset_old
    if iY == iX || strcmp('ioffset_old', varNames{iY}) ...
        || strcmp('ioffset_old', varNames{iX})
        % Do nothing
    else
        % Create figure with a unique number and don't show it
        h = figure(iX*nVars+iY);
        h.Visible = 'off';

        % Define the vectors to correlate
        xVec = allVecs{iX};
        yVec = allVecs{iY};

        % Remove traces from both vectors if at least one of them 
        %   has NaN as a value for that trace
        neitherNaN = ~isnan(xVec) & ~isnan(yVec);
        xVecNoNaN = xVec(neitherNaN);
        yVecNoNaN = yVec(neitherNaN);
        peakclassNoNaN = peakclass(neitherNaN);

        % Find all the peak classes that remains
        classRemains = unique(peakclassNoNaN);
        thisPeakclassLabels = peakclassLabels(classRemains);
        thisPeakclassColor = peakclassColor(classRemains, :);

        % Plot a scatter plot with each peak class in a different color
        gscatter(xVecNoNaN, yVecNoNaN, peakclassNoNaN, thisPeakclassColor, 'o');
        legend(thisPeakclassLabels, 'Location', 'eastoutside');

        % Compute and show the correlation coefficient
        correlation = corr2(xVecNoNaN, yVecNoNaN);
        corrMat(iAll) = correlation;
        text(1.1, 0.95, ...
            ['Correlation coefficient: ', num2str(correlation, 3)], ...
            'Units', 'normalized'); 

        % Titles and axis labels
        title(strjoin({'Correlation of', varLabels{iY}, ...
                        'vs.', varLabels{iX}, titleMod}));
        xlabel(varLabels{iX});
        ylabel(varLabels{iY});

        % Save the figure
        figName = [varNames{iY}, '_vs_', varNames{iX}, suffix, '.png'];
        figFile = fullfile(outfolderCorr, figName);
        saveas(h, figFile);

        % Copy the figure file to a subdirectory if the absolute value 
        %   of the correlation coefficient is above threshold
        if abs(correlation) > corr_thr
            subFigFile = fullfile(subfolderCorr, figName);
            copyfile(figFile, subFigFile);
        end
    end
    close all force hidden
end

%% Save correlation matrix to a mat file
corrMat_file = fullfile(outfolder, ['correlation_matrix', suffix, '.mat']);
save(corrMat_file, 'varLabels', 'varNames', 'corrMat', '-v7.3');
%}

%% Perform multi-linear regression with an interaction term 
%   on maxslopval and maxslopetime
% Find indices of variables to regress on
idxX1 = find_in_strings('maxslopetime', varNames);
idxX2 = find_in_strings('maxslopeval', varNames);

% Extract vectors and construct X matrix
x1 = allVecs{idxX1}(:);
x2 = allVecs{idxX2}(:);
X = [ones(size(x1)), x1, x2, x1.*x2];

% Find indices of output variables
idxY{1} = find_in_strings('spikelatency', varNames);
idxY{2} = find_in_strings('spikesperpeak', varNames);
idxY{3} = find_in_strings('spikefrequency', varNames);

% Compute the regression coefficients y = b0 + b1 * x1 + b2 * x2 + b3 * x1 * x2
y = cell(3);            % stores output vectors
b = cell(3);            % stores regression coefficients
bInt = cell(3);         % stores confidence intervals of regression coefficients
for iOut = 1:3          % for each output variable
    % Extract output vector
    y{iOut} = allVecs{idxY{iOut}}(:);

    % Perform regression
    [b{iOut}, bInt{iOut}] = regress(y{iOut}, X);
end

% Compare the confidence intervals of each regression coefficient
% TODO

% Create a meshgrid of maxslopetime and maxslopeval values
x1fit = linspace(min(x1), max(x1), 100);
x2fit = linspace(min(x2), max(x2), 100);
[X1FIT, X2FIT] = meshgrid(x1fit, x2fit);

% Plot the data and the model together
YFIT = cell(3);         % stores output variable values
parfor iOut = 1:3       % for each output variable
    h = figure('Visible', 'off');
    scatter3(x1, x2, y{iOut}, 'filled');
    hold on
    YFIT{iOut} = b{iOut}(1) + b{iOut}(2)*X1FIT + ...
                    b{iOut}(3)*X2FIT + b{iOut}(4)*X1FIT.*X2FIT;
    mesh(X1FIT, X2FIT, YFIT{iOut})
    xlabel(varLabels{idxX1});
    ylabel(varLabels{idxX2});
    zlabel(varLabels{idxY{iOut}})
    figName = [varNames{idxY{iOut}}, '_vs_', ...
                varNames{idxX1}, '_', varNames{idxX2}, suffix, '.png'];
    figFile = fullfile(outfolderCorr, figName);
    saveas(h, figFile);
    close(h);
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
%% OLD CODE
peakclass = m.peakclass;                % peakclass for color grouping
    %% Update peakclass for color grouping
    peakclass = peakclass(indtofit);

%% Generate file names
counter = 1;
filenames = {};
for g = 1:nVars-1
    for n = g+1:nVars
        if g ~= n && strcmp('ioffset_old',varNames(g)) == 0 && strcmp('ioffset_old',varNames(n)) == 0
            filenames(g*nVars+n) = strcat(varNames(g),'_',varNames(n), suffix, '.png');
        end
    end
end

%% Plot correlations
counter = 1;
for k = 1:nVars-1
    for j = k+1:nVars
        if j ~= k
            h = figure(k*nVars+j);
            h.Visible = 'off';
            gscatter(allVecs{k}, allVecs{j}, peakclass, peakclassColor, 'o');
            legend(peakclassLabels, 'Location', 'eastoutside');
            title(strjoin(['Correlation of', varLabels(k), 'vs.', varLabels(j), titleMod]));
            xlabel(varLabels(k));
            ylabel(varLabels(j));            
            figname = fullfile(outfolderCorr, filenames{k*nVars+j});
            saveas(h, figname);            
            close(h);
        end
    end
end

            legend(peakclassLabels, 'Location', 'eastoutside');
ltspeaktime_ind = find(strcmp(varNames, 'ltspeaktime'));

%            close(h);

    if sum(peakclass == cl) > 0

% matched to existing scheme of histg.m
peakclassColor = [0.2118, 0.1647, 0.5294;
           0.0157, 0.4235, 0.8784;
           0.0235, 0.6078, 0.8118;
           0.2157, 0.7216, 0.6118;
           0.6510, 0.7490, 0.4196;
           1.0000, 0.7490, 0.2392;
           1.0000, 1.0000, 0.0000];

            filenames{g*nVars+n} = ...
                [varNames{g}, '_', varNames{n}, suffix, '.png'];
            title(strjoin({'Correlation of', varLabels{k}, 'vs.', varLabels{j}, titleMod}));

            noNanCorrelK = [];
            noNanCorrelJ = [];
            if ltsOnly
                % TODO: Was this working? Why not plot it with a different filename?
                goodlts_ind = find(~isnan(allVecs{ltspeaktime_ind}));
                gscatter(xVec(goodlts_ind), yVec(goodlts_ind), ...
                        peakclass(goodlts_ind), peakclassColor(4:7,:), 'o');
                classExists_sub = classExists(classExists >= 4);
                legend(peakclassLabels(classExists_sub), 'Location', 'eastoutside');
                % Removes NaN values and corresponding values from both arrays (needed for corr2) 
                noNanCorrelK = xVec(goodlts_ind);
                noNanCorrelJ = yVec(goodlts_ind);
                noNanIndicesK = find(~isnan(noNanCorrelK));
                noNanCorrelK = noNanCorrelK(noNanIndicesK); 
                noNanCorrelJ = noNanCorrelJ(noNanIndicesK);
                noNanIndicesJ = find(~isnan(noNanCorrelJ));
                noNanCorrelJ = noNanCorrelJ(noNanIndicesJ);
                noNanCorrelK = noNanCorrelK(noNanIndicesJ);
            else
                gscatter(xVec, yVec, peakclass, ...
                            peakclassColor(classExists, :), 'o');
                legend(peakclassLabels(classExists), 'Location', 'eastoutside');
                % Removes NaN values and corresponding values from both arrays (needed for corr2) 
                noNanIndicesK = find(~isnan(xVec));
                noNanCorrelK = xVec(noNanIndicesK); 
                noNanCorrelJ = yVec(noNanIndicesK);
                noNanIndicesJ = find(~isnan(noNanCorrelJ));
                noNanCorrelJ = noNanCorrelJ(noNanIndicesJ);
                noNanCorrelK = noNanCorrelK(noNanIndicesJ);
            end

        if iY >= iX || strcmp('ioffset_old', varNames{iX}) ...
            || strcmp('ioffset_old', varNames{iY})
            % Don't plot both ways and skip ioffset_old

%% Generate file names
filenames = cell(1, nVars^2);            % preallocate memory
parfor iX = 1:nVars
    for iY = 1:nVars
        if iY ~= iX || strcmp('ioffset_old', varNames{iX}) ...
            || strcmp('ioffset_old', varNames{iY})
            % Don't autocorrelate and skip ioffset_old
        else
            filenames{iX*nVars+iY} = [varNames{iY}, '_vs_', ...
                                        varNames{iX}, suffix, '.png'];
        end
    end
end
figname = fullfile(outfolderCorr, filenames{iX*nVars+iY});

for iX = 1:nVars
    for iY = 1:nVars
    end
end

corrMat = zeros(nVars, nVars);          % stores correlation coefficients

    ltspeaktime = allVecs{find_in_strings('ltspeaktime', varNames)};
    indGoodLts = find(~isnan(ltspeaktime));
    xVec = xVec(indGoodLts);
    yVec = yVec(indGoodLts);

        classExists = [];                       % peak classes that exist
        for cl = 1:nPeakclass                   % for every peak class
            if any(peakclassNoNaN == cl)        % if it exists in peakclassNoNaN
                % Add the class to classExists
                classExists = [classExists, cl];
            end
        end

    switch iOut
    case 1
        zlabel('Spike Latency (ms)');
    case 2
        zlabel('Spikes Per Peak');
    case 3
        zlabel('Spike Frequency (Hz)');
    otherwise
    end

%}
