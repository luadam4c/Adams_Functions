function [h, g] = plot_grouped_scatter(figname, X, Y, grouping, grouping_labels, xLabel, xUnits, yLabel, yUnits, titleStr, varargin)
%% Plot and save a grouped scatter plot with 95% confidence ellipses
% Usage: [h, g] = plot_grouped_scatter(figname, X, Y, grouping, grouping_labels, xLabel, xUnits, yLabel, yUnits, titleStr, varargin)
%
% Requires:
%       /home/Matlab/Adams_Functions/plot_ellipse.m
%
% Used by:
%		/media/adamX/Paula_IEIs/paula_iei4.m
%
% 2017-12-13 - Modified from plot_grouped_histogram.m

%% TODO: make the following optional arguments with given default
ellipseNPoints = 1000;              % 1000 points
ellipseLineStyle = '-';             % dashed line
ellipseLineWidth = 1;

%% Default values for optional arguments
plotEllipseDefault = true;          % whether to plot ellipses by default
percentCLDefault = 95;              % default confidence level (%)
xScaleDefault = 'linear';
yScaleDefault = 'linear';
xLimitsDefault = [];
yLimitsDefault = [];
markerSizeDefault = 6;
markerTypeDefault = 'o';
markerLineWidthDefault = 0.5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'plot_grouped_scatter';

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'PlotEllipse', plotEllipseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'PercentCL', percentCLDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));
addParameter(iP, 'XScale', xScaleDefault, ...
    @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(iP, 'YScale', yScaleDefault, ...
    @(x) any(validatestring(x, {'linear', 'log'})));
addParameter(iP, 'XLimits', xLimitsDefault, ...
    @(x) validateattributes(x, {'numeric', 'categorical', ...
        'datetime', 'duration'}, {'vector', 'numel', 2}));
addParameter(iP, 'YLimits', yLimitsDefault, ...
    @(x) validateattributes(x, {'numeric', 'categorical', ...
        'datetime', 'duration'}, {'vector', 'numel', 2}));
addParameter(iP, 'MarkerSize', markerSizeDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'MarkerType', markerTypeDefault, ...
    @(x) any(validatestring(x, {'o', '+', '*', '.', 'x', 's', 'd', ...
                                '^', 'v', '>', '<', 'p', 'h'})));
addParameter(iP, 'MarkerLineWidth', markerLineWidthDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive'}));

% Read from the Input Parser
parse(iP, varargin{:});
plotEllipse = iP.Results.PlotEllipse;
percentCL = iP.Results.PercentCL;
xScale = iP.Results.XScale;
yScale = iP.Results.YScale;
xLimits = iP.Results.XLimits;
yLimits = iP.Results.YLimits;
markerSize = iP.Results.MarkerSize;
markerType = iP.Results.MarkerType;
markerLineWidth = iP.Results.MarkerLineWidth;

%% Determine colors
nGroups = max(grouping);
cm = colormap(parula(nGroups));

%% Fit each group with a bivariate Gaussian
if plotEllipse
    % Compute the critical value for the chi-squared distribution
    %   to have a cumulative probability of percentCL % confidence level
    criticalValue = chi2inv(percentCL/100, 2);

    % Transform data to column vectors
    X = X(:);
    Y = Y(:);

    % Transform data to scale used for plotting and combine into a matrix
    if strcmp(xScale, 'log')
        Xscaled = log10(X);
    else
        Xscaled = X;
    end
    if strcmp(yScale, 'log')
        Yscaled = log10(Y);
    else
        Yscaled = Y;
    end
    data = [Xscaled, Yscaled];

    % Use the sample means and sample covariance matrices, 
    %   the maximum likelihood estimates of the corresponding parameters 
    %   in a bivariate Gaussian fit, to compute the center, 
    %   half lengths, and angles of the percentCL % confidence ellipse
    mus = cell(nGroups, 1);         % stores mean vectors for each group
    sigmas = cell(nGroups, 1);      % stores covariance matrices for each group
    eigenvalues = cell(nGroups, 1); % stores eigenvalues for each group
    eigenvectors = cell(nGroups, 1);% stores eigenvectors for each group
    centers = cell(nGroups, 1);     % stores centers of ellipses for each group
    halflengths = cell(nGroups, 1); % stores half lengths of ellipses for each group
    theta0s = zeros(nGroups, 1);    % stores angles of rotation for each group
    xValues = cell(nGroups, 1);     % stores x values for ellipse
    yValues = cell(nGroups, 1);     % stores y values for ellipse
    for iGroup = 1:nGroups
        % Compute the sample mean for non-NaN data
        mus{iGroup} = nanmean(data(grouping == iGroup, :)); % a row vector

        % Compute the sample covariance for non-NaN data
        sigmas{iGroup} = nancov(data(grouping == iGroup, :));

        % Compute the eigenvalues of the covariance matrix
        eigenvalues{iGroup} = eig(sigmas{iGroup});

        if any(eigenvalues{iGroup} <= 0)
            fprintf('There seems to be a direct correlation!\n');
            fprintf('No ellipse will be plotted for Group #%d!\n', iGroup);
            fprintf('Covariances: \n');
            disp(sigmas{iGroup});
            fprintf('Eigenvalues: \n');
            disp(eigenvalues{iGroup});
            fprintf('\n');
        end

        % Find the half lengths and angle of rotation for the ellipse
        %   Note: if an eigenvalue is zero, there is direct correlation and
        %           no ellipse will be plotted
        if length(eigenvalues{iGroup}) >= 2 && ...
            all(eigenvalues{iGroup} > 0)

            % Compute half lengths for the 95% confidence ellipse
            halflengths{iGroup} = sqrt(eigenvalues{iGroup} .* criticalValue);

            % Compute the eigenvectors of the covariance
            [eigenvectors{iGroup}, ~] = eig(sigmas{iGroup});

            % Compute the angle of rotation in radians
            %   This is the angle between the x axis and the first eigenvector
            theta0s(iGroup) = atan(eigenvectors{iGroup}(2, 1)/...
                                        eigenvectors{iGroup}(1, 1));

        end

        if length(mus{iGroup}) >= 2 && ~isempty(halflengths{iGroup}) && ...
            ~isempty(theta0s(iGroup))
            % Obtain the x and y values of the ellipse on the scaled plot
            [~, xPlot, yPlot] = ...
                plot_ellipse(mus{iGroup}, halflengths{iGroup}, ...
                            theta0s(iGroup), 'NPoints', ellipseNPoints, ...
                            'ToPlot', false);
            if strcmp(xScale, 'log')
                xValues{iGroup} = 10.^(xPlot);
            else
                xValues{iGroup} = xPlot;
            end
            if strcmp(yScale, 'log')
                yValues{iGroup} = 10.^(yPlot);
            else
                yValues{iGroup} = yPlot;
            end
        end
    end
end

%% Plot and save scatter plot
h = figure('Visible', 'off');
clf(h);
g = gscatter(X, Y, grouping, cm, markerType);
set(g, 'MarkerSize', markerSize);
set(g, 'LineWidth', markerLineWidth);
if plotEllipse
    % Plot a 95% confidence ellipse for each group
    hold on;
    for iGroup = 1:nGroups
        if ~isempty(xValues{iGroup}) && ~isempty(yValues{iGroup})
                plot(xValues{iGroup}, yValues{iGroup}, ...
                    'Color', cm(iGroup, :), ...
                    'LineStyle', ellipseLineStyle, ...
                    'LineWidth', ellipseLineWidth);
        end
    end    
    hold off;
end
set(gca, 'XScale', xScale, 'YScale', yScale);
if ~isempty(xLimits)
    xlim(xLimits);
end
if ~isempty(yLimits)
    ylim(yLimits);
end
legend(grouping_labels, 'Interpreter', 'none', 'location', 'eastoutside');    
if ~isempty(xUnits)
    xlabel([xLabel, ' (', xUnits, ')']);
else
    xlabel(xLabel);
end
if ~isempty(yUnits)
    ylabel([yLabel, ' (', yUnits, ')']);
else
    ylabel(yLabel);
end
title(titleStr, 'Interpreter', 'none');
saveas(h, figname, 'png');
close(h);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Compute the center position of the ellipse
if length(mus{iGroup}) >= 2 && ~isnan(mus{iGroup}(1))
    if strcmp(xScale, 'log')
        centers{iGroup}(1) = 10^(mus{iGroup}(1));
    else
        centers{iGroup}(1) = mus{iGroup}(1);
    end
end
if length(mus{iGroup}) >= 2 && ~isnan(mus{iGroup}(2))
    if strcmp(yScale, 'log')
        centers{iGroup}(2) = 10^(mus{iGroup}(2));
    else
        centers{iGroup}(2) = mus{iGroup}(2);
    end
end

%}
