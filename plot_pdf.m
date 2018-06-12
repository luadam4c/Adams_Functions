function [xValues, pdfValues] = plot_pdf (X, pdfModel, varargin)
%% Plots pdf fit of data X and return vectors for the plots
% Usage: [xValues, pdfValues] = plot_pdf (X, pdfModel, varargin)
% Example(s):
%       plot_pdf(data, @(x) normpdf(x));
%       model = fitdist(data, 'Kernel');
%       [xValues, pdfValues] = plot_pdf(data, @(x) pdf(model, x), ...
%                                           'NPointsToPlot', 1000, ...
%                                           'PlotFlag', false);
%       Value = [2; 3; 4];
%       Color = {'Red'; 'Green'; 'Purple'};
%       LineStyle = {'--'; '--'; '--'};
%       Label = {'Peak 1'; 'Peak 2'; 'Threshold'};
%       linesToPlot = table(Value, Color, LineStyle, Label);
%       plot_pdf(data, @(x) pdf(model, x), 'LinesToPlot', linesToPlot);
%
% Outputs:
%       xValues     - x values for the pdf plot
%                   specified as a numeric column vector
%       pdfValues   - y values for the pdf plot
%                   specified as a numeric column vector
% Side Effects:
%       Plots a scaled probability density function
% Arguments:    
%       X           - data to distribute among bins
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       pdfModel    - pdf function for a fitted model
%                   must be a function handle
%       varargin    - 'NPointsToPlot': number of points for plotting the pdf
%                   must be a positive integer scalar
%                   default == 500
%                   - 'PlotFlag': whether to plot the pdf
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'LinesToPlot': table with fields 
%                   must be a table with the following columns:
%                       'Value', 'Color', 'LineStyle', 'Label' 
%                   default == []
%
% Requires:
%       /home/Matlab/Adams_Functions/histproperties.m
%       /home/Matlab/Downloaded_Functions/rgb.m
%
% Used by:    
%       /home/Matlab/Adams_Functions/fit_kernel.m
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%
% File History:
% 2018-06-12 Adapted from fit_IEI.m 
% 

%% Default values for optional arguments
nPointsToPlotDefault = 500;     % default number of points for plotting the pdf
plotFlagDefault = true;         % whether to plot a pdf by default
linesToPlotDefault = [];        % do not plot lines by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions across servers
if ~isdeployed
    if exist(fullfile(pwd, 'Downloaded_Functions'), 'dir') == 7
        functionsdirectory = pwd;
    elseif exist('/home/Matlab/', 'dir') == 7
        functionsDirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsDirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsDirectory does not exist!');
    end
    addpath(fullfile(functionsDirectory, 'Downloaded_Functions')); 
                                            % for rgb.m
end

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'X', ...                    % data to distribute among bins
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'nonempty'}));
addRequired(iP, 'pdf');

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NPointsToPlot', nPointsToPlotDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'PlotFlag', plotFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'LinesToPlot', linesToPlotDefault);
    % TODO: istable?

% Read from the Input Parser
parse(iP, X, pdfModel, varargin{:});
nPointsToPlot = iP.Results.NPointsToPlot;
plotFlag = iP.Results.PlotFlag;
linesToPlot = iP.Results.LinesToPlot;

% Compute the area and edges for the histogram
[harea, edges] = histproperties(X);

% Construct x values for the pdf plots
xValues = linspace(edges(1), edges(end), nPointsToPlot)';

% Construct y values for the pdf plots
pdfValues = harea * pdfModel(xValues);

%% Plot the pdf
if plotFlag
    plot(xValues, pdfValues, 'DisplayName', 'fit');
    if ~isempty(linesToPlot)
        hold on
        ylimits = get(gca, 'YLim');
        for iRow = 1:height(linesToPlot)
            line(linesToPlot.Value(iRow) * ones(1, 2), ylimits, ...
                 'Color', rgb(linesToPlot.Color{iRow}), ...
                 'LineStyle', linesToPlot.LineStyle{iRow}, ...
                 'DisplayName', linesToPlot.Label{iRow});
        end
        hold off
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}
