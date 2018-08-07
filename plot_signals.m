function plot_signals(timeVec, data, varargin)

%% Default values for optional arguments
titleDefault = 'Some data';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'timeVec', ... 
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'data', ... 
    @(x) validateattributes(x, {'numeric'}, {'nonempty'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Title', titleDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, timeVec, data, varargin{:});
figTitle = iP.Results.Title;

% Get the number of channels in the data
nChannels = size(data, 2);

% Decide on the colormap
cm = colormap(jet(nChannels));

% Find the minimum time
minTime = min(timeVec);

% Find the maximum time
maxTime = max(timeVec);

% Plot each channel as a different subplot
for iPlot = 1:nChannels
    % Create a subplot
    subplot(nChannels, 1, iPlot);
  
    % Plot the signal against the time vector
    plot(timeVec, data(:, iPlot), 'Color', cm(iPlot, :));

    % Plot the signal against the time vector
    xlim([minTime, maxTime]);

    % Create a label for the Y axis
    ylabel('EEG Amp (uV)');

    % Create a title for the first subplot
    if iPlot == 1
        title(figTitle);
    end
    
    % Create a label for the X axis only for the last subplot
    if iPlot == nChannels
        xlabel('Time (s)');
    end
end
