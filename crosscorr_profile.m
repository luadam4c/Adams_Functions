function corrProf = crosscorr_profile(data, samplingRate, outFolder, fileBase)
%% data: each channel is a column

%% Hard-coded parameter
windowSizeSeconds = 1;              % window size in seconds
iChannelToPlot = 2;                 % the channel number to plot

%% Create an output folder for cross correlation figures
figFolder = fullfile(outFolder, 'crosscorr');
if exist(figFolder, 'dir') ~= 7
    mkdir(figFolder);
end

%% Preparation
% Count the number of samples
nSamples = size(data, 1);

% Count the number of channels
nChannels = size(data, 2);

% Determine the window size (samples) for the correlation profile
windowSize = windowSizeSeconds * samplingRate;

% Determine the number of windows for the correlation profile
nWindows = ceil(nSamples / windowSize);

% Determine the number of pairs
nPairs = nChannels * (nChannels - 1) / 2;

% Create a time vector in seconds for the windows
timeWindows = (0.5:(nWindows - 0.5))' * windowSizeSeconds;

% Create a time vector in seconds for the entire channel
timeVector = (1:nSamples)' * 1/samplingRate;

% Get the maximum time
maxTime = timeVector(end);

%% Raw data cross correlation with corr2
% Iterate over each channel
indPairs = zeros(nPairs, 2);
ct = 0;
for i = 1:nChannels
    % Only iterate up to current channel number because lower triangle of
    % matrix is same as upper
    for j = 1:i
        if i ~= j
            % Increment the count
            ct = ct + 1;
            
            % Store the index
            indPairs(ct, 1) = i;
            indPairs(ct, 2) = j;            
        end
    end
end
firstOfPairs = indPairs(:, 1);
secondOfPairs = indPairs(:, 2);

% Compute all correlation profiles over time
corrProf = cell(nPairs, 1);
parfor iPair = 1:nPairs
    % Get the indices in corrProf for this pair
    i = firstOfPairs(iPair);
    j = secondOfPairs(iPair);
    
    % Initialize a cross correlation vector over time
    corrThis = zeros(nWindows, 1);

    for k = 1:nWindows
        % Get the starting index of the time window
        idxStart = (k - 1) * windowSize + 1;

        % Get the ending index of the time window
        if k * windowSize <= nSamples
            idxEnd = k * windowSize;
        else
            idxEnd = nSamples;
        end

        % Calculate the cross-correlation coefficient for the
        % two channels in this time window
        corrThis(k) = corr2(data(idxStart:idxEnd, i), ...
                            data(idxStart:idxEnd, j));


    end

    % Store the correlation profile in the cell array
    corrProf{iPair} = corrThis;
end

% Decide on the colormap
cm = colormap(jet(nChannels - 1));

% Plot EEG and correlation coefficient
h = figure;
% Channnel 1 EEG plot
ax1 = subplot(2, 1, 1);
hold on
plot(timeVector, data(:, iChannelToPlot));
ylabel('Voltage(mV)')
title(['EEG for channel ', num2str(iChannelToPlot)]);

% Plot the correlation profiles of each pair of consecutive channels
ax2 = subplot(2, 1, 2);
hold on
for i = 1:(nChannels - 1)
    % Find the corresponding index in indPairs for the pair [i, i+1]
    idxPair = find(firstOfPairs == i + 1 & secondOfPairs == i);
    
    % Correlation profile plot
    corrLabel = ['ch', num2str(i), '-ch', num2str(i + 1)];
    plot(timeWindows, corrProf{idxPair}, ...
         'Color', cm(i, :), 'DisplayName', corrLabel);
    xlabel('Time (seconds)')
    ylabel('Correlation Coefficients')
end
legend('location', 'best');
title('Cross correlation profile');

% Link the x axes
linkaxes([ax1, ax2], 'x');

% Set the x axis limits
xlim([0, maxTime]);

% Create an overarching title
suptitle(fileBase);

% Save the figure
figname = fullfile(figFolder, [fileBase, '_corrprofile']);
saveas(h, figname, 'jpg');

%% Raw data cross correlation with xcorr
for i = 1:(nChannels-1)
    % Only iterate up to current channel number because lower triangle of
    % matrix is same as upper
    [acor, lag] = xcorr(data(:, i), data(:, i+1));

    % Find the lag difference with the largest correlation
    [~, I] = max(abs(acor));
    lagDiff = lag(I);
    
    % Compute the lag in time units
    timeDiff = lagDiff/samplingRate;

    % Create a figure
    h = figure;

    % Plot the the cross correlation
    corrLabel = ['ch', num2str(i), '-ch', num2str(i+1)];
    plot(lag, acor, 'Color', cm(i, :), 'DisplayName', corrLabel);

    % Mark the lagDiff
    a3 = gca;
    a3.XTick = unique([a3.XTick, lagDiff]);

    % Create a text for the lagDiff
    diffSampl = ['lagDiff = ', num2str(lagDiff), ' samples'];
    text(0.55, 0.55, diffSampl, 'Units', 'normalized');

    % Create a legend
    legend('location', 'northeast');
    
    % Create an x label
    xlabel('Lag (samples)')
    
    % Create a y label
    ylabel('Correlation')
    
    % Create a title
    title(['Cross correlation between channel ', num2str(i), ...
            ' and channel ', num2str(i+1)]);
        
    % Save the figure
    figname = fullfile(figFolder, [fileBase, '_', corrLabel]);
    saveas(h, figname, 'jpg');
end

%{
OLD CODE:

channelCompar = {};
    channelCompar = [channelCompar, strcat('ch', num2str(i), ...
                                            '-ch', num2str(i+1))];

legend(channelCompar);

   a3.XTick = sort([-30000000:500000:30000000, lagDiff]);

channelNames = {};
            channelNames = [channelNames, ['ch', num2str(i), ...
                                                '-ch', num2str(j)]];
legend(channelNames)

for i = 1:nChannels
    for j = 1:i
        if i ~= j
            % Correlation profile plot
            corrLabel = ['ch', num2str(i), '-ch', num2str(j)];
            plot(corrProf{i,j}, 'DisplayName', corrLabel);
            xlabel('Time(second)')
            ylabel('Correlation Coefficients')
        end
    end
end

% Creating an empty cell with dimensions nChannels x nChannels
corrProf = cell(nChannels, nChannels);
    % Store the correlation profile in the cell array
    corrProf{i, j} = corrThis;
    plot(corrProf{i, i+1}, 'DisplayName', corrLabel);

%% Bandpass filter the signal (300-4000Hz) and cross correlation with cross2

%% Bandpass filter the signal (0.5-300Hz) and cross correlation with cross2

%% Bandpass filter the signal (300-4000Hz) and cross correlation with xcorr


%% Bandpass filter the signal (0.5-300Hz) and cross correlation with xcorr

hold off

%}