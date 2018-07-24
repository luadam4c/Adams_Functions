function [data, file, dataAll] = read_adicht(fileName)

% fileName = 'C:\Users\Pinn Analysis\Desktop\Shinnosuke\data\cobalt multi 07192018.adicht';

%% Add the path to the adi functions
addpath('C:\Users\Pinn Analysis\Desktop\Shinnosuke\JimHokanson-adinstruments_sdk_matlab-393033f')

%% Read the .adicht file
file = adi.readFile(fileName);

% Get the number of records
nRecords = file.n_records;

% Get the number of channels
nChannels = file.n_channels;

%% Read all the data
dataAll = cell(nChannels, nRecords);
for iChannel = 1:nChannels
    % Get the channel name
    channelName = ['Channel ', num2str(iChannel)];
    for iRecord = 1:nRecords       
        % Get the data for this channel and record
        dataAll{iChannel, iRecord} = ...
            file.getChannelByName(channelName).getData(iRecord);
    end
end

%%
% Find the total number of samples for each channel
allLengths = cellfun(@length, dataAll);
allNSamples = sum(allLengths, 2);

% Check if all channels have the same number of samples
if length(unique(allNSamples)) > 1
    error('The channels do not have the same number of data points!');
end

% Find the number of samples for each channel
nSamples = allNSamples(1);

% Concatenate data from each channel
data = zeros(nSamples, nChannels);
for iChannel = 1:nChannels
    % Preallocate data for this channel
    dataThis = zeros(nSamples, 1);

    % Place the records in the vector in turn
    idxStart = 1;           	% starting index for the next record
    for iRecord = 1:nRecords
        % Get the data vector for this record
        dataThisRecord = dataAll{iChannel, iRecord};
        
        % Get the length of the data vector for this record
        lengthThisRecord = allLengths(iChannel, iRecord);
        
        % Get the ending index
        idxEnd = (idxStart - 1) + lengthThisRecord;
        
        % Place the data vector for this record in the data
        %   vector for this channel
        dataThis(idxStart:idxEnd, 1) = dataThisRecord;
        
        % Update the starting index for the next record
        idxStart = idxEnd + 1;
    end

    % Store the data vector for this channel in the cell array
    data(:, iChannel) = dataThis;
end
