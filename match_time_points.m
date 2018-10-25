function dataNew = match_time_points (dataOld, timeNew)
%% Interpolates data (containing a time column) to match the time points of a new time vector
% Usage: dataNew = match_time_points (dataOld, timeNew)
%
% Outputs:
%       dataNew     - new data array
%                   specified as a numeric array
%                       with the first column being the time column
% Arguments:
%       dataOld     - original data array
%                   first column must be the time column
%                   must be a numeric array
%       timeNew     - new time vector
%                   must be a numeric vector
%
% Requires:
%       cd/print_or_show_message.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/run_neuron_once_4compgabab.m

%
% File History:
% 2018-08-10 Adapted from code in run_neuron_once_4compgabab.m
% TODO: Input parser
%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Make sure the new time vector is a column
timeNew = timeNew(:);

% Count the number of samples in the new time vector
nSamples = length(timeNew);

% Get the number of columns in the old data
nCols = size(dataOld, 2);

% Preallocate a new data array
dataNew = zeros(nSamples, nCols);

% Place the new time vector in the first column of the new data array
dataNew(:, 1) = timeNew;

% No need to proceed if there are no other columns
if nCols < 2
    message = 'There are no columns to interpolate other than the time column!';
    mTitle = 'Match Time Points Warning';
    icon = 'warn';
    print_or_show_message(message, 'MTitle', mTitle, 'Icon', icon);
    return;
end

%% Interpolate
% Extract the old time vector
timeOld = dataOld(:, 1);

% Find the unique time points and the corresponding indices in timeOld
[timeUnique, indUnique] = unique(timeOld);

% Interpolate all other columns
dataNew(:, 2:end) = interp1q(timeUnique, dataOld(indUnique, 2:end), timeNew);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Note: This is actually slower:
% Do for each other column of the old data array
parfor iCol = 2:nCols
    % Get the old column vector
    vecOld = dataOld(indUnique, iCol);

    % Interpolate to make a new column vector
    vecNew = interp1q(timeUnique, vecOld, timeNew);

    % Place the new column vector in the new data array
    dataNew(:, iCol) = vecNew;
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%