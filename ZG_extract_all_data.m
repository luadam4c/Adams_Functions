function ZG_extract_all_data (varargin)
%% Extract each cell and concatenate sweeps into a single file under a subdirectory of its own
% Usage: ZG_extract_all_data (varargin)
%
% Requires:
%       cd/find_in_strings.m
%
% Used by:
%       /home/barrettlab/detect_with_minEASE/ZG_extract_all_data.sh
%
% File History:
%   2018-02-09 Created by Adam Lu
%   2018-02-09 sliceNum changed to sliceStr
%   2018-02-13 Now skips a subdirectory if the slice is already processed
%   2018-02-21 Now replaces values more than 10 standard deviations 
%               away with an averaged value
%   2018-02-26 Changed outFolder to all_data
%   2018-08-01 Moved groupLabel extraction from master list here from ZG_extract_all_IEIs.m

%% Hard-coded constants
%   Note: Must be consistent with ZG_remedy_file_names.sh
masterListFile = '/home/barrettlab/holySheet/CaImagingExperiments_MasterList.xlsx';
masterListTab = '4mark';
dataDirNamesColNum = 3;     % dataDirNames column number in the master list tab
groupNamesColNum = 5;       % groupNames column number in the master list tab
groupNameDefault = 'Unnamed';                  % group label used if 
                                                %   not provided in master list
inFolder = '/home/barrettlab/imaging_data/';    % input data folder
outFolder = '/home/barrettlab/detect_with_minEASE/all_data/';
                                                % output folder
logFolder = '/home/barrettlab/detect_with_minEASE/';
                                                % folder for placing log files
subdirSuffix = '_caltracer';                    % suffix for subdirectories
subdirPrefix = 'Data';                          % prefix for subdirectories
datafileSuffix = '_exp.mat';                    % suffix for data files

% Define the pattern for subdirectories
patternSubDir = [subdirPrefix, '*', subdirSuffix];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Create a diary for the standard output
currentDateTime = datestr(clock, 'yyyymmddTHHMM');
logFileName = ['ZG_extract_all_data_', currentDateTime, '.log'];
diary(fullfile(logFolder, logFileName));

% Find all the subdirectories (each is a slice) to extract
allSubDir = dir(fullfile(inFolder, patternSubDir));
nSubDir = length(allSubDir);                    % number of data subdirectories

% Create output directory if doesn't exist
if exist(outFolder, 'dir') ~= 7
    mkdir(outFolder);
    fprintf('New directory made: %s\n\n', outFolder);
end

% Read slice groupings
[~, ~, masterList] = xlsread(masterListFile, masterListTab);
dataDirNamesOrig = masterList(2:end, dataDirNamesColNum);
groupLabelsOrig = masterList(2:end, groupNamesColNum);

% Convert the contents to strings if not already so
dataDirNames = cellfun(@num2str, dataDirNamesOrig, 'UniformOutput', false);
groupLabels = cellfun(@num2str, groupLabelsOrig, 'UniformOutput', false);

% Replace any empty group label with the default
isNotLabelled = cellfun(@(x) strcmp(x, 'NaN'), groupLabels);
nNotLabelled = sum(isNotLabelled);
groupLabels(isNotLabelled) = repmat({groupNameDefault}, nNotLabelled, 1);

% Add an 'x' to the beginning of group labels so that they can be 
%   valid field names of structures
groupLabels = cellfun(@(x) ['x', x], groupLabels, 'UniformOutput', false);

% Go through all data subdirectories
for iSubDir = 1:nSubDir
    % Get name of this data subdirectory
    subDirName = allSubDir(iSubDir).name;

    % Extract slice number 
    splitted = strsplit(subDirName, '_');
    scanned = textscan(splitted{1}, 'Data%s');
    sliceStr = scanned{1}{1};

    % Create the slice label
    sliceLabel = ['data', sliceStr];

    % Find the matching slice label in the dataDirNames cell array
    %   Note: Must search for substrings and ignore case because 
    %           dataDirNames has elements of the form DataXXX_caltracer_ORAMA
    idxSlice = find_in_strings(sliceLabel, dataDirNames, ...
                                    'SearchMode', 'substrings', ...
                                    'IgnoreCase', true);

    if isempty(idxSlice)
        fprintf('Cannot find the slice label %s in dataDirNames!!\n', sliceLabel);
        % Use the default group label (prepended by 'x')
        groupLabel = ['x', groupNameDefault];
    else
        % Get the corresponding group label
        groupLabel = groupLabels{idxSlice};
    end

    % Construct a groupSlice label
    groupSliceLabel = [groupLabel, '_', sliceLabel];

    % Skip this directory if output files already exist
    slicePattern = fullfile(outFolder, [groupSliceLabel, '_*']);
    if ~isempty(dir(slicePattern)) 
        continue;
    end

    % Find the pattern for data files in this subdirectory
    patternFile = ['*(*)', datafileSuffix];

    % Find all the data files in this directory
    dataFiles = dir(fullfile(inFolder, subDirName, patternFile));
                    
    % Get the total number of data files
    %   This should be the number of sweeps
    nFiles = length(dataFiles);
    nSwps = nFiles;

    % Print warning if no files found
    if nFiles == 0
        fprintf(['There are no files matching the pattern %s ', ...
                'in the subdirectory %s!!\n\n'], patternFile, subDirName);
        continue;
    end

    % For each sweep, extract the data traces while checking 
    %   whether the number of cells (number of rows) is consistent
    nCells = [];                % initialize as empty
    allData = cell(nSwps, 1);  % stores all sweeps and all cells
    for iSwp = 1:nSwps
        % Find the pattern for the data file for this sweep
        patternSwp = ['*(', num2str(iSwp), ')', datafileSuffix];

        % Find all the files for this sweep
        swpFiles = dir(fullfile(inFolder, subDirName, patternSwp));
        
        % Check if there is exactly one file
        nSwpFiles = length(swpFiles);
        if nSwpFiles < 1
            fprintf(['There should be at least one file with the pattern ', ...
                    '%s in the subdirectory %s!!\n\n'], patternSwp, subDirName);
            nCells = [];
            break;
        elseif nSwpFiles > 1
            fprintf(['There should be only one file with the pattern %s ', ...
                    'in the subdirectory %s!!\n\n'], patternSwp, subDirName);
            nCells = [];
            break;
        end

        % Get the name of the matfile for this sweep
        matfileThisSwp = swpFiles(1).name;

        % Load the contents of the matfile
        m = load(fullfile(inFolder, subDirName, matfileThisSwp));

        % Check if any data exists for this sweep
        if ~isfield(m.handles.app.experiment, 'traces')
            fprintf(['Trace data is missing for the file %s ', ...
                    'in the subdirectory %s!!\n\n'], ...
                    matfileThisSwp, subDirName);
            nCells = [];
            break;
        end

        % Extract the data for this sweep
        %   and transpose it so that each column is a cell
        dataThisSwp = m.handles.app.experiment.traces';
        
        % Check the number of cells (now columns)
        nCellsThisSwp = size(dataThisSwp, 2);
        if isempty(nCells)
            nCells = nCellsThisSwp;
        else
            if nCellsThisSwp ~= nCells
                fprintf(['The number of cells in the file %s', ...
                            ' is inconsistent!!\n\n'], matfileName);
                nCells = [];
                break;
            end
        end

        % Put the data for this sweep in allData
        allData{iSwp} = dataThisSwp;
    end

    % Check if process was terminated
    if isempty(nCells)
        fprintf(['The files in the subdirectory %s ', ...
                    'are inconsistent!!\n\n'], subDirName);
        continue;
    end

    % Concatenate the data from all sweeps
    dataAllCells = cell2mat(allData);

    % Extract the data for each cell
    for iCell = 1:nCells
        % Create a cell label
        cellLabel = ['cell', num2str(iCell)];

        % Construct groupSliceCell label for this cell
        groupSliceCell = [groupLabel, '_', sliceLabel, '_', cellLabel];

        % Construct subdirectory for this cell
        subdirThisCell = fullfile(outFolder, groupSliceCell);
        if exist(subdirThisCell, 'dir') ~= 7
            mkdir(subdirThisCell);
            fprintf('New directory made: %s\n\n', subdirThisCell);
        end
                            
        % Construct file name for this cell
        matfileNameThisCell = [groupSliceCell, '.mat'];

        % Extract the data for this cell
        dataThisCell = dataAllCells(:, iCell);

        % Get the number of data values
        nValues = length(dataThisCell);
        
        % Get the mean and standard deviation of all data values
        meanThis = mean(dataThisCell);
        stdThis = std(dataThisCell);

        % Compute the successive value differences for thie cell
        diffData = diff(dataThisCell);

        % Compute the standard deviation of all absolute differences
        stdAbsDiff = std(abs(diffData));

        % Check if there are any indices with abnormal values:
        %   Let the two preceding differences be dp2, dp1
        %   and let the two trailing differences be dt1, dt2
        %   then the value is abnormal if 
        %       dp1 > ratio * mean([abs(dp2), abs(dt2)]) &
        %       dt1 < -ratio * mean([abs(dp2), abs(dt2)])
        %   OR
        %       dp1 < -ratio * mean([abs(dp2), abs(dt2)]) &
        %       dt1 > ratio * mean([abs(dp2), abs(dt2)])
        %   AND if abs(dp1) is more than ratio times stdAbsDiff
        ratio = 10;
        baseDiff = mean(abs([diffData(1:(end-3)), diffData(4:end)]), 2);
        baseDiff = [abs(diffData(3)); baseDiff; abs(diffData(end-2))];
        thresDiff = ratio * baseDiff;
        indTooLarge = find(((diffData(1:(end-1)) > thresDiff & ...
                             diffData(2:end) < -thresDiff) | ...
                            (diffData(1:(end-1)) < -thresDiff & ...
                             diffData(2:end) > thresDiff)) & ...
                           abs(diffData(1:(end-1))) > ratio * stdAbsDiff ...
                          ) + 1;

        % If such a value exist, replace it with the average of the 
        %   previous and next values
        if ~isempty(indTooLarge)
            % Find the total number of indices that are too large
            nIndTooLarge = length(indTooLarge);

            % Replace the data values that are too large and print
            %   message
            for iInd = 1:nIndTooLarge
                % Get the index of dataThisCell to change
                idx = indTooLarge(iInd);

                % Store the value at that index
                valueFrom = dataThisCell(idx);

                % Decide on the value to change to
                if idx == 1
                    % If there is no previous value, use the next value
                    valueTo = dataThisCell(idx + 1);
                elseif idx == nValues
                    % If there is no next value, use the previous value
                    valueTo = dataThisCell(idx - 1);
                else
                    % Otherwise, use the mean of the previous and next values
                    valueTo = mean([dataThisCell(idx - 1), ...
                                    dataThisCell(idx + 1)]);
                end

                % Change the value
                dataThisCell(idx) = valueTo;

                % Print message to standard output
                fprintf(['Value at index %d of data for %s was ', ...
                         'changed from %g to %g!!\n\n'], ...
                        idx, groupSliceCell, valueFrom, valueTo);
            end
        end

        % Save the file in its own directory
        save(fullfile(subdirThisCell, matfileNameThisCell), ...
                'dataThisCell');
    end
end

% Close the diary for standard output
diary off;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

    for iFile = 1:nFiles
        % Get the name of the file
        matfileName = dataFiles(iFile).name;

currentDateTime = datestr(datetime('now', 'Format', 'yyyyMMdd''T''HHmm'));

    patternFile = ['*', num2str(sliceNum), '(*)', datafileSuffix];
        patternSwp = ['*', num2str(sliceNum), ...
                        '(', num2str(iSwp), ')', datafileSuffix];

    scanned = textscan(subDirName, 'Data%f');   
    sliceNum = scanned{1};
        sliceCell = ['data', num2str(sliceStr), '_cell', num2str(iCell)];

        thresThis = meanThis + 10 * stdThis;

        % Check if there is any data value more than 100 times 
        %   the standard deviation away from the mean
        thresThis = meanThis + 100 * stdThis;
        indTooLarge = find(dataThisCell > thresThis);
        
        % Check if there is any data value more than 1e4 
        indTooLarge = find(dataThisCell > 1e4);

        % Check if there are successive value differences more than 10 times 
        %   the standard deviation of all value differences 
        thresDiff = 10 * stdDiff;
        indTooLarge = find(diffData(1:(end-1)) > thresDiff & ...
                            diffData(2:end) < -thresDiff) + 1;
        if abs(diffData(1) < -thresDiff)
            indTooLarge = [1, indTooLarge];
        end
        if abs(diffData(end) > thresDiff)
            indTooLarge = [indTooLarge, nValues];
        end

        % Compute the standard deviation of all value differences
        stdDiff = std(diffData);

        % Check if there are any indices with abnormally large values:
        %   Let the two preceding differences be dp2, dp1
        %   and let the two trailing differences be dt1, dt2
        %   then the value is abnormally large if 
        %       dp1 > ratio * mean([abs(dp2), abs(dt2)]) &
        %       dt1 < -ratio * mean([abs(dp2), abs(dt2)])
        %   and if the value is more than 2 standard deviations away
        %       from the mean
        ratio = 500;
        baseDiff = mean(abs([diffData(1:(end-3)), diffData(4:end)]), 2);
        baseDiff = [abs(diffData(3)); baseDiff; abs(diffData(end-2))];
        thresDiff = ratio * baseDiff;
        indTooLarge = find(diffData(1:(end-1)) > thresDiff & ...
                            diffData(2:end) < -thresDiff & ...
                            dataThisCell(2:(end-1)) > meanThis + 2 * stdThis ...
                            ) + 1;

        % Check if there are any indices with abnormal values:
        %   Let the two preceding differences be dp2, dp1
        %   and let the two trailing differences be dt1, dt2
        %   then the value is abnormal if 
        %       dp1 > ratio * mean([abs(dp2), abs(dt2)]) &
        %       dt1 < -ratio * mean([abs(dp2), abs(dt2)])
        %   OR
        %       dp1 < -ratio * mean([abs(dp2), abs(dt2)]) &
        %       dt1 > ratio * mean([abs(dp2), abs(dt2)])
        %   AND if abs(dp1) is more than 10 times stdAbsDiff
        ratio = 10;
        baseDiff = mean(abs([diffData(1:(end-3)), diffData(4:end)]), 2);
        baseDiff = [abs(diffData(3)); baseDiff; abs(diffData(end-2))];
        thresDiff = ratio * baseDiff;
        indTooLarge = find(((diffData(1:(end-1)) > thresDiff & ...
                             diffData(2:end) < -thresDiff) | ...
                            (diffData(1:(end-1)) < -thresDiff & ...
                             diffData(2:end) > thresDiff)) & ...
                           abs(diffData(1:(end-1))) > 10 * stdAbsDiff ...
                          ) + 1;
                          
outFolder = '/home/barrettlab/detect_with_minEASE/all_output/';

% Construct sliceCell string for this cell
sliceCell = ['data', sliceStr, '_cell', num2str(iCell)];

% Construct subdirectory for this cell
subdirThisCell = fullfile(outFolder, sliceCell);
if exist(subdirThisCell, 'dir') ~= 7
    mkdir(subdirThisCell);
    fprintf('New directory made: %s\n\n', subdirThisCell);
end
                    
% Construct file name for this cell
matfileNameThisCell = [sliceCell, '.mat'];
        % Print message to standard output
        fprintf(['Value at index %d of data for %s was ', ...
                 'changed from %g to %g!!\n\n'], ...
                idx, sliceCell, valueFrom, valueTo);

% Skip this directory if output files already exist
slicePattern = fullfile(outFolder, [sliceLabel, '_*']);
if ~isempty(dir(slicePattern)) 
    continue;
end

%}
