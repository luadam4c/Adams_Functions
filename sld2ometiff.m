function tiffPaths = sld2ometiff (varargin)
%% Converts each SlideBook file to a directory of OME-TIFF files
% Usage: tiffPaths = sld2ometiff (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       [~, fullPaths] = all_files('Extension', 'sld', 'Recursive', true);
%       if numel(fullPaths) > 1; sld2ometiff('FileNumber', 1); end
%
% Outputs:
%       tiffPaths   - paths to created OME-TIFF files
%                   specified as a cell array 
%                       of cell arrays of character vectors
%
% Arguments:
%       varargin    - 'OutputFormat': output format
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'separated' - separated tiffs for each image plane
%                       'stacked'   - stacked multi-plane tiffs
%                   default == 'separated'
%                   - 'FileNumber': File number to convert
%                   must be a numeric scalar
%                   default == []
%                   - Any other parameter-value pair for TODO()
%
% Requires:
% TODO
%       ~/bfmatlab/
%       cd/all_files.m
%       cd/check_dir.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fileparts.m
%       cd/locate_functionsdir.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-09-09 Created by Adam Lu (adapted from https://docs.openmicroscopy.org/bio-formats/5.7.3/developers/matlab-dev.html)
% 2020-03-31 Now defaults to the creation of multi-page tiffs ('stacked')

%% Hard-coded parameters
validOutputFormats = {'separated', 'stacked'};

%% Default values for optional arguments
outputFormatDefault = 'stacked';
fileNumberDefault = [];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
% iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputFormat', outputFormatDefault, ...
    @(x) any(validatestring(x, validOutputFormats)));
addParameter(iP, 'fileNumber', fileNumberDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));

% Read from the Input Parser
parse(iP, varargin{:});
outputFormat = validatestring(iP.Results.OutputFormat, validOutputFormats);
fileNumber = iP.Results.fileNumber;

% Keep unmatched arguments for the TODO() function
% otherArguments = iP.Unmatched;

%% If not compiled, add directories to search path for required functions
%   Note: If addpath is used, adding `if ~isdeployed` is important to 
%          avoid errors if the function is used as part of a compiled program
%   Note: addpath takes a long time, so use addpath_custom for checking
if ~isdeployed
    % Locate the functions directory
    functionsDirectory = locate_functionsdir;

    % Add path for bfopen.m , etc.
    addpath_custom(fullfile(functionsDirectory, 'bfmatlab'));
end

%% Preparation
% Get all .sld files in the directory
[~, fullPaths] = all_files('Extension', 'sld', 'Recursive', false, ...
                            'ForceCellOutput', true);

% Restrict the file paths if requested
if ~isempty(fileNumber)
    fullPaths = fullPaths(fileNumber);
end

%% Do the job
% Convert each .sld file
tiffPaths = cellfun(@(a) sld2ometiff_helper(a, outputFormat), ...
                    fullPaths, 'UniformOutput', false);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function tiffPaths = sld2ometiff_helper(slidePath, outputFormat)
%% Converts a .sld file to a directory of OME-TIFFs

%% Hard-coded parameters
outExtension = '.tiff';
frameStr = 'frame';
dimensionPlanes = 4;
% TODO: outExtension = '.ome.tiff';

% Print message
fprintf('Converting %s to TIFFs ...\n', slidePath);

% Decide on an output directory for the OME-TIFFs
outDirPath = extract_fileparts(slidePath, 'pathbase');

% Extract the file base for the slidebook file
slideBase = extract_fileparts(slidePath, 'base');

% Make sure that directory exists
check_dir(outDirPath);

% Create a Bio-Formats reader for the SlideBook file
try
    reader1 = bfGetReader(slidePath);
catch
    error('%s could not be read!', slidePath);
end

% Iterate the series until it reaches a series that has more than one frame
iSeries = 1;
reader1.setSeries(1);
while reader1.getImageCount() <= 1
    iSeries = iSeries + 1;
    reader1.setSeries(iSeries);
end

% Get the frame count for this series
nFrames = reader1.getImageCount();

% Read in the OME metadata stored
omeMeta = reader1.getMetadataStore();

% Create frame numbers
frameNumbers = transpose(1:100);
% frameNumbers = transpose(1:nFrames);

switch outputFormat
    case 'separated'
        % Create file bases for the OME-TIFF files
        tiffFileBases = create_labels_from_numbers(frameNumbers, ...
                            'Prefix', [slideBase, '_', frameStr], ...
                            'Suffix', outExtension);

        % Create full paths for the OME-TIFF files
        tiffPaths = fullfile(outDirPath, tiffFileBases);

        % Read and save each frame
        arrayfun(@(x) read_then_save(reader1, x, tiffPaths{x}, omeMeta), ...
                frameNumbers);
    case 'stacked'
        % Create a full path for the OME-TIFF file
        tiffPaths = fullfile(outDirPath, [slideBase, outExtension]);

        % Read each frame
        allFramesCell = arrayfun(@(i) bfGetPlane(reader1, i), frameNumbers, ...
                                'UniformOutput', false);

        % Concatenate along the dimension for multiple image planes
        allFrames = cat(dimensionPlanes, allFramesCell{:});

        % Save the frames as a stacked TIFF file
        % https://www.mathworks.com/help/matlab/ref/tiff.html
        % https://www.awaresystems.be/imaging/tiff/tifftags/baseline.html
        tagStruct.Photometric = Tiff.Photometric.MinIsBlack;
        tagStruct.PlanarConfiguration = Tiff.PlanarConfiguration.Chunky;
        tagStruct.Compression = Tiff.Compression.None;
        tagStruct.SampleFormat = Tiff.SampleFormat.UInt;
        tagStruct.BitsPerSample = 16;
        tagStruct.SamplesPerPixel = 1;
        tagStruct.RowsPerStrip = 1024;
        tagStruct.ImageLength = size(allFrames, 1);
        tagStruct.ImageWidth = size(allFrames, 2);
        save_as_tiff(allFrames, tiffPaths, tagStruct);

        % Save the frames as a stacked OME-TIFF file
        % TODO: bfsave(allFrames, tiffPaths, 'metadata', omeMeta);
    otherwise
        error('outputFormat unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function frame = read_then_save(reader, iFrame, tiffPath, omeMeta)
% Read a frame and save it

% Read the frame
frame = bfGetPlane(reader, iFrame);

% Save the frame as a TIFF file
imwrite(frame, tiffPath);

% Save the frame as an OME-TIFF file
% TODO: bfsave(frame, tiffPath, 'metadata', omeMeta);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function save_as_tiff (imgData, tiffPath, tagStruct)
%% Save image data as a TIFF file
% TODO: Pull out as its own function
%
% File History:
% Adapted from https://www.mathworks.com/matlabcentral/fileexchange/35684-multipage-tiff-stack

% Count the number of planes
nPlanes = size(imgData, 4);

% Create a new Tiff object
tiffObj = Tiff(tiffPath, 'w');

% Write all planes
for iPlane = 1:nPlanes
    % Set the tag info for this image plane
    tiffObj.setTag(tagStruct);

    % Write the data for this image plane
    tiffObj.write(imgData(:, :, :, iPlane));

    % Create a new image directory for the next plane
    if iPlane ~= nPlanes
       tiffObj.writeDirectory();
    end
end

% Close the Tiff object
tiffObj.close();

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Load all the data in the SlideBook file
data = bfopen(slidePath);

% Construct a Bio-Formats reader decorated with the Memoizer wrapper
reader1 = loci.formats.Memoizer(bfGetReader(), 0);

% Initialize the reader2 with an input file to cache the reader2
reader1.setId(slidePath);

% Close reader1
reader1.close()

% Enter parallel loop
parfor iFrame = 1:nFrames
    % Create a path to this OME-TIFF
    tiffPathThis = tiffPaths{iFrame};

    % Initialize logging at INFO level
    bfInitLogging('INFO');

    % Initialize a new reader2 per worker as Bio-Formats is not thread safe
    % Note: Initialization should use the memo file cached before entering the
    %           parallel loop
    reader2 = javaObject('loci.formats.Memoizer', bfGetReader(), 0);

    % Set to the desired input file
    reader2.setId(slidePath);

    % Set to the desired series
    reader2.setSeries(iSeries);

    % Read in the OME metadata stored
    omeMeta = reader2.getMetadataStore();

    % Close reader2
    reader2.close()

    % Get the frame of interest
    frame = bfGetPlane(reader2, iFrame);

    % Save as OME-TIFF
    bfsave(frame, tiffPathThis, 'metadata', omeMeta);
end

% Grab all frames
framesAll = arrayfun(@(x) bfGetPlane(reader1, x), 1:nFrames, ...
                        'UniformOutput', false);

% Save all frames as individual OME-TIFF files
cellfun(@(x, y) bfsave(x, y, 'metadata', omeMeta), framesAll, tiffPaths);

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
