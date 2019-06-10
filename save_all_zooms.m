function save_all_zooms (fig, figPath, zoomWins, varargin)
%% Save a figure in various zoom windows
% Usage: save_all_zooms (fig, figPath, zoomWins, varargin)
% Arguments:
%       fig         - figure to save
%                   must be a a figure object handle or a Simulink block diagram
%       figPath     - figure path or name
%                   must be a string scalar or a character vector
%       zoomWins    - windows to zoom to
%                   must be empty or a numeric vector with 2 elements,
%                       or a numeric array with 2 rows
%                       or a cell array of numeric arrays
%                   default == []
%       varargin    - 'OutFolder': the name of the directory that 
%                                       plots will be placed
%                   must be a string scalar or a character vector
%                   default == detected from figPath or pwd
%                   - 'Figtypes': figure type(s) for saving; 
%                       e.g., 'png', 'fig', or {'png', 'fig'}, etc.
%                   could be anything recognised by the built-in 
%                       saveas() function 
%                   (see isfigtype.m under Adams_Functions)
%                   default == 'png'
%                   
% Requires:
%       cd/check_subdir.m
%       cd/count_vectors.m
%       cd/create_labels_from_numbers.m
%       cd/extract_fileparts.m
%       cd/force_column_cell.m
%       cd/force_string_end.m
%       cd/isfigtype.m
%       cd/save_all_figtypes.m
%
% Used by:
%       cd/parse_multiunit.m

% File History:
% 2019-06-03 Pulled from parse_multiunit.m
% 

%% Hard-coded parameters
fullLabel = 'full';

%% Default values for optional arguments
outFolderDefault = '';              % set later
figTypesDefault = 'png';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 3
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfigPath);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an Input Parser
addRequired(iP, 'fig');
    % TODO: validation function %);
addRequired(iP, 'figPath', ...
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
addRequired(iP, 'zoomWins', ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['zoomWins must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));

% Add optional inputs to the Input Parser
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Figtypes', figTypesDefault, ...
    @(x) min(isfigtype(x, 'ValidateMode', true)));

% Read from the Input Parser
parse(iP, fig, figPath, zoomWins, varargin{:});
outFolder = iP.Results.OutFolder;
[~, figtypes] = isfigtype(iP.Results.Figtypes, 'ValidateMode', true);

%% Preparation

% Make sure figPath ends with the extension
figExt = extract_fileparts(figPath, 'extension');
figPath = force_string_end(figPath, ['.', figExt]);

% Decide on the output directory
if isempty(outFolder)
    outFolder = extract_fileparts(figPath, 'directory');

    if isempty(outFolder)
        outFolder = pwd;
    end
end

% Extract the file base
figBase = extract_fileparts(figPath, 'base');

% Force as a column cell array of column vectors
zoomWins = force_column_cell(zoomWins);

% Count the number of zoom windows
nWins = count_vectors(zoomWins);

% Create labels for zoom windows
zoomLabels = create_labels_from_numbers(1:nWins, 'Prefix', 'zoom');

% Create subdirectories for different zooms
check_subdir(outFolder, {fullLabel, zoomLabels{:}});

%% Do the job
% Get the figure
figure(fig);

% Save the full figure
figPathFull = fullfile(outFolder, fullLabel, [figBase, '_', fullLabel]);
save_all_figtypes(fig, figPathFull, 'png');

% Save zoomed figures
for iWin = 1:nWins
    % Get the current zoom label
    zoomLabelThis = zoomLabels{iWin};

    % Change x limits
    xlim(zoomWins{iWin});

    % Create full path to figure
    figPath = fullfile(outFolder, zoomLabelThis, [figBase, '_', zoomLabelThis]);

    % Save figure
    save_all_figtypes(fig, figPath, 'png');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
