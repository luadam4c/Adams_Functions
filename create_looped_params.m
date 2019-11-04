function [pchnames, pchvalues, nSims, nParams, pValues, nSimsEachParam] = ...
                create_looped_params (loopMode, pNames, pLabels, pIsLog, ...
                                    pMin, pMax, pInc, varargin)
%% Construct parameters to change for each trial from loopMode, pNames, pIsLog, pMin, pMax, pInc 
% Usage: [pchnames, pchvalues, nSims, nParams, pValues, nSimsEachParam] = ...
%               create_looped_params (loopMode, pNames, pLabels, pIsLog, ...
%                                   pMin, pMax, pInc, varargin)
% Outputs:    
%       pchnames    - a cell array of parameter names or ordered pairs 
%                       of parameter names for each trial
%       pchvalues   - a numeric array of parameter values or a cell array 
%                       of ordered paris of parameter values for each trial
%       nSims       - total number of simulations
%       nParams     - number of different parameters
%       pValues     - a cell array of arrays of parameter values
%       nSimsEachParam  - number of parameter values for each parameter
%
% Arguments:    
%       loopMode    - how to loop through parameters: 'cross' or 'grid'
%                   must be an unambiguous, case-insensitive match 
%                       to one of the following: 
%                       'cross' - Loop through each parameter 
%                                   while fixing others
%                       'grid'  - Loop through all possible 
%                                   combinations of parameters
%       pNames      - names of parameters to loop through
%                   must be a cell array of strings or character arrays
%       pLabels     - labels of parameters to loop through
%                   must be a cell array of strings or character arrays
%       pIsLog      - whether increments of parameters is in log
%                   must be a vector of logicals or 0s/1s
%       pMin        - minimum values of parameters to loop through
%                   must be a numeric vector
%       pMax        - maximum values of parameters to loop through
%                   must be a numeric vector
%       pInc        - increments of parameters to loop through
%                   must be a numeric vector
%       varargin    - 'OutFolder': directory to place outputs, 
%                                   e.g. '/media/shareX/share/'
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'FileLabel': label for output files
%                   must be a string scalar or a character vector
%                   default == datestr(clock, 30);
%                   - 'NCells': number of cells
%                   must be a positive integer
%                   default == 0;
%                   - 'ActMode': activation mode
%                   must be a nonegative integer
%                   default == 0;
%
% Requires:
%       cd/all_ordered_pairs.m
%       cd/argfun.m
%       cd/force_column_vector.m
%
% Used by:
%       cd/m3ha_network_launch.m
%       /media/adamX/RTCl/neuronlaunch.m

% File History:
% 2017-05-03 Moved from /media/adamX/RTCl/neuronlaunch.m
% 2018-05-08 Changed tabs to spaces and limited width to 80


%% Hard-coded parameters
possibleLoopmodes = {'cross', 'grid'};
loopedparamsfile = 'loopedparams.mat';
            % must be consistent with get_loopparams.m & combine_loopparams.m

%% Default values for optional arguments
defaultFileLabel = datestr(clock, 30);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 7
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to an input Parser
addRequired(iP, 'loopMode', ... % how to loop through parameters
    @(x) any(validatestring(x, possibleLoopmodes)));
addRequired(iP, 'pNames', ...   % names of parameters to loop through
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'pLabels', ...  % labels of parameters to loop through
    @(x) assert(iscell(x) && ...
                (min(cellfun(@ischar, x)) || min(cellfun(@isstring, x))), ...
        'Second input must be a cell array of strings or character arrays!'));
addRequired(iP, 'pIsLog', ...   % whether increments of parameters is in log
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary', 'vector'}));
addRequired(iP, 'pMin', ...     % minimum values of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'pMax', ...     % maximum values of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));
addRequired(iP, 'pInc', ...     % increments of parameters to loop through
    @(x) validateattributes(x, {'numeric'}, {'vector'}));

% Add parameter-value pairs to the input Parser
addParameter(iP, 'OutFolder', '', ...   % directory to place outputs
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));    
                                                    % introduced after R2016B
addParameter(iP, 'FileLabel', defaultFileLabel, ...     % label for output files
    @(x) validateattributes(x, {'char', 'string'}, {'nonempty'}));
%    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                    % introduced after R2016B
addParameter(iP, 'NCells', 0, ...       % number of cells
    @(x) validateattributes(x, {'numeric'}, {'positive', 'integer'}));
addParameter(iP, 'ActMode', 0, ...      % activation mode
    @(x) validateattributes(x, {'numeric'}, {'nonnegative', 'integer'}));

% Read from the input Parser
parse(iP, loopMode, pNames, pLabels, pIsLog, pMin, pMax, pInc, varargin{:});
loopMode = validatestring(loopMode, possibleLoopmodes);
outfolder = iP.Results.OutFolder;
filelabel = iP.Results.FileLabel;
nCells = iP.Results.NCells;
actMode = iP.Results.ActMode;

% Set default arguments
if isempty(outfolder)
    outfolder = pwd;
end

% Check relationships between arguments
if numel(unique([numel(pNames), numel(pLabels), numel(pIsLog), ...
                numel(pMin), numel(pMax), numel(pInc)])) ~= 1
    error('pNames, pLabels, pIsLog, pMin, pMax, pInc not all equal length!!');
end

%% Calculate number of parameters to loop through
nParams = numel(pNames);
pValues = cell(1, nParams);     % parameter values that will be used
nSimsEachParam = zeros(1, nParams);      % number of simulations per parameter type
for p = 1:nParams
    % Generate parameter values that will be used
    if pIsLog(p)
        pValues{p} = exp(log(pMin(p)):log(pInc(p)):log(pMax(p)))';
    else
        pValues{p} = (pMin(p):pInc(p):pMax(p))';
    end
    nSimsEachParam(p) = length(pValues{p});
end

%% Find total number of simulations
switch loopMode
case 'cross'
    % Total number of simulations in 'cross' mode
    nSims = sum(nSimsEachParam);
case 'grid'
    % Total number of simulations in 'grid' mode
    nSims = prod(nSimsEachParam);
end

%% Create pchnames & pchvalues
switch loopMode
case 'cross'
    pchnames = cell(nSims, 1);    % stores parameter name for each trial
    pchvalues = zeros(nSims, 1);  % stores parameter value for each trial
    ct = 0;                         % counts number of simulations used
    for p = 1:nParams
        % Store parameter names and values for each trial 
        %   (to keep things in order)
        indices = (ct + 1):(ct + nSimsEachParam(p));     % indices for this parameter
        pchnames(indices) = repmat(pNames(p), 1, length(indices));
                                                % parameter name for each trial
        pchvalues(indices) = (pValues{p})';     % parameter value for each trial 

        % Update count of number of files used
        ct = ct + nSimsEachParam(p);
    end
case 'grid'
    % Repeat each parameter name nSimsEachParam times
    pNamesRepeated = cell(1, nParams);     % stores each parameter name 
                                    %   repeated nSimsEachParam times
    for p = 1:nParams
        pNamesRepeated{p} = repmat(pNames(p), nSimsEachParam(p), 1);
    end

    % Construct all possible ordered pairs of values and corresponding names
    pchvalues = all_ordered_pairs(pValues); 
                                    % a cell array of all ordered pairs
    pchnames = all_ordered_pairs(pNamesRepeated);
                                    % a cell array of ordered pairs of 
                                    %   corresponding parameter names
end

% Force as columns
[pchnames, pchvalues] = ...
    argfun(@(x) force_column_vector(x, 'TreatCellAsArray', true), ...
            pchnames, pchvalues);

%% Save looping variables in a mat file named by the date & time
save(fullfile(outfolder, sprintf('%s_%s', filelabel, loopedparamsfile)), ...
    'loopMode', 'nSims', 'nParams', 'nSimsEachParam', ...
    'pNames', 'pLabels', 'pIsLog', 'pMin', 'pMax', 'pInc', ...
    'pValues', 'pchnames', 'pchvalues', ...
    'nCells', 'actMode', '-v7.3');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
