function [nump, pNames, pLabels, pIsLog, pValues, nperp, pchnames, pchvalues, nCells, actMode, loopMode] = extract_looped_params (infolder)
%% Extracts parameters that were looped in the simulation from loopedparams.mat
% USAGE: [nump, pNames, pLabels, pIsLog, pValues, nperp, pchnames, pchvalues, nCells, actMode, loopMode] = extract_looped_params (infolder)
% Arguments:
%     TODO
%
% Requires:
%       infolder/*loopedparams.mat
%       cd/all_ordered_pairs.m
%
% Used by:
%       cd/combine_looped_params.m
%       cd/m3ha_network_tuning_curves.m
%       cd/m3ha_network_tuning_maps.m
%       cd/m3ha_network_raster_plot.m
%       /media/adamX/RTCl/tuning_curves.m

% File History:
% 2017-04-14 Moved from m3ha_network_raster_plot.m
% 2017-04-17 Now extracts pValues if it exists
% 2017-04-17 Removed latency_cells_to_plot
% 2017-05-03 Moved to Adams_Functions
% 2017-05-03 Renamed function get_loopparams -> extract_looped_params
% 2017-05-03 Changed loopedParamsFile 'loopvariables.mat' -> 'loopedparams.mat'
% 2017-05-04 Added loopedParamsFileOld for compatibility with old data
% 2018-05-08 Changed tabs to spaces and limited width to 80

%% Must be consistent with create_looped_params.m & combine_loopparams.m
loopedParamsFile = 'loopedparams.mat';
loopedParamsFileOld = 'loopvariables.mat';

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% TODO: Add required inputs to an Input Parser

%% Find the loopedparams.mat in infolder
loopFiles = dirr(infolder, loopedParamsFile);
loopFilesOld = dirr(infolder, loopedParamsFileOld);
allLoopFiles = [loopFiles; loopFilesOld];
if length(allLoopFiles) > 1
    error('Too many versions of %s or %s!', ...
            loopedParamsFile, loopedParamsFileOld);
elseif length(allLoopFiles) < 1
    error('%s or %s is required!', loopedParamsFile, loopedParamsFileOld);
else
    % There should be only one loopedparams.mat or loopvariables.mat
    m = matfile(fullfile(infolder, allLoopFiles(1).name), ...  
        'Writable', true);          % make writable for pValues and nperp
end

%% Extract pNames, pLabels, pIsLog, nCells, actMode, loopMode
pNames = m.pNames;
pLabels = m.pLabels;
pIsLog = m.pIsLog;
nCells = m.nCells;
if isempty(whos(m, 'actMode'))
                    % for compatibility with older versions of loopedparams.mat
    actMode = 1;
else
    actMode = m.actMode;
end
if isempty(whos(m, 'loopMode'))
                    % for compatibility with older versions of loopedparams.mat
    loopMode = 'cross';
else
    loopMode = m.loopMode;
end

%% Compute nump
nump = numel(pNames);                   % number of parameters changed

%% Extract pValues if it exists and compute nperp
% Otherwise, generate parameter values used from pmin, pmax, pinc
if ~isempty(whos(m, 'pValues'))         % pValues already exists
    % Extract pValues and compute nperp
    pValues = m.pValues;
    nperp = cellfun(@length, pValues);

    % Check if number of parameter values is consistent
    if ~isempty(whos(m, 'nperp'))
        if ~isequal(nperp, m.nperp)
            error('Number of parameter values is inconsistent!');
        end
    elseif ~isempty(whos(m, 'ntrperp')) % older version of loopvariables.m
        if ~isequal(nperp, m.ntrperp)
            error('Number of parameter values is inconsistent!');
        end
    end
else                % pValues don't exist yet
    % Extract pmin, pmax & pinc
    pmin = m.pmin;
    pmax = m.pmax;
    pinc = m.pinc;

    % Generate parameter values used and compute nperp
    pValues = cell(1, nump);
    nperp = zeros(1, nump);
    for p = 1:nump
        
        if pIsLog(p)
            pValues{p} = exp(log(pmin(p)):log(pinc(p)):log(pmax(p)))';
                                            % parameters used as a column vector
        else
            pValues{p} = (pmin(p):pinc(p):pmax(p))';
                                            % parameters used as a column vector
        end
        nperp(p) = length(pValues{p});
    end

    % Save pValues and nperp
    m.pValues = pValues;
    m.nperp = nperp;
end

%% Store parameter names and values for each file (to keep things in order)
if ~isempty(whos(m, 'pchnames')) && ...
    ~isempty(whos(m, 'pchvalues'))  % pchnames & pchvalues already exist
    % Extract pchnames & pchvalues
    pchnames = m.pchnames;          % changed parameter name for each file
    pchvalues = m.pchvalues;        % changed parameter value for each file
else
    switch loopMode
    case 'cross'
        % Create pchnames & pchvalues
        pchnames = cell(1, nfiles);     % changed parameter name for each file
        pchvalues = zeros(1, nfiles);   % changed parameter value for each file
        ct = 0;                         % counts number of files used
        for p = 1:nump
            % Store parameter names and values for each file 
            %   (to keep things in order)
            indices = (ct + 1):(ct + nperp(p));     % indices for this parameter
            pchnames(indices) = repmat(pNames(p), 1, length(indices));
            pchvalues(indices) = (pValues{p})';

            % Update count of number of files used
            ct = ct + nperp(p);
        end
    case 'grid'
        % Repeat each parameter name nperp times
        pNames_rep = cell(1, nump);
                            % stores each parameter name repeated nperp times
        for p = 1:nump
            pNames_rep{p} = repmat(pNames(p), nperp(p), 1);
        end

        % Construct all possible ordered pairs of values and corresponding names
        pchvalues = all_ordered_pairs(pValues);
                                            % a cell array of all ordered pairs
        pchnames = all_ordered_pairs(pNames_rep);
                                            % a cell array of ordered pairs of 
                                            %   corresponding parameter names
    end    

    % Save pchnames & pchvalues
    m.pchnames = pchnames;
    m.pchvalues = pchvalues;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
