function [startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
    find_initial_slopes(tvecCprAll, ivecCprAll, vvecCprAll, allNSamples, varargin)
%% Find all initial slopes from a set of current pulse responses
% Usage: [startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
%   find_initial_slopes(tvecCprAll, ivecCprAll, vvecCprAll, allNSamples, varargin);
% Arguments:
%       TODO
%       varargin    - 'UseCurrentFlag': whether to use the current trace
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'NSamplesForPlot': the number of samples to average when 
%                                   plotting CPR
%                   must be a positive integer scalar
%                   default == 2
%
% Requires:
%       cd/compute_initial_slopes.m
%
% Used by:
%       cd/m3ha_initial_slopes.m
%

% File History:
% 2018-09-11 Moved from cd/m3ha_initial_slopes.m
% 2018-09-12 Added useCurrentFlag

%% Default values for optional arguments
useCurrentFlagDefault = true;      % use the current trace by default
nSamplesForPlotDefault = 2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'NSamplesForPlot', nSamplesForPlotDefault, ...       
    @(x) validateattributes(x, {'numeric'}, {'scalar', 'positive', 'integer'}));
addParameter(iP, 'UseCurrentFlag', useCurrentFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
nSamplesForPlot = iP.Results.NSamplesForPlot;
useCurrentFlag = iP.Results.UseCurrentFlag;

% Count the number of sweeps
nSwps = numel(tvecCprAll);

% Count the number of possible nSamples to use
nNSamples = length(allNSamples);

% Initialize vector to store slopes
startSlopes = zeros(nSwps, nNSamples);
endSlopes = zeros(nSwps, nNSamples);
avgSlopes = zeros(nSwps, nNSamples);
isUnbalancedAll = zeros(nSwps, nNSamples);

% Initialize vector to store the indices used
indUsedForPlot = zeros(nSwps, 4);

% Find slopes
parfor iSwp = 1:nSwps
% for iSwp = 1:nSwps
%     if iSwp ~= 285 && iSwp ~= 585 && iSwp ~= 297
%         continue;
%     end
    % Read vectors
    tvecCpr = tvecCprAll{iSwp};
    vvecCpr = vvecCprAll{iSwp};
    ivecCpr = ivecCprAll{iSwp};

    % Compute average slopes for all possible nSamples
    startSlopesThis = zeros(1, nNSamples);
    endSlopesThis = zeros(1, nNSamples);
    avgSlopesThis = zeros(1, nNSamples);
    isUnbalancedThis = zeros(1, nNSamples);
    for iNSample = 1:nNSamples
        % Get the current nSamples
        nSamples = allNSamples(iNSample);

        % Compute slopes
        if useCurrentFlag
            [avgSlope, startSlope, endSlope, ...
                indsUsedForPlotSwp, isUnbalanced] = ...
                compute_initial_slopes(tvecCpr, vvecCpr, ...
                                                'IvecCpr', ivecCpr, ...
                                                'NSamples', nSamples);
        else
            [avgSlope, startSlope, endSlope, ...
                indsUsedForPlotSwp, isUnbalanced] = ...
                compute_initial_slopes(tvecCpr, vvecCpr, ...
                                                'NSamples', nSamples);
        end

        % Store slopes
        startSlopesThis(iNSample) = startSlope;
        endSlopesThis(iNSample) = endSlope;
        avgSlopesThis(iNSample) = avgSlope;
        isUnbalancedThis(iNSample) = isUnbalanced;

        % Store indices for plotting
        if nSamples == nSamplesForPlot
            indUsedForPlot(iSwp, :) = indsUsedForPlotSwp;
        end
    end

    % Store the slopes for this sweep
    startSlopes(iSwp, :) = startSlopesThis;
    endSlopes(iSwp, :) = endSlopesThis;
    avgSlopes(iSwp, :) = avgSlopesThis;
    isUnbalancedAll(iSwp, :) = isUnbalancedThis;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%