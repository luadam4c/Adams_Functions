function [startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
    find_initial_slopes(tvecCprAll, ivecCprAll, vvecCprAll, allNSamples);
%% Find all initial slopes from a current pulse response
% Usage: [startSlopes, endSlopes, avgSlopes, isUnbalancedAll, indUsedForPlot] = ...
%   find_initial_slopes(tvecCprAll, ivecCprAll, vvecCprAll, allNSamples);
% Arguments:
%       TODO
%
% Requires:
%       /home/Matlab/Adams_Functions/compute_average_initial_slopes.m
%
% Used by:
%       /media/adamX/m3ha/data_dclamp/initial_slopes.m
%

% File History:
% 2018-09-11 Moved from /media/adamX/m3ha/data_dclamp/initial_slopes.m

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 4
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
% TODO

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
                compute_average_initial_slopes(tvecCpr, vvecCpr, ...
                                                'IvecCpr', ivecCpr, ...
                                                'NSamples', nSamples);
        else
            [avgSlope, startSlope, endSlope, ...
                indsUsedForPlotSwp, isUnbalanced] = ...
                compute_average_initial_slopes(tvecCpr, vvecCpr, ...
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