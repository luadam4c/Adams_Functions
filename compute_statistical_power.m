function output = compute_statistical_power (varargin)
%% Computes statistical power from either raw data or estimated parameters
% Usage: output = compute_statistical_power (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       data0 = 2 * randn(100, 1) + 1;
%       data1 = 2 * randn(50, 1) + 2;
%       compute_statistical_power('DataNull', data0, 'DataAlt', data1)
%       compute_statistical_power('PNull', [1, 2], 'PAlt', 2, 'NSamples', 5:10:100)
%       compute_statistical_power('PNull', [1, 2], 'PAlt', 2, 'NSamples', 5:10:100, 'FullOutput', true)
%
% Outputs:
%       output   - statistical power or everything in a table
%                   specified as a numeric column vector or a table
% Arguments:
%       varargin    - 'TestType': type of test to use
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'z'     - z test for data with known std
%                       't'     - paired or one-sample t test
%                       't2'    - two-sample t test
%                       'var'   - Chi-squared test of variance
%                       'p'     - test for proportion (binomial distribution)
%                   default == 't2'
%                   - 'DataNull': null hypothesis data
%                   must be a numeric array or a cell array of numeric vectors
%                   default == []
%                   - 'DataAlt': alternative hypothesis data
%                   must be a numeric array or a cell array of numeric vectors
%                   default == []
%                   - 'MeanNull': null hypothesis mean
%                   must be a numeric vector
%                   default == estimated from data
%                   - 'MeanAlt': alternative hypothesis mean
%                   must be a numeric vector
%                   default == estimated from data
%                   - 'Stdev': common standard deviation
%                   must be a numeric vector
%                   default == estimated from data
%                   - 'PNull': null hypothesis parameters
%                               Note: For 'z', 't' and 't2' tests
%                                      these are [mean, std] pairs 
%                   must be a numeric array or a cell array of numeric vectors
%                   default == estimated from data
%                   - 'PAlt': alternative hypothesis parameters
%                   must be a numeric vector
%                   default == estimated from data
%                   - 'NSamples': number of samples for both hypotheses
%                   must be a numeric vector
%                   default == computed from data
%                   - 'NSamplesNull': null hypothesis number of samples
%                   must be a numeric vector
%                   default == computed from data
%                   - 'NSamplesAlt': alternative hypothesis number of samples
%                   must be a numeric vector
%                   default == computed from data
%                   - 'FullOutput': whether to output everything in a table
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - Any other parameter-value pair for the sampsizepwr() function
%
% Requires:
%       cd/argfun.m
%       cd/cell2num.m
%       cd/compute_stats.m
%       cd/count_samples.m
%       cd/count_vectors.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%       cd/match_row_count.m
%       cd/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-20 Created by Adam Lu
% 

%% Hard-coded parameters
validTestTypes = {'z', 't', 't2', 'var', 'p'};

%% Default values for optional arguments
testTypeDefault = 't2';
dataNullDefault = [];
dataAltDefault = [];
meanNullDefault = [];
meanAltDefault = [];
stdevDefault = [];
pNullDefault = [];
pAltDefault = [];
nSamplesDefault = [];
nSamplesRatioDefault = [];
nSamplesNullDefault = [];
nSamplesAltDefault = [];
fullOutputDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'TestType', testTypeDefault, ...
    @(x) any(validatestring(x, validTestTypes)));
addParameter(iP, 'DataNull', dataNullDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['DataNull must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'DataAlt', dataAltDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['DataAlt must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'MeanNull', meanNullDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'MeanAlt', meanAltDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'Stdev', stdevDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'PNull', pNullDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['PNull must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'PAlt', pAltDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamples', nSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesRatio', nSamplesRatioDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesNull', nSamplesNullDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesAlt', nSamplesAltDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'FullOutput', fullOutputDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
testType = validatestring(iP.Results.TestType, validTestTypes);
dataNull = iP.Results.DataNull;
dataAlt = iP.Results.DataAlt;
meanNull = iP.Results.MeanNull;
meanAlt = iP.Results.MeanAlt;
stdev = iP.Results.Stdev;
pNull = iP.Results.PNull;
pAlt = iP.Results.PAlt;
nSamples = iP.Results.NSamples;
nSamplesRatio = iP.Results.NSamplesRatio;
nSamplesNull = iP.Results.NSamplesNull;
nSamplesAlt = iP.Results.NSamplesAlt;
fullOutput = iP.Results.FullOutput;

% Keep unmatched arguments for the sampsizepwr() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Match the formats of the two vector sets
[dataNull, dataAlt] = ...
    match_format_vector_sets(dataNull, dataAlt, 'ForceCellOutputs', true);

% Count the number of different values to test
[nVecs(1), nVecs(2), nVecs(3)] = ...
    argfun(@count_vectors, dataNull, dataAlt, pNull);
[nVecs(4), nVecs(5), nVecs(6), nVecs(7), ...
    nVecs(8), nVecs(9), nVecs(10), nVecs(11)] = ...
    argfun(@count_samples, meanNull, meanAlt, stdev, pAlt, ...
          nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio);

% Count the maximum number of values to test
nVectors = max(nVecs);

% Force as column cell arrays
[testType, pNull] = argfun(@force_column_cell, testType, pNull);

% Force as column cell arrays
[meanNull, meanAlt, stdev, pAlt, ...
    nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio] = ...
    argfun(@(x) force_column_cell(x, 'TreatVectorAsElement', false), ...
            meanNull, meanAlt, stdev, pAlt, ...
            nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio);

% Match the row counts of all inputs
[testType, dataNull, dataAlt, meanNull, meanAlt, stdev, pNull, pAlt, ...
    nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio] = ...
    argfun(@(x) match_row_count(x, nVectors), ...
            testType, dataNull, dataAlt, meanNull, meanAlt, stdev, pNull, pAlt, ...
            nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio);

%% Do the job
statPower = ...
    cellfun(@(a, b, c, d, e, f, g, h, i, j, k, l) ...
            compute_statistical_power_helper(a, b, c, d, e, f, ...
                                        g, h, i, j, k, l, otherArguments), ...
            testType, dataNull, dataAlt, meanNull, meanAlt, ...
            stdev, pNull, pAlt, ...
            nSamplesNull, nSamplesAlt, nSamples, nSamplesRatio);

%% Output results in a table
if fullOutput
    % Convert back to numeric vectors
    [pAlt, nSamplesNull, nSamplesAlt] = ...
        argfun(@cell2num, pAlt, nSamplesNull, nSamplesAlt);

    % Place together in a table
    output = table(dataNull, dataAlt, pNull, pAlt, ...
                    nSamplesNull, nSamplesAlt, statPower);
else
    output = statPower;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function statPower = compute_statistical_power_helper (testType, ...
                                dataNull, dataAlt, ...
                                meanNull, meanAlt, stdev, ...
                                pNull, pAlt, nSamplesNull, nSamplesAlt, ...
                                nSamples, nSamplesRatio, otherArguments)
%% Compute the statistical power for one set of data

% Count the number of samples for the null hypothesis
if isempty(nSamplesNull)
    if ~isempty(dataNull)
        nSamplesNull = count_samples(dataNull);
    elseif ~isempty(nSamples)
        nSamplesNull = nSamples;
    else
        fprintf('Number of samples not provided!\n');
    end
end

% Count the number of samples for the alternative hypothesis
if isempty(nSamplesAlt)
    if ~isempty(dataAlt)
        nSamplesAlt = count_samples(dataAlt);
    elseif ~isempty(nSamples) && ~isempty(nSamplesRatio)
        nSamplesAlt = nSamples * nSamplesRatio;
    elseif ~isempty(nSamples)
        nSamplesAlt = nSamples;
    else
        fprintf('Number of samples not provided!\n');
    end
end

% Compute the sample size
%   Note: use the smaller of the two samples for sampsizepwr()
if isempty(nSamples)
    nSamples = min([nSamplesNull, nSamplesAlt]);
end

% Compute the sample size ratio
if isempty(nSamplesRatio)
    % Compute the larger of the two samples
    nSamplesMax = max([nSamplesNull, nSamplesAlt]);

    % Compute the sample size ratio
    nSamplesRatio = nSamplesMax ./ nSamples;
end

% Estimate a standard deviation if not provided
if isempty(stdev)
    stdNull = compute_stats(dataNull, 'std');
    stdAlt = compute_stats(dataAlt, 'std');
    stdev = mean([stdNull, stdAlt]);
end

% Estimate parameters if not provided
if isempty(pNull)
    pNull = estimate_dist_params(dataNull, testType, meanNull, stdev, 'null');
end
if isempty(pAlt)
    pAlt = estimate_dist_params(dataAlt, testType, meanAlt, stdev, 'alt');
end

% Compute the statistical power of the test
statPower = sampsizepwr(testType, pNull, pAlt, [], nSamples, ...
                        'Ratio', nSamplesRatio, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = estimate_dist_params (data, testType, ...
                                        meanValue, stdValue, paramType)
%% Estimates parameters for the underlying distribution

switch testType
    case {'z', 't', 't2'}
        % Compute the mean if not provided
        if isempty(meanValue)
            meanValue = compute_stats(data, 'mean');
        end

        % Put the normal distribution parameters together
        switch paramType
            case 'null'
                params = [meanValue, stdValue];
            case 'alt'
                params = meanValue;
            otherwise
                error('paramType unrecognized!');
         end 
    case 'var'
        % Compute the sample variance
        params = compute_stats(data, 'var');
    case 'p'
        % Estimate the proportion of ones from a Bernoulli trial
        params = sum(data)/numel(data);
    otherwise
        error('testType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%