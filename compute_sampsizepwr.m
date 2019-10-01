function output = compute_sampsizepwr (varargin)
%% Computes the sample size needed, the statistical power or the alternative hypothesis parameter from either raw data or estimated parameters
% Usage: output = compute_sampsizepwr (varargin)
% Explanation:
%       This function uses sampsizepwr() but accepts data vectors as
%           arguments and allows the user to input multiple parameters at once.
%       The default output is a table that contains everything,
%           but this can be changed with the 'OutputType' optional argument
%       The default test type is a two-sample t-test,
%           but this can be changed with the 'TestType' optional argument, 
%           which can be a cell vector as well
%       The default statistical power desired is 0.80, 
%           but this can be changed with the 'StatPower' optional argument, 
%           which can be a numeric vector as well
%       Basically, all optional parameters can be passed as cell vectors,
%           and the program will automatically compute the maximum number
%           of rows that will accommodate all distinct inputs. 
%       Not yet implemented: (TODO) outer product of two input vectors with
%           different lengths that are greater than 1
%       To save your results, use the 'SaveSheetFlag' or 'SaveMatFlag'
%           options.
%
% Example(s):
%       data0 = 2 * randn(100, 1) + 1;
%       data1 = 2 * randn(50, 1) + 2;
%       compute_sampsizepwr('DataCon', data0, 'DataExp', data1)
%       compute_sampsizepwr('DataCon', data0, 'DataExp', data1, 'TestType', {'z', 't', 't2'})
%       compute_sampsizepwr('DataCon', data0, 'PAlt', 1.5:0.1:2)
%       compute_sampsizepwr('PNull', [1, 2], 'PAlt', 2, 'NSamples', 5:10:100)
%       sampsize = compute_sampsizepwr('PNull', [1, 2], 'PAlt', 2, 'OutputType', 'sampsize')
%       sampsize = compute_sampsizepwr('PNull', [1, 2], 'PAlt', 2, 'StatPower', 0.8:0.01:0.95, 'OutputType', 'sampsize')
%       power = compute_sampsizepwr('PNull', [1, 2], 'PAlt', 2, 'NSamples', 5:10:100, 'OutputType', 'power')
%       pAlt = compute_sampsizepwr('PNull', [1, 2], 'NSamples', 5:10:100, 'OutputType', 'pAlt')
%
% Outputs:
%       output      - output specified by outputType or everything in a table
%                   specified as a numeric column vector or a table
%
% Arguments:
%       varargin    - 'OutputType': type of output
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'all'       - everything in a table
%                       'scalars'   - everything except data vectors
%                       'sampsize'  - sample size needed
%                       'statpower' - computed statistical power
%                       'pAlt'      - parameter for alternative hypothesis
%                   default == 'all'
%                   - 'TestType': type of test to use
%                   must be an unambiguous, case-insensitive match to one of: 
%                       'z'     - paired z test
%                       't'     - paired t test
%                       't2'    - two-sample t test
%                       'var'   - Chi-squared test of variance
%                       'p'     - test for proportion (binomial distribution)
%                   or a cell array of such strings
%                   default == 't2'
%                   - 'DataCon': control group data
%                   must be a numeric array or a cell array of numeric vectors
%                   default == []
%                   - 'DataExp': experimental group data
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
%                   - 'StatPowerDesired': statistical power desired
%                   must be a numeric vector
%                   default == 0.80
%                   - 'NSamplesCon': null hypothesis number of samples
%                   must be a numeric vector
%                   default == computed from data
%                   - 'NSamplesExp': alternative hypothesis number of samples
%                   must be a numeric vector
%                   default == computed from data
%                   - 'NSamples': number of samples for both hypotheses
%                   must be a numeric vector
%                   default == computed from data
%                   - 'NSamplesRatio': sample size ratio
%                   must be a numeric vector
%                   default == computed from data
%                   - 'SaveSheetFlag': whether to save values in a spreadsheet
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveMatFlag': whether to save values in a .mat file
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
%       cd/extract_elements.m
%       cd/force_column_cell.m
%       cd/iscellnumeric.m
%       cd/match_format_vector_sets.m
%       cd/match_row_count.m
%       cd/print_cellstr.m
%       cd/struct2arglist.m
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2019-08-20 Created by Adam Lu
% 2019-09-20 Added 'scalars' as an outputType
% 2019-10-01 Renamed dataNull, dataAlt as dataCon, dataExp
% 

%% Hard-coded parameters
validOutputTypes = {'all', 'scalars', 'sampsize', 'power', 'pAlt'};
validTestTypes = {'z', 't', 't2', 'var', 'p'};

% TODO: Make optional arguments
fileSuffix = 'sample_size_power_analysis';

%% Default values for optional arguments
outputTypeDefault = 'all';
testTypeDefault = 't2';
dataConDefault = [];
dataExpDefault = [];
meanNullDefault = [];
meanAltDefault = [];
stdevDefault = [];
pNullDefault = [];
pAltDefault = [];
statPowerDesiredDefault = 0.8;
nSamplesConDefault = [];
nSamplesExpDefault = [];
nSamplesDefault = [];
nSamplesRatioDefault = [];
saveSheetFlagDefault = false;   % don't save values in a spreadsheet by default
saveMatFlagDefault = false;     % don't save values in a .mat file by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OutputType', outputTypeDefault, ...
    @(x) any(validatestring(x, validOutputTypes)));
addParameter(iP, 'TestType', testTypeDefault, ...
    @(x) assert(iscell(x) || any(validatestring(x, validTestTypes)), ...
                ['TestType must be either a cell array of character vectors', ...
                    'or one of ', print_cellstr(validOutputTypes, 'ToPrint', false)]));
addParameter(iP, 'DataCon', dataConDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['DataCon must be either a numeric array ', ...
                    'or a cell array of numeric arrays!']));
addParameter(iP, 'DataExp', dataExpDefault, ...
    @(x) assert(isnumeric(x) || iscellnumeric(x), ...
                ['DataExp must be either a numeric array ', ...
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
addParameter(iP, 'StatPowerDesired', statPowerDesiredDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesCon', nSamplesConDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesExp', nSamplesExpDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamples', nSamplesDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'NSamplesRatio', nSamplesRatioDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'SaveSheetFlag', saveSheetFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveMatFlag', saveMatFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
outputType = validatestring(iP.Results.OutputType, validOutputTypes);
testType = iP.Results.TestType;
dataCon = iP.Results.DataCon;
dataExp = iP.Results.DataExp;
meanNull = iP.Results.MeanNull;
meanAlt = iP.Results.MeanAlt;
stdev = iP.Results.Stdev;
pNull = iP.Results.PNull;
pAlt = iP.Results.PAlt;
statPowerDesired = iP.Results.StatPowerDesired;
nSamplesCon = iP.Results.NSamplesCon;
nSamplesExp = iP.Results.NSamplesExp;
nSamples = iP.Results.NSamples;
nSamplesRatio = iP.Results.NSamplesRatio;
saveSheetFlag = iP.Results.SaveSheetFlag;
saveMatFlag = iP.Results.SaveMatFlag;

% Keep unmatched arguments for the sampsizepwr() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Match the formats of the two vector sets
[dataCon, dataExp] = ...
    match_format_vector_sets(dataCon, dataExp, 'ForceCellOutputs', true);

% Force as column cell arrays
[testType, pNull] = argfun(@force_column_cell, testType, pNull);

% Force as column cell arrays with each scalar in a cell
[meanNull, meanAlt, stdev, pAlt, statPowerDesired, ...
    nSamplesCon, nSamplesExp, nSamples, nSamplesRatio] = ...
    argfun(@(x) force_column_cell(x, 'TreatVectorAsElement', false), ...
            meanNull, meanAlt, stdev, pAlt, statPowerDesired, ...
            nSamplesCon, nSamplesExp, nSamples, nSamplesRatio);

% Count the number of different values to test
[nVecs(1), nVecs(2), nVecs(3), nVecs(4), nVecs(5), nVecs(6), nVecs(7), ...
    nVecs(8), nVecs(9), nVecs(10), nVecs(11), nVecs(12), nVecs(13)] = ...
    argfun(@numel, testType, dataCon, dataExp, pNull, ...
            meanNull, meanAlt, stdev, pAlt, statPowerDesired, ...
            nSamplesCon, nSamplesExp, nSamples, nSamplesRatio);

% Check input lengths
if numel(unique(nVecs)) > 2
    fprintf('There are more than two input lengths greater than 1!\n');
    return
end
        
% Count the maximum number of values to test
nVectors = max(nVecs);

% Match the row counts of all inputs
[testType, dataCon, dataExp, meanNull, meanAlt, stdev, pNull, pAlt, ...
    statPowerDesired, nSamplesCon, nSamplesExp, nSamples, nSamplesRatio] = ...
    argfun(@(x) match_row_count(x, nVectors), ...
            testType, dataCon, dataExp, meanNull, meanAlt, stdev, pNull, pAlt, ...
            statPowerDesired, nSamplesCon, nSamplesExp, nSamples, nSamplesRatio);

%% Do the job
% Compute the sample size or power
[nSamplesNeeded, statPower, pNull, pAlt, ...
    nSamplesCon, nSamplesExp, statPowerDesired] = ...
    cellfun(@(a, b, c, d, e, f, g, h, i, j, k, l, m) ...
            compute_sampsizepwr_helper(a, b, c, d, e, f, ...
                                    g, h, i, j, k, l, m, otherArguments), ...
            testType, dataCon, dataExp, meanNull, meanAlt, ...
            stdev, pNull, pAlt, statPowerDesired, ...
            nSamplesCon, nSamplesExp, nSamples, nSamplesRatio, ...
            'UniformOutput', false);

%% Output results in a table
switch outputType
    case {'all', 'scalars'}
        % Convert back to numeric vectors
        [nSamplesNeeded, statPower, pAlt, ...
            nSamplesCon, nSamplesExp, statPowerDesired] = ...
            argfun(@cell2num, nSamplesNeeded, statPower, pAlt, ...
                nSamplesCon, nSamplesExp, statPowerDesired);

        % Place together in a table
        sampsizeTable = table(testType, dataCon, dataExp, ...
                    nSamplesCon, nSamplesExp, ...
                    pNull, pAlt, statPower, statPowerDesired, nSamplesNeeded);

        % Save the table as a .mat file
        if saveMatFlag 
            matPath = [create_time_stamp, '_', fileSuffix, '.mat'];
            save(matPath, 'sampsizeTable');
        end

        % Save the table as a spreadsheet file
        if saveSheetFlag 
            sheetPath = [create_time_stamp, '_', fileSuffix, '.csv'];
            writetable(sampsizeTable, sheetPath);
        end
        
        % Remove raw data if requested
        if strcmpi(outputType, 'scalars')
            % Compute scalar variables
            meanNull = extract_elements(pNull, 'specific', 'Index', 1);
            stdNull = extract_elements(pNull, 'specific', 'Index', 2);
            meanAlt = pAlt;

            % Add scalar variables
            sampsizeTable = ...
                addvars(sampsizeTable, meanNull, stdNull, meanAlt, ...
                        'Before', 'pNull');

            % Remove vector variables
            sampsizeTable = ...
                removevars(sampsizeTable, {'dataCon', 'dataExp', ...
                                            'pNull', 'pAlt'});
        end
        
        % Return as output
        output = sampsizeTable;
    case 'sampsize'
        output = cell2num(nSamplesNeeded);
    case 'power'
        output = cell2num(statPower);
    case 'pAlt'
        output = cell2num(pAlt);
    otherwise
        error('outputType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [nSamplesNeeded, statPower, pNull, pAlt, ...
                nSamplesCon, nSamplesExp, statPowerDesired] = ...
                compute_sampsizepwr_helper (testType, dataCon, dataExp, ...
                                meanNull, meanAlt, stdev, ...
                                pNull, pAlt, statPowerDesired, ...
                                nSamplesCon, nSamplesExp, ...
                                nSamples, nSamplesRatio, otherArguments)
%% Compute the statistical power for one set of data

% Decide on the number of samples for the control group
if isempty(nSamplesCon)
    if ~isempty(dataCon)
        % Count from data if provided
        nSamplesCon = count_samples(dataCon);
    elseif ~isempty(nSamples)
        % Just use number of samples
        nSamplesCon = nSamples;
    elseif ~isempty(nSamplesExp)
        % Use the same number of samples as the experimental group
        nSamplesCon = nSamplesExp;
    end
end

% Decide on the number of samples for the experimental group
if isempty(nSamplesExp)
    if ~isempty(dataExp)
        % Count from data if provided
        nSamplesExp = count_samples(dataExp);
    elseif ~isempty(nSamples)
        % Just use number of samples
        if ~isempty(nSamplesRatio)
            % Multiply by the ratio if provided
            nSamplesExp = nSamples * nSamplesRatio;
        else
            % The ratio is 1 by default
            nSamplesExp = nSamples;
        end
    elseif ~isempty(nSamplesCon)
        % Use the same number of samples as the control group
        nSamplesExp = nSamplesCon;
    end
end

% Decide on the sample size
%   Note: use the smaller of the two samples for sampsizepwr()
if isempty(nSamples)
    if ~isempty(nSamplesCon) || ~isempty(nSamplesExp)
        % Use the smaller number of samples
        nSamples = min([nSamplesCon, nSamplesExp]);
    else
        % Keep empty
    end
end

% Decide on the sample size ratio
if isempty(nSamplesRatio)
    if ~isempty(nSamplesCon) || ~isempty(nSamplesExp)
        % Compute the larger of the two samples
        nSamplesMax = max([nSamplesCon, nSamplesExp]);

        % Compute the sample size ratio
        nSamplesRatio = nSamplesMax ./ nSamples;
    else
        nSamplesRatio = 1;
    end
end

% Estimate a standard deviation if not provided
%   Note: this might be NaN
if isempty(stdev)
    if strcmpi(testType, 'z') || strcmpi(testType, 't')
        if ~isempty(dataExp) && ~isempty(dataCon)
            stdev = compute_stats(dataExp - dataCon, 'std');
        elseif ~isempty(dataExp)
            stdev = compute_stats(dataExp, 'std');
        else
            stdev = compute_stats(dataCon, 'std');
        end
    else
        stdNull = compute_stats(dataCon, 'std');
        stdAlt = compute_stats(dataExp, 'std');
        stdev = nanmean([stdNull, stdAlt]);
    end
end

% Estimate null hypothesis parameters if not provided
%   Note: this might be NaN
if isempty(pNull)
    pNull = estimate_dist_params(dataCon, testType, meanNull, stdev, 'null');
end


% Estimate alternative hypothesis parameters if not provided
if isempty(pAlt)
    if ~isempty(dataExp)
        if strcmpi(testType, 'z') || strcmpi(testType, 't')
            pAlt = estimate_dist_params(dataExp - dataCon, testType, ...
                                        meanAlt, stdev, 'alt');
        else
            pAlt = estimate_dist_params(dataExp, testType, ...
                                        meanAlt, stdev, 'alt');
        end
    elseif ~isempty(meanAlt)
        pAlt = meanAlt;
    elseif ~isempty(nSamples)
        pAlt = sampsizepwr(testType, pNull, [], statPowerDesired, nSamples, ...
                            'Ratio', nSamplesRatio, otherArguments{:});
    else
        pAlt = NaN;
    end
end

% Compute the statistical power of the test
if ~isempty(nSamples)
    statPower = sampsizepwr(testType, pNull, pAlt, [], nSamples, ...
                            'Ratio', nSamplesRatio, otherArguments{:});
else
    statPower = NaN;
end

% Compute the number of samples needed to reach that statistical power desired
nSamplesNeeded = sampsizepwr(testType, pNull, pAlt, statPowerDesired, [], ...
                            'Ratio', nSamplesRatio, otherArguments{:});

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function params = estimate_dist_params (data, testType, ...
                                        meanValue, stdValue, paramType)
%% Estimates parameters for the underlying distribution

% Compute the mean parameter value if not provided
if isempty(meanValue)
    switch testType
        case {'z', 't'}
            switch paramType
                case 'null'
                    meanValue = 0;
                case 'alt'
                    meanValue = compute_stats(data, 'mean');
                otherwise
                    error('paramType unrecognized!');
            end
        case 't2'
            meanValue = compute_stats(data, 'mean');
        case 'var'
            % Compute the sample variance
            meanValue = compute_stats(data, 'var');
        case 'p'
            % Estimate the proportion of ones from a Bernoulli trial
            meanValue = sum(data)/numel(data);
        otherwise
            error('testType unrecognized!');
    end
end

switch testType
    case {'z', 't', 't2'}
        % Put the normal distribution parameters together
        switch paramType
            case 'null'
                params = [meanValue, stdValue];
            case 'alt'
                params = meanValue;
            otherwise
                error('paramType unrecognized!');
         end 
    case {'var', 'p'}
        params = meanValue;
    otherwise
        error('testType unrecognized!');
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

if ~any(isemptycell(dataCon))
    dataCon = cellfun(@(x, y) x - y, dataCon, dataCon, ...
                        'UniformOutput', false);
end
if ~any(isemptycell(dataCon)) && ~any(isemptycell(dataExp))
    dataExp = cellfun(@(x, y) x - y, dataExp, dataCon, ...
                        'UniformOutput', false);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%