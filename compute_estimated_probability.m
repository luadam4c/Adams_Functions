function [estProb, lowerCI, upperCI] = compute_estimated_probability(posCounts, totalCounts)
%% Fits a mixed-effects logistic regression model to positive and total counts across many subjects
% Usage: [estProb, lowerCI, upperCI] = compute_estimated_probability(posCounts, totalCounts)
% Explanation:
%       This function reformats count data into binary responses and fits 
%       a binomial generalized linear mixed-effects model (GLME) to estimate
%       the overall probability, lower confidence interval, and upper confidence interval.
%
% Outputs:
%       estProb     - estimated probability calculated from the model coefficient
%                   specified as a numeric scalar
%       lowerCI     - lower confidence interval of the estimated probability
%                   specified as a numeric scalar
%       upperCI     - upper confidence interval of the estimated probability
%                   specified as a numeric scalar
%
% Arguments:
%       posCounts   - vector of number of positive counts per subject
%                   must be a numeric vector
%       totalCounts - vector of totalCounts counts (posCounts+neg) per subject
%                   must be a numeric vector
%
% Requires:
%       fitglme() from the Statistics and Machine Learning Toolbox
%
% Used by:
%       scAAV/analyze_qupath_figure2.m
%       scAAV/analyze_qupath_ddio.m

% File History:
% 2020-05-29 Created by Jeff Moore as stats_cell_counts_melogit.m
% 2026-03-04 Annotated by Gemini and renamed by Adam
% 2026-03-04 Added data validation by Gemini

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Preparation
% Filter out invalid data points (NaNs or zero total counts)
isValidData = ~isnan(posCounts) & ~isnan(totalCounts) & (totalCounts > 0);
posCounts = posCounts(isValidData);
totalCounts = totalCounts(isValidData);

% Return NaNs if there are no valid data points left
if isempty(posCounts)
    estProb = NaN;
    lowerCI = NaN;
    upperCI = NaN;
    return;
end

% Generate a cell array of binary responses (1 for positive, 0 for negative) for each subject
binaryResponseCell = arrayfun(@(p, t) [ones(p, 1); zeros(t - p, 1)], ...
                                    posCounts(:), totalCounts(:), 'UniformOutput', false);

% Generate a cell array of repeated subject IDs matching the totalCounts counts for each subject
subjectIdCell = arrayfun(@(t, id) ones(t, 1) * id, ...
                            totalCounts(:), (1:numel(totalCounts))', 'UniformOutput', false);

% Concatenate the binary response cell array into a single numerical vector
y = cell2mat(binaryResponseCell);

% Concatenate the subject ID cell array into a single numerical vector
subjectId = cell2mat(subjectIdCell);

%% Do the job
% Create a table from the response and ID vectors
tbl = table(y, subjectId);

% Define the mixed-effects formula
formula = 'y ~ 1 + (1|subjectId)';

% Fit the binomial GLME model
glme = fitglme(tbl, formula, 'Distribution', 'Binomial');

% Extract the estimated coefficients and their confidence intervals
coef = glme.Coefficients.Estimate;
coefLci = glme.Coefficients.Lower;
coefUci = glme.Coefficients.Upper;

% Convert the logit coefficients back to probabilities
estProb = glmval(coef, 1, 'logit', 'constant', 'off');
lowerCI = glmval(coefLci, 1, 'logit', 'constant', 'off');
upperCI = glmval(coefUci, 1, 'logit', 'constant', 'off');

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%