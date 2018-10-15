function [fitResults] = parse_fitObject (fitObject, varargin)
%% Extract information from a cfit or sfit object
% Usage: [fitResults] = parse_fitObject (fitObject, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       fitResults      - parsed results, including the fields:
%                           nCoeffs
%                           equationStr
%                           fitFormula
%                           coeffNames
%                           coeffValues
%                           confInt
%                       specified as a structure
%                       - TODO: a table containing the parsed results
%                       specified as a table
% Arguments:    
%       fitObject   - a cfit or sfit object returned by fit()
%                   must be a cfit or sfit object
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:    
%       cd/fit_2exp.m

% File History:
% 2018-10-10 Created by Adam Lu, adapted some code from 
%               https://www.mathworks.com/matlabcentral/answers/99155-how-do-i-extract-the-output-equation-from-a-cfit-object-in-curve-fitting-toolbox-2-0-r2009a
% 

%% Hard-coded parameters
nSigFig = 2;    % maximum number of significant figures in the equation string

%% Default values for optional arguments

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

% Add required inputs to the Input Parser
addRequired(iP, 'fitObject');
    % TODO: validation function %);

% Add parameter-value pairs to the Input Parser

% Read from the Input Parser
parse(iP, fitObject, varargin{:});

% Check relationships between arguments
% TODO

%% Extract information from the fit object
% Apply built-in functions
fitFormula = formula(fitObject);
coeffNames = coeffnames(fitObject);
coeffValues = coeffvalues(fitObject);
confInt = confint(fitObject);

% Count the number of coefficients
nCoeffs = numel(coeffNames);

%% Generate equation string
% Initialize with the formula
equationStr = fitFormula;

% Replace each coefficient with its value
for iCoeff = 1:nCoeffs
    % Get the parameter name
    param = coeffNames{iCoeff};

    % Get the parameter value
    value = coeffValues(iCoeff);

    % Replace all occurrences of param with its value
    equationStr = strrep(equationStr, param, num2str(value, nSigFig));
end

%% Store results in output
fitResults.nCoeffs = nCoeffs;
fitResults.equationStr = equationStr;
fitResults.fitFormula = fitFormula;
fitResults.coeffNames = coeffNames;
fitResults.coeffValues = coeffValues;
fitResults.confInt = confInt;

% TODO
% fitResults = table(nCoeffs, equationStr, ...
%                     fitFormula, coeffNames, coeffValues, confInt);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Count the number of characters in the parameter name
nChars = length(param);

% Find the starting index of the first matching parameter 
%   within the equation string
idxStart = regexp(equationStr, param);

% Substitute matching parameters with its value
%   while there is still a matching parameter
while ~isempty(idxStart) 
    % Substitute the parameter by its value
    equationStr = [equationStr(1:idxStart-1), ...
                    num2str(value), ...
                    equationStr(idxStart+nChars:end)];

    % Find the starting index of the next matching parameter 
    %   within the equation string
    idxStart = regexp(equationStr, param);
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%