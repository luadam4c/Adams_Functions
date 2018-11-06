function equationLatex = create_latex_string (equation, varargin)
%% Creates a LaTeX string from an equation used for fitting
% Usage: equationLatex = create_latex_string (equation, varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       equationLatex   - LaTeX form of equation
%                   specified as a character vector
% Arguments:
%       equation    - equation form
%                   must be a character vector
%
% Used by:
%       cd/plot_cfit_pulse_response.m

% File History:
% 2018-11-01 Created by Adam Lu
% 

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
addRequired(iP, 'equation', ...
    @(x) validateattributes(x, {'char'}, {'2d'}));

% Read from the Input Parser
parse(iP, equation, varargin{:});

%% Do the job
% Replace the exponentials
equation = replace(equation, '*exp(', '*e^{');
equation = replace(equation, 'exp(', 'e^{');

% Remove the * if unnecessary
equation = replace(equation, '*(', '(');

% Replace <= with \leq
equation = replace(equation, '<=', '\leq');

% Count the number of characters
nChars = numel(equation);

% Replace the matching ')' for each '{' with '}'
nLeft = 0;          % the number of left parenthesis encountered so far
toDetect = false;   % whether a '{' was encountered
for iChar = 1:nChars
    % Get this character 
    thisChar = equation(iChar);

    % Test the character
    if thisChar == '{'
        % Start detection
        toDetect = true;

        % Increment nLeft and start detection
        nLeft = nLeft + 1;
    elseif thisChar == '(' && toDetect
        % Increment nLeft and start detection
        nLeft = nLeft + 1;
    elseif thisChar == ')' && toDetect
        % Decrement nLeft
        nLeft = nLeft - 1;

        % If there is no nLeft left, then we need to replace the character
        %   and end detection
        if nLeft == 0
            % Replace with '}'
            equation(iChar) = '}';

            % Set detect flag to be false
            toDetect = false;
        end
    end
end

% Prepend and append with '$$'
equation = ['$$', equation, '$$'];

%% Output results
equationLatex = equation;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%