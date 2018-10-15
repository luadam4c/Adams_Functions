function nSigFig = sigfig(number, varargin)
%% Returns the number of significant figures from a number (numeric or string)
% Usage: nSigFig = sigfig(number, varargin)
% Explanation:
%       Returns 0 for 0
%       It seems like the maximum number of significant figures
%           possible is 16. And it still doesn't do well 
%           with very small numbers (2.3e-400) or very large numbers (2.3e400)
%                                                           -AL 2/14/2018
%
%       https://www.mathworks.com/matlabcentral/answers/52068-machine-precision-problem-in-code
%       For floating point numbers, the number of digits is not well defined, 
%       because the displayed output use the base 10 while the internal storage 
%       uses the base 2. 
%       The conversion of the bases effects the number of digits due to the 
%       limited precision of the internal representation of doubles. 
%       In consequence, there cannot be a program, which replies the 
%       "number of digits after the decimal point" for "any rational number".
% Outputs:
%       nSigFig     - number of significant figures
%                   specified as an integer
% Arguments:    
%       number      - the number of interest
%                   must be a numeric scalar or a character vector or a string
%       varargin    - TODO
%
% Used by:    
%
% File History:
% 2018-02-14 Created by Adam Lu
% 2018-02-15 Fixed the infinite loop when the argument is 0
% TODO: Easy: the case where the input is a character vector or a string
% TODO: Hard: Can bi2de and de2bi solve the precision problem?

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
addRequired(iP, 'number', @isnumeric);
    % TODO: Change validation function once the character array version is implemented

% Read from the Input Parser
parse(iP, number, varargin{:});

%% Do the job
% If the number is within eps of 0, return 0
if number < eps(0)
    nSigFig = 0;
    return
end

% Find the first place of the first nonzero digit
firstNonzeroPlace = -floor(log10(abs(number)));

% Find the place of the last nonzero digit after the decimal point
temp = abs(number);
lastAfterPoint = 0;
difference1 = abs(floor(temp) - temp);
difference2 = abs(ceil(temp) - temp);
precision = eps(temp);
while difference1 >= 2*precision && difference2 >= 2*precision
            % above needed for very small numbers
    lastAfterPoint = lastAfterPoint + 1;
    temp = abs(number * 10^lastAfterPoint);
    difference1 = abs(floor(temp) - temp);
    difference2 = abs(ceil(temp) - temp);
    precision = eps(temp);
end

% Find the place before the last nonzero digit before the decimal point
temp = abs(number);
beforeLastBeforePoint = 0;
difference1 = abs(floor(temp) - temp);
difference2 = abs(ceil(temp) - temp);
precision = eps(temp);
while difference1 < 2*precision || difference2 < 2*precision
            % above needed for very large numbers
    beforeLastBeforePoint = beforeLastBeforePoint - 1;
    temp = abs(number * 10^beforeLastBeforePoint);
    difference1 = abs(floor(temp) - temp);
    difference2 = abs(ceil(temp) - temp);
    precision = eps(temp);
end
lastBeforePoint = beforeLastBeforePoint + 1;

% The last decimal place is lastAfterPoint if it is greater than zero
%   Otherwise, the number is an integer so we can apply num2str
%   and remove the trailing zeros
if lastAfterPoint > 0            % number is not an integer
    % The last nonzero decimal place is 
    %   the place of the last nonzero digit after the decimal point
    lastNonzeroPlace = lastAfterPoint;
elseif beforeLastBeforePoint < 0 % number is an integer
    % The last nonzero decimal place is 
    %   the place of the last nonzero digit before the decimal point
    lastNonzeroPlace = lastBeforePoint;
else
    error('Something wrong with code!');
end

% Compute the number of significant figures
nSigFig = lastNonzeroPlace - firstNonzeroPlace + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

% Print integer as string
% This doesn't work in Matlab for decimal values 
%   because floating point is not exact
numberStr = num2str(number);

% Find the first decimal place
firstNonzeroPlace = -floor(log10(number));

% Find the last decimal place
place = firstNonzeroPlace;
temp = abs(number);
while temp ~= 0
    temp = temp - round(temp, place);
    place = place + 1;
end
lastNonzeroPlace = place - 1;

while floor(temp) ~= temp


numberStr = num2str(number, maxSigFig);

    % Print number as string with the maximum number of significant figures
    numberStr = num2str(number);
%    numberStr = num2str(number, sprintf('%%%d.%df', maxSigFig, maxSigFig));

    % Remove decimal point
    numberStr = strrep(numberStr, '.', '');

    % Remove negative sign
    numberStr = strrep(numberStr, '-', '');

    % Find the first non-zero digit
    first = find(numberStr ~= '0', 1, 'first');

    % Find the last non-zero digit
    last = find(numberStr ~= '0', 1, 'last');

    % Remove leading and trailing zeros
    numberStr = numberStr(first:last);

    % Get the number of significant figures
    nSigFig = length(numberStr);

% TODO: Make 'MaxSigFig' an optional argument
%% Default values for optional arguments
maxSigFigDefault = 15;      % maximum number of significant figures by default
% TODO: Change this
maxSigFig = maxSigFigDefault;

%       /home/Matlab/minEASE/minEASE.m

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%