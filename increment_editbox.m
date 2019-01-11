function newStr = increment_editbox(hEditBox, minNo, maxNo, incr, strings, varargin)
%% Increment or decrement editbox value based on direction
% Usage: newStr = increment_editbox(hEditBox, minNo, maxNo, incr, strings, varargin)
% Explanation:
%   Note: (1) strings is a cell array of character arrays or strings
%         (2) direction can be 'up' or 'down' based on sign of incr
%         (3) Goes circular with strings and linear if strings is empty
%               Order of 'up' if linear: minNo -> maxNo
%               Order of 'up' if circular: minNo -> maxNo -> strings -> minNo
%   TODO
% Outputs:    
%   TODO
% Arguments:    
%   TODO
%
%
% Used by:    
%       /home/Matlab/minEASE/gui_examine_events.m
%
% File History:
% 2017-06-09 Moved from gui_examine_events.m
% 2017-06-10 Now uses built-in strcmp() instead of find_in_strings.m
% 2017-06-12 Now allows incr to be nonzero and sets the direction 
%               according to the sign of incr
% 2017-07-31 Added newStr and return it as an output
% 2017-08-03 Fixed bug by adding newStr = prevStr;
% 2017-10-17 Added parameter 'RankNo' to allow incrementing according to rank
% 2017-10-18 Added at_lower_limit() & at_upper_limit()
% 2018-02-03 Added parameter 'Operation' to allow incrementing by multiples
%

possibleOperations = {'add', 'multiply'};

valueNoDefault = [];                    % default value of each number
operationDefault = 'add';               % default operation for incrementation

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% TODO: Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = 'increment_editbox';

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ValueNo', valueNoDefault);
addParameter(iP, 'Operation', operationDefault, ...
    @(x) any(validatestring(x, possibleOperations)));

% Read from the Input Parser
parse(iP, varargin{:});
valueNo = iP.Results.ValueNo;                   % value of each number
operation = validatestring(iP.Results.Operation, possibleOperations);

% Extract number of strings
nStrs = numel(strings);                         % number of strings

% Extract this string from edit box
prevStr = get(hEditBox, 'String');              % previous string

% Convert to number
%   Note: NaN is returned if not a number
prevNo = str2double(prevStr);

% Check direction
if incr > 0
    direction = 'up';
elseif incr < 0
    direction = 'down';
end

% Get absolute value of increment
incrAbs = abs(incr);

% If valueNo is used, build full vector of numbers
%   and rank the numbers by value and find the corresponding increment of rank
if ~isempty(valueNo)
    nNo = length(valueNo);                      % number of numbers to rank
    if nNo == 1 || minNo == maxNo
        allNo = minNo;

        % Find the corresponding increment of rank
        incrRank = incrAbs;
    else
        spacing = (maxNo - minNo)/(nNo - 1);
        allNo = minNo:spacing:maxNo;

        % Find the corresponding increment of rank
        incrRank = round(incrAbs/spacing);
    end
    [~, origInd] = sort(valueNo);
    
    % Give a rank to each event
    rankNo(origInd) = (1:nNo)';
else
    nNo = [];
    allNo = [];
    rankNo = [];
    incrRank = [];
end

if nStrs == 0                                   % no strings exist, linear
    switch direction
    case 'down'
        if isnan(prevNo)
            error('The string %s is invalid!', prevStr);
        elseif at_lower_limit(prevNo, minNo, incrAbs, allNo, rankNo)
            % Do nothing
            newStr = prevStr;
        else                                    % not at lower limit
            % Decrement number
            newStr = decrement(prevNo, incrAbs, nNo, allNo, ...
                                rankNo, incrRank, operation);
        end       
    case 'up'
        if isnan(prevNo)
            error('The string %s is invalid!', prevStr);
        elseif at_upper_limit(prevNo, maxNo, incrAbs, allNo, rankNo, nNo)
            % Do nothing
            newStr = prevStr;
        else                                    % not at upper limit
            % Increment number
            newStr = increment(prevNo, incrAbs, nNo, allNo, ...
                                rankNo, incrRank, operation);
        end
    otherwise
        error('direction undefined!');
    end

elseif nStrs > 0                                % strings exist, circular
    % Check if prevStr is in strings
    %   Note: (1) An empty matrix is returned if not found
    %         (2) The 'exact' search mode is necessary to include the case
    %             where there are empty strings in strings
    idx = find(strcmp(prevStr, strings));

    switch direction
    case 'down'
        if idx == 1                         % first of strings
            % Set to maximum number or number with maximum value
            if ~isempty(rankNo)
                newStr = num2str(allNo(rankNo == nNo));
            else
                newStr = num2str(maxNo);
            end
        elseif ~isempty(idx)                % in strings
            % Set to previous string
            newStr = strings{idx - 1};
        elseif isnan(prevNo)
            error('The string %s is invalid!', prevStr);
        elseif at_lower_limit(prevNo, minNo, incrAbs, allNo, rankNo)
            % Set to last of strings
            newStr = strings{nStrs};
        else                                % not at lower limit
            % Decrement number
            newStr = decrement(prevNo, incrAbs, nNo, allNo, ...
                                rankNo, incrRank, operation);
        end
    case 'up'
        if idx == nStrs                     % last of strings
            % Set to minimum number or number with maximum value
            if ~isempty(rankNo)
                newStr = num2str(allNo(rankNo == 1));
            else
                newStr = num2str(1);
            end
        elseif ~isempty(idx)                % in strings
            % Set to next string
            newStr = strings{idx + 1};
        elseif isnan(prevNo)
            error('The string %s is invalid!', prevStr);
        elseif at_upper_limit(prevNo, maxNo, incrAbs, allNo, rankNo, nNo)
            % Set string to first of strings
            newStr = strings{1};
        else                                % not at upper limit
            % Increment number
            newStr = increment(prevNo, incrAbs, nNo, allNo, ...
                                rankNo, incrRank, operation);
        end
    otherwise
        error('direction undefined!');
    end
end

% Update string in edit box
set(hEditBox, 'String', newStr);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newStr = decrement(prevNo, incrAbs, nNo, allNo, rankNo, incrRank, operation)
%% Decrement number by rank or not

if ~isempty(rankNo)
    % Find the rank of the previous number
    prevRank = rankNo(allNo == prevNo);

    % Find the next rank smaller than the previous rank
    newRank = prevRank - incrRank;
    if newRank >= 1
        newNo = allNo(rankNo == newRank);
    else
        newNo = allNo(rankNo == 1);
    end
else
    % Just decrement the number
    switch operation
    case 'add'
        newNo = prevNo - incrAbs;    
    case 'multiply'
        newNo = prevNo / incrAbs;
    otherwise
        error('Incrementing operation undefined! Error with code!');
    end
end

% Generate new string
newStr = num2str(newNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function newStr = increment(prevNo, incrAbs, nNo, allNo, rankNo, incrRank, operation)
%% Increment number by rank or not

if ~isempty(rankNo)
    % Find the rank of the previous number
    prevRank = rankNo(allNo == prevNo);

    % Find the next rank greater than the previous rank
    newRank = prevRank + incrRank;
    if newRank <= nNo
        newNo = allNo(rankNo == newRank);
    else
        newNo = allNo(rankNo == nNo);
    end   
else
    % Just increment the number  
    switch operation
    case 'add'
        newNo = prevNo + incrAbs;    
    case 'multiply'
        newNo = prevNo * incrAbs;
    otherwise
        error('Incrementing operation undefined! Error with code!');
    end
end

% Generate new string
newStr = num2str(newNo);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function atLowerLimit = at_lower_limit(prevNo, minNo, incrAbs, allNo, rankNo)
%% Whether previous number is at lower limit

if ~isempty(rankNo)
    % Find the rank of the previous number
    prevRank = rankNo(allNo == prevNo);

    % At lower limit if the rank is 1
    atLowerLimit = prevRank == 1;
else
    % At lower limit if the next number would be smaller than the minimum
    atLowerLimit = prevNo - incrAbs < minNo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function atUpperLimit = at_upper_limit(prevNo, maxNo, incrAbs, allNo, rankNo, nNo)
%% Whether previous number is at upper limit

if ~isempty(rankNo)
    % Find the rank of the previous number
    prevRank = rankNo(allNo == prevNo);

    % At upper limit if the rank is 1
    atUpperLimit = prevRank == nNo;
else
    % At upper limit if the next number would be larger than the maximum
    atUpperLimit = prevNo + incrAbs > maxNo;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

idx = find_in_strings(prevStr, strings, 'SearchMode', 'exact');

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
