function label = create_label_from_sequence (integers, varargin)
%% Creates a single label from a sequence of integers with an optional prefix or suffix
% Usage: label = create_label_from_sequence (integers, varargin)
% Explanation:
%       TODO
% Example(s):
%       label = create_label_from_sequence([-1, 1, 2, 5, 6, 7])
%       label = create_label_from_sequence([5; 3; 1; 8])
%       label = create_label_from_sequence(magic(3))
%       label = create_label_from_sequence([-1, 1, 2, 5], 'Suffix', ' Mississippi')
%       label = create_label_from_sequence([-1, 1, 2, 5], 'Prefix', 'Traces ')
% Outputs:
%       label       - label created
%                   specified as a cell array of character vectors
% Arguments:
%       integers    - integers
%                   must be an array of one the following types:
%                       'numeric', 'logical', 'datetime', 'duration'
%       varargin    - 'Prefix': string to place before each number
%                   must be a string scalar or a character vector
%                   default == ''
%                   - 'Suffix': string to place after each number
%                   must be a string scalar or a character vector
%                   default == ''
%
%
% Requires:
%       cd/create_error_for_nargin.m
%       cd/convert_to_char.m
%
% Used by:
%       cd/plot_measures.m
%       cd/plot_swd_psth.m
%       /home/Matlab/plethR01/plethR01_analyze.m

% File History:
% 2019-08-07 Created by Adam Lu
% TODO: Add 'DelimiterThrough' and 'DelimiterSeparate'
% 

%% Hard-coded parameters
delimiterThrough = '-';
delimiterSeparate = ',';

%% Default values for optional arguments
prefixDefault = '';     % no string to place before each number by default
suffixDefault = '';     % no string to place after each number by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
addRequired(iP, 'integers', ...
    @(x) validateattributes(x, {'numeric', 'logical', ...
                                'datetime', 'duration'}, {'3d'}));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Suffix', suffixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));

% Read from the Input Parser
parse(iP, integers, varargin{:});
prefix = iP.Results.Prefix;
suffix = iP.Results.Suffix;

%% Preparation
% Force the integers as a column vector
integers = force_column_vector(integers, 'ToLinearize', true);

%% Do the job
if isempty(integers)
    integersStr = '';
else
    % Sort the unique integers in ascending order
    uniqueIntegers = unique(integers, 'sorted');

    % Count the number of unique integers
    nNums = numel(uniqueIntegers);

    if nNums == 1
        integersStr = num2str(uniqueIntegers);
    else
        % Test whether each integer has a consecutive partner following it
        hasNext = [(uniqueIntegers(2:end) - uniqueIntegers(1:end-1) == 1); false];

        % Test whether each integer has a consecutive partner preceding it
        hasPrev = [false; hasNext(1:end-1)];

        % Decide whether each integer is a first of a run
        isFirstOfRun = hasNext & ~hasPrev;

        % Decide whether each integer is a last of a run
        isLastOfRun = ~hasNext & hasPrev;

        % Decide whether each integer is a standalone
        isStandalone = ~hasNext & ~hasPrev;

        % Create snippets for each integer
        snippets = cell(nNums, 1);
        for iNum = 1:nNums
            % Get the current integer
            intThis = uniqueIntegers(iNum);

            % Create snippet
            if iNum == 1
                snippets{iNum} = num2str(intThis);
            elseif isStandalone(iNum) || isFirstOfRun(iNum)
                snippets{iNum} = [delimiterSeparate, num2str(intThis)];
            elseif isLastOfRun(iNum)
                snippets{iNum} = [delimiterThrough, num2str(intThis)];
            else
                snippets{iNum} = '';
            end
        end

        % Concatenate snippets together
        integersStr = strjoin(snippets, '');
    end
end

% Create the final label
label = [prefix, integersStr, suffix];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
