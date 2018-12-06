function product = outer_product (X, Y)
%% Returns the outer product of two vectors (could be cell arrays)
% Usage: product = outer_product (X, Y)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       product     - TODO: Description of product
%                   specified as a TODO
% Arguments:
%       reqarg1     - TODO: Description of reqarg1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Used by:
%       cd/m3ha_organize_sweep_indices.m

% File History:
% 2018-12-06 Created by Adam Lu
% 

%% Hard-coded parameters

%% Default values for optional arguments
param1Default   = [];                   % default TODO: Description of param1

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 2
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
end

% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add required inputs to the Input Parser
% addRequired(iP, 'reqarg1', ...                  % TODO: Description of reqarg1
%     % TODO: validation function %);

% Add parameter-value pairs to the Input Parser
% addParameter(iP, 'param1', param1Default, ...
%     % TODO: validation function %);

% Read from the Input Parser
% parse(iP, reqarg1, varargin{:});
% param1 = iP.Results.param1;

% Check relationships between arguments
% TODO: X and Y must have the same lengths

%% Preparation
% TODO

%% Do the job
if iscell(X) && iscell(Y)
    % Count the number of elements
    nX = numel(X);
    nY = numel(Y);

    % Force X as a column vector
    X = reshape(X, [], 1);

    % Force Y as a row vector
    Y = reshape(Y, 1, []);

    % Repeat X nY times
    Xrepeated = repmat(X, [1, nY]);

    % Repeat Y nX times
    Yrepeated = repmat(Y, [nX, 1]);

    % Join the strings
    product = cellfun(@(x, y) [x, ', ', y], Xrepeated, Yrepeated, ...
                    'UniformOutput', false);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%