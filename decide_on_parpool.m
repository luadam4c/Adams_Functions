function poolObj = decide_on_parpool (varargin)
%% Creates or modifies a parallel pool object
% Usage: poolObj = decide_on_parpool (varargin)
% Explanation:
%       TODO
%
% Example(s):
%       poolObj = decide_on_parpool;
%       poolObj = decide_on_parpool([5, 10]);
%       poolObj = decide_on_parpool('MaxNumWorkers', Inf);
%       poolObj = decide_on_parpool('RenewParpoolFlag', true);
%
% Outputs:
%       poolObj     - parallel Pool object decided on
%                   specified as a parallel Pool object
%
% Arguments:
%       poolSize    - (opt) range of number of workers for the parallel pool
%                   must be a 1 or 2 element positive integer vector
%                   default == [1, maxNumWorkers]
%       varargin    - 'ProfileName': profile name for the parallel pool
%                   must be a string scalar or a character vector
%                   default == 'local'
%                   - 'RenewParpoolFlag': whether to refresh parallel pool
%                                           every batch
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNumWorkers': maximum number of workers for parfor
%                                   set to 0 if parfor not to be used
%                                   set to Inf to use maximum number of workers
%                   must be a nonnegative integer or Inf
%                   default == Inf
%                   - Any other parameter-value pair for parpool()
%
% Requires:
%       cd/struct2arglist.m
%       cd/is_in_parallel.m
%
% Used by:
%       cd/run_neuron.m

% File History:
% 2019-10-31 Created by Adam Lu
% 

%% Default values for optional arguments
poolSizeDefault = [];           % set later
profileNameDefault = 'local';
renewParpoolFlagDefault = false;% don't refresh parallel pool every batch
maxNumWorkersDefault = Inf;     % no limits by default

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;
iP.KeepUnmatched = true;                        % allow extraneous options

% Add optional inputs to the Input Parser
addOptional(iP, 'poolSize', poolSizeDefault);

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'ProfileName', profileNameDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'RenewParpoolFlag', renewParpoolFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNumWorkers', maxNumWorkersDefault, ...
    @(x) assert(isinf(x) || isscalar(x) && isaninteger(x) && x >= 0, ...
                'MaxNumWorkers must be Inf or a nonnegative integer!'));

% Read from the Input Parser
parse(iP, varargin{:});
poolSize = iP.Results.poolSize;
profileName = iP.Results.ProfileName;
renewParpoolFlag = iP.Results.RenewParpoolFlag;
maxNumWorkers = iP.Results.MaxNumWorkers;

% Keep unmatched arguments for the parpool() function
otherArguments = struct2arglist(iP.Unmatched);

%% Preparation
% Decid on the poolsize
if isempty(poolSize)
    poolSize = [1, maxNumWorkers];
else
    if numel(poolSize) > 1
        poolSize(2) = min(poolSize(2), maxNumWorkers);
    else
        poolSize = min(poolSize, maxNumWorkers);
    end
end

%% Do the job
% Get current parallel pool object without creating a new one
poolObj = gcp('nocreate');

% Return if in parallel loop
if is_in_parallel
    return
end

% Create a default parallel pool object if it doesn't exist
if isempty(poolObj)
    poolObj = parpool(profileName, poolSize, otherArguments{:});
end

% Count the number of workers in the current parallel pool object
oldNumWorkers = poolObj.NumWorkers;

% If the number of workers is out of range, 
%   or if requested, recreate the parallel pool object
if renewParpoolFlag || ...
        oldNumWorkers < poolSize(1) || oldNumWorkers > poolSize(end)
    delete(poolObj);
    poolObj = parpool(profileName, poolSize, otherArguments{:});
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%