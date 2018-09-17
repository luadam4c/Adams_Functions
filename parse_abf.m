function [data, siUs, abfParams] = parse_abf(fileName, varargin)
%% Loads and parses an abf file
% Usage: [data, siUs, abfParams] = parse_abf(fileName, varargin)
%
% Outputs:
%       data        - full data
%       siUs        - sampling interval in microseconds
%       abfParams   - TODO
% Arguments:
%       fileName    - file name could be either the full path or 
%                       a relative path in current directory
%                       .abf is not needed (e.g. 'B20160908_0004')
%                   must be a string scalar or a character vector
%       varargin    - 'Verbose': whether to output parsed results
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%
% Requires:
%       /home/Matlab/Adams_Functions/construct_abffilename.m
%       /home/Matlab/Downloaded_Functions/abf2load.m or abfload.m
%       /home/Matlab/Brians_Functions/identify_channels.m
%
% Used by:
%       cd/plot_traces_abf.m
%
% File history: 
% 2018-09-17 - Moved from plot_traces_abf.m
% 2018-09-17 - Added 'Verbose' as a parameter

%% Hard-coded constants
US_PER_MS = 1e3;            % number of microseconds per millisecond
US_PER_S = 1e6;             % number of microseconds per second

%% Default values for optional arguments
verboseDefault = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Add directories to search path for required functions
if ~isdeployed
    if exist('/home/Matlab/', 'dir') == 7
        functionsdirectory = '/home/Matlab/';
    elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
        functionsdirectory = '/scratch/al4ng/Matlab/';
    else
        error('Valid functionsdirectory does not exist!');
    end
    addpath(fullfile(functionsdirectory, '/Downloaded_Functions/'));
                                            % for abf2load.m or abfload.m
    addpath(fullfile(functionsdirectory, '/Brians_Functions/'));
                                            % for identify_channels.m
end

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
addRequired(iP, 'fileName', ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
                                                % introduced after R2016b
% Add parameter-value pairs to the Input Parser
addParameter(iP, 'Verbose', verboseDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, fileName, varargin{:});
verbose = iP.Results.Verbose;

%% Do the job
% Create the full path to .abf file robustly
abfFullFileName = construct_abffilename(fileName);

% Load abf file, si is in us
if exist('abf2load', 'file') == 2
    try
        [data, siUs] = abf2load(abfFullFileName);
    catch ME
        printf('The file %s cannot be read!\n', abfFullFileName);
        rethrow(ME)
        return
    end
elseif exist('abfload', 'file') == 2
    [data, siUs] = abfload(abfFullFileName);
end

% Find data dimensions and make sure it is <= 3
ndim = ndims(data);        % number of dimensions in data
if ndim > 3
    error('Cannot parse data with more than 3 dimensions!');
end

% Query data dimensions
nSamples = size(data, 1);          % number of samples
nChannels = size(data, 2);         % number of channels
if ndim == 3
    nSweeps = size(data, 3);       % number of sweeps
end

% Identify proper channel units and labels
[~, channelUnits, channelLabels] = identify_channels(data);

% Convert sampling interval to other units
siMs = siUs / US_PER_MS;
siSeconds = siUs / US_PER_S;

% Store in abfParams
abfParams.channelUnits = channelUnits;
abfParams.channelLabels = channelLabels;
abfParams.ndim = ndim;
abfParams.nSamples = nSamples;
abfParams.nChannels = nChannels;
abfParams.nSweeps = nSweeps;

% Write results to standard output
if verbose
    fprintf('Sampling Interval (us) = %g\n', siUs);
    fprintf('Sampling Interval (ms) = %g\n', siMs);
    fprintf('Sampling Interval (s) = %g\n', siSeconds);
    fprintf('Channel Units = %s\n', strjoin(channelUnits, ', '));
    fprintf('Channel Labels = %s\n', strjoin(channelLabels, ', '));
    fprintf('Number of data dimensions = %d\n', ndim);
    fprintf('Number of samples = %d\n', nSamples);
    fprintf('Number of channels = %d\n', nChannels);
    fprintf('Number of sweeps = %d\n\n', nSweeps);
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

