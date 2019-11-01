function output = run_neuron (hocFile, varargin)
%% Runs NEURON using a .hoc file and a cell array of simulation commands/files
% Usage: output = run_neuron (hocFile, varargin)
% Explanation:
%       TODO
%
% Example(s):
%       output = run_neuron(hocFile, 'SimCommands', simCommands);
%
% Outputs:
%       output      - a table with the following fields:
%                       simStatus - simulation returned status
%                       simStdOut - simulation standard output
%                       hasError  - whether an error was encountered
%                   specified as a 2d table
%
% Arguments:    
%       hocFile     - .hoc file containing functions to use by NEURON
%                   must be a valid file
%       varargin    - 'SimCommands': simulation commands written in hoc
%                       can be either the actual commands or file names
%                   must be a cell array of strings or character vectors
%                   default == {}
%                   - 'SimNumbers': simulation numbers to run
%                       Note: this corresponds to the index in SimCommands
%                   must be empty or a positive integer vector
%                   default == []
%                   - 'OpenGuiFlag': whether to open GUI for NEURON
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'RunCommand': the command used to run NEURON
%                   must be a string scalar or a character vector
%                   default == 'nrngui' if OpenGui is true, 'nrniv' otherwise
%                   - 'Prefix': prefix to prepend to file names
%                   must be a character array
%                   default == ''
%                   - 'OutFolder': the directory where outputs will be placed
%                   must be a string scalar or a character vector
%                   default == pwd
%                   - 'DebugFlag': whether debugging
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'OnHpcFlag': whether on a high performance 
%                                   computing server
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'SaveStdOutFlag': whether to save standard outputs
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == true
%                   - 'RenewParpoolFlag': whether to refresh parallel pool
%                                           every batch
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%                   - 'MaxNumWorkers': maximum number of workers for parfor
%                                   set to 0 if parfor not to be used
%                                   set to Inf to use maximum number of workers
%                   must be a nonnegative integer or Inf
%                   default == Inf
%
% Requires:
%       cd/decide_on_parpool.m
%       cd/files2contents.m
%       cd/isaninteger.m
%       cd/ispositiveintegervector.m
%       cd/force_column_cell.m
%       cd/force_row_vector.m
%       cd/force_string_end.m
%
% Used by:    
%       cd/m3ha_network_launch.m
%       cd/m3ha_neuron_run_and_analyze.m

% File History:
% 2018-10-19 Adapted from code in ~/m3ha/run_neuron_once_4compgabab.m
% 2018-10-19 SimCommands can now be file names
% 2018-11-06 Added 'RunCommand' as an optional argument
% 2019-10-31 Added 'RenewParpoolFlag' as an optional argument
% 2019-10-31 Added 'MaxNumWorkers' as an optional argument
% 2019-10-31 Added 'SimNumbers' as an optional argument
% 

%% Hard-coded parameters
nrnivCommand = 'nrniv';
nrnguiCommand = 'nrngui';
specialFile = 'x86_64/special';
noErrorsStr = 'No_Errors!';

% The command on a high performance computing server 
%   for loading modules required for NEURON to work when the code
%   is compiled
moduleLoadCommandsHpc = sprintf([...
    '# Refresh bash shell to load environmental variables\n', ...    
    'source /etc/bashrc\n', ...
    '# Purge all loaded modules\n', ...
    'module purge\n', ...
    '# Load the newest version of the Intel compiler (16.0)\n', ...
    '#   This is necessary before loading the Open MPI library\n', ...
    'module load intel\n', ...
    '# Load the newest version of the Open MPI library (2.1.1)\n', ...
    'module load openmpi\n', ...
    '# Load the newest version of NEURON (7.4)\n', ...
    'module load neuron\n']);

%% Default values for optional arguments
simCommandsDefault = {};        % no simulation commands outside the .hoc file
                                %   by default
simNumbersDefault = [];         % run all simulations by default
openGuiFlagDefault = false;     % don't open GUI by default
runCommandDefault = '';         % set later
prefixDefault = '';             % prepend nothing to file names by default
outFolderDefault = pwd;         % use the present working directory for outputs
                                %   by default
debugFlagDefault = false;       % not in debug mode by default
onHpcFlagDefault = false;       % not on a high performance computing
                                %   server by default
saveStdOutFlagDefault = true;   % save standard outputs by default
renewParpoolFlagDefault = false;% don't refresh parallel pool every batch
maxNumWorkersDefault = Inf;     % no limits by default

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
addRequired(iP, 'hocFile', ...
    @(x) isfile(x));

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'SimCommands', simCommandsDefault, ...
    @(x) iscellstr(x) || isstring(x) || ischar(x));
addParameter(iP, 'SimNumbers', simNumbersDefault, ...
    @(x) assert(isempty(x) || ispositiveintegervector(x), ...
                'SimNumbers must be empty or a positive integer vector!'));
addParameter(iP, 'OpenGuiFlag', openGuiFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RunCommand', runCommandDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'Prefix', prefixDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'OutFolder', outFolderDefault, ...
    @(x) validateattributes(x, {'char', 'string'}, {'scalartext'}));
addParameter(iP, 'DebugFlag', debugFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'OnHpcFlag', onHpcFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'SaveStdOutFlag', saveStdOutFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'RenewParpoolFlag', renewParpoolFlagDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));
addParameter(iP, 'MaxNumWorkers', maxNumWorkersDefault, ...
    @(x) assert(isinf(x) || isscalar(x) && isaninteger(x) && x >= 0, ...
                'MaxNumWorkers must be Inf or a nonnegative integer!'));

% Read from the Input Parser
parse(iP, hocFile, varargin{:});
simCommands = iP.Results.SimCommands;
simNumbers = iP.Results.SimNumbers;
openGuiFlag = iP.Results.OpenGuiFlag;
runCommand = iP.Results.RunCommand;
prefix = iP.Results.Prefix;
outFolder = iP.Results.OutFolder;
debugFlag = iP.Results.DebugFlag;
onHpcFlag = iP.Results.OnHpcFlag;
saveStdOutFlag = iP.Results.SaveStdOutFlag;
renewParpoolFlag = iP.Results.RenewParpoolFlag;
maxNumWorkers = iP.Results.MaxNumWorkers;

% Make sure simCommands is a column cell array
simCommands = force_column_cell(simCommands);

%% Preparation
% Decide on the command for running NEURON
if isempty(runCommand)
    if openGuiFlag
        % Use nrngui
        runCommand = nrnguiCommand;
    else
        % Use nrniv
        %   Note: Dr. Carnevale said this is better than running 
        %           the special file
        runCommand = nrnivCommand;
        % runCommand = specialFile;
    end
end

% Decide on the moduleLoadCommands
if onHpcFlag && isdeployed      % if on a high performance computing server
                                %   and the code is compiled
    % Need to load modules
    moduleLoadCommands = moduleLoadCommandsHpc;
else                            % if on a local machine
    % Don't need to load modules
    moduleLoadCommands = '';
end

% Compile .mod files if specialFile does not exist and custom .mod files exist
if ~isfile(specialFile)
    % Check if .mod files exist
    modFiles = dir('*.mod');

    % Compile custom .mod files if they exist
    if ~isempty(modFiles)
        % Note: this only works for unix systems
        if isunix
            unix('nrnivmodl');
        end
    end
end

% Start a timer
if debugFlag
    timerSimulations = tic();
end

% Create a dummy command if empty
if isempty(simCommands)
    simCommands = {''};
end

% Count the number of possible simulations to run
nSims = numel(simCommands);

% Decide on the simulation numbers to run
%   Note: parfor takes a row vector
if isempty(simNumbers)
    simNumbers = 1:nSims;
else
    simNumbers = force_row_vector(simNumbers);
end

% Count actual number of simulations
nSimsActual = numel(simNumbers);

% Replace file names with file contents for simCommands
simCommandsStr = files2contents(simCommands);

% Make sure prefix ends with a '_'
prefix = force_string_end(prefix, '_', 'OnlyIfNonempty', true);

%% Run simulations
% Initialize results
simStatus = zeros(nSims, 1);    % stores simulation statuses
simStdOut = repmat({'No_Errors!'}, nSims, 1);     % stores simulation standard outputs
timeTaken = zeros(nSims, 1);    % stores simulation times

% Decide on the parallel Pool object, recreating it if necessary
poolObj = decide_on_parpool('MaxNumWorkers', maxNumWorkers, ...
                            'RenewParpoolFlag', renewParpoolFlag);

% Retrieve the number of workers
numWorkers = poolObj.NumWorkers;

% Create a counter for the number of simulations completed
nSimsCompleted = 0;                         

% Simulate while not all simulations are completed yet
while nSimsCompleted < nSimsActual
    % Determine the simulation numbers for this batch
    if renewParpoolFlag && nSimsActual - nSimsCompleted > numWorkers
        simNumbersThisBatch = simNumbers(nSimsCompleted + (1:numWorkers));
    else
        simNumbersThisBatch = simNumbers((nSimsCompleted + 1):nSimsActual);
    end

    % Count the number of simulations in this batch
    nSimsThisBatch = numel(simNumbersThisBatch);

    % Run the simulations for this batch
    parfor iSim = simNumbersThisBatch
        % Start a timer
        timerOneSim = tic();

        % Run NEURON, using a here document to append simulation commands
        %   Note: Any commands present in the .hoc file will also be run
        [simStatus(iSim), simStdOut{iSim}] = ...
            unix(sprintf(['%s\n', '%s %s - << here\n', ...
                          '%s\n', 'print "%s"\n', 'here'], ...
                         moduleLoadCommands, runCommand, ...
                         hocFile, simCommandsStr{iSim}, noErrorsStr));

        % End the timer
        timeTaken(iSim) = toc(timerOneSim);
    end

    % Recreate a parallel pool object if requested to clear leaked memory
    if renewParpoolFlag
        delete(poolObj);
        poolObj = parpool('local', numWorkers); 
    end

    % Update number of simulations completed
    nSimsCompleted = nSimsCompleted + nSimsThisBatch;
end

%% Analyze simulation standard outputs
% Check if each simulation encountered an error
hasError = ~contains(simStdOut, noErrorsStr);

% If saveStdOutFlag is true or if there is an error,
%   save simulation standard output in a text file
parfor iSim = 1:nSims
    if saveStdOutFlag || hasError(iSim)
        % Construct file path
        simOutFilePath = fullfile(outFolder, ...
            [prefix, 'simulation_', num2str(iSim), '_stdout.txt']);

        % Open the file for writing
        fid = fopen(simOutFilePath, 'w');

        % Print status and output to file
        fprintf(fid, ['Return status was: %d\n\n', ...
                      'Simulation output was:\n\n%s\n'], ...
                      simStatus(iSim), simStdOut{iSim});

        % Close the file
        fclose(fid);
    end
end

% Print to standard output
if any(hasError)
    % Look for all indices of simulations with error
    indProblematic = find(hasError);

    % Print all standard outputs for simulations with error
    for iSim = indProblematic
        fprintf(['Simulation for Sweep #%d', ...
                ' ran into errors with standard output:\n'], iSim);
        fprintf('%s\n', simStdOut{iSim});
    end
else
    % Display simulation outputs for debugging purposes 
    %   Note: these are also saved as text files if saveStdOutFlag is true
    if debugFlag
        % Print no errors
        fprintf('No Errors in NEURON!\n');

        %{
        fprintf('\n');
        for iSim = 1:nSims
            disp(simStdOut{iSim});
            fprintf('\n');
        end
        %}
    end
end

%% Display time taken
if debugFlag
    timeTakenSimulations = toc(timerSimulations);
    fprintf(['It took %3.3g seconds to run all', ...
                ' %d simulations with NEURON!!\n'], ...
            timeTakenSimulations, nSims);
    fprintf('\n');
end

%% Store outputs
output = table(simStatus, hasError, timeTaken, simStdOut);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

indProblematic = find(hasError > 0, 1);

% SOMETIMES DOESN'T PRINT STUFF FOR NO REASON!
%    error(simStdOut{indProblematic});

if saveStdOutFlag || isempty(strfind(simStdOut{iSim}, 'No_Errors!'))

hasError = isemptycell(strfind(simStdOut, 'No_Errors!'));

if ~isempty(prefix) && prefix(end) ~= '_'
    prefix = [prefix, '_'];
end

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
