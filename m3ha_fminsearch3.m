function [simplexOut, exitFlag] = m3ha_fminsearch3 (outparams)
%% Applies the Nelder-Mead simplex algorithm to optimize parameters (modified version of fminsearch for the m3ha project)
% Usage: [simplexOut, exitFlag] = m3ha_fminsearch3 (outparams)
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
%   FMINSEARCH3  Modified version of FMINSEARCH to allow use by OPTIMIZERGUI.m
%   Modified further from FMINSEARCH2 to implement parameter bounds, based
%   on "fminsearchbnd.m" (John D'Errico, 7/23/06)
%
%   FMINSEARCH Multidimensional unconstrained nonlinear minimization (Nelder-Mead).
%   X = FMINSEARCH(FUN,X0) starts at X0 and attempts to find a local minimizer 
%   X of the function FUN.  FUN is a function handle.  FUN accepts input X and 
%   returns a scalar function value F evaluated at X. X0 can be a scalar, vector 
%   or matrix.
%
%   X = FMINSEARCH(FUN,X0,OPTIONS)  minimizes with the default optimization
%   parameters replaced by values in the structure OPTIONS, created
%   with the OPTIMSET function.  See OPTIMSET for details.  FMINSEARCH uses
%   these options: Display, TolX, TolFun, MaxFunEvals, MaxIter, FunValCheck,
%   PlotFcns, and OutputFcn.
%
%   X = FMINSEARCH(PROBLEM) finds the minimum for PROBLEM. PROBLEM is a
%   structure with the function FUN in PROBLEM.objective, the start point
%   in PROBLEM.x0, the options structure in PROBLEM.options, and solver
%   name 'fminsearch' in PROBLEM.solver. The PROBLEM structure must have
%   all the fields.
%
%   [X,FVAL]= FMINSEARCH(...) returns the value of the objective function,
%   described in FUN, at X.
%
%   [X,FVAL,EXITFLAG] = FMINSEARCH(...) returns an EXITFLAG that describes 
%   the exit condition of FMINSEARCH. Possible values of EXITFLAG and the 
%   corresponding exit conditions are
%
%    1  Maximum coordinate difference between current best point and other
%       points in simplex is less than or equal to TolX, and corresponding 
%       difference in function values is less than or equal to TolFun.
%    0  Maximum number of function evaluations or iterations reached.
%   -1  Algorithm terminated by the output function.
%
%   [X,FVAL,EXITFLAG,OUTPUT] = FMINSEARCH(...) returns a structure
%   OUTPUT with the number of iterations taken in OUTPUT.iterations, the
%   number of function evaluations in OUTPUT.funcCount, the algorithm name 
%   in OUTPUT.algorithm, and the exit message in OUTPUT.message.
%
%   Examples
%     FUN can be specified using @:
%        X = fminsearch(@sin,3)
%     finds a minimum of the SIN function near 3.
%     In this case, SIN is a function that returns a scalar function value
%     SIN evaluated at X.
%
%     FUN can be an anonymous function:
%        X = fminsearch(@(x) norm(x),[1;2;3])
%     returns a point near the minimizer [0;0;0].
%
%     FUN can be a parameterized function. Use an anonymous function to
%     capture the problem-dependent parameters:
%        f = @(x,c) x(1).^2+c.*x(2).^2;  % The parameterized function.
%        c = 1.5;                        % The parameter.
%        X = fminsearch(@(x) f(x,c),[0.3;1])
%        
%   FMINSEARCH uses the Nelder-Mead simplex (direct search) method.
%
%   See also OPTIMSET, FMINBND, FUNCTION_HANDLE.
%
%   Reference: Jeffrey C. Lagarias, James A. Reeds, Margaret H. Wright,
%   Paul E. Wright, "Convergence Properties of the Nelder-Mead Simplex
%   Method in Low Dimensions", SIAM Journal of Optimization, 9(1):
%   p.112-147, 1998.
%
%   Copyright 1984-2010 The MathWorks, Inc.
%   $Revision: 1.21.4.17 $  $Date: 2010/05/13 17:39:00 $
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%
% Arguments: TODO 
%       varargin    - 'OnHpcFlag': whether on a high performance computing server
%                   must be numeric/logical 1 (true) or 0 (false)
%                   default == false
%
% Requires:
%       cd/m3ha_log_errors_params.m
%       cd/m3ha_neuron_run_and_analyze.m
%       cd/save_params.m
%       cd/set_default_flag.m
%       cd/set_fields_zero.m
%       cd/set_figure_properties.m
%       cd/restore_fields.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/m3ha_optimizer_4compgabab.m

% File History:
% 2011-01-30 - Modified by CLK
% 2016-07-21 - Added old_err
% 2017-01-17 - Now saves the error figure and params as .p & .mat files
% 2017-01-17 - Now logs everything using m3ha_log_errors_params.m
% 2017-01-17 - Changed outparams.runnum_auto to be current number 
%                   (removed "+ 1")
% 2017-01-21 - Cleaned up code
% 2017-01-21 - Changed the number of parameters to compare against to n 
%                   (from min(2, n ))
% 2017-01-21 - Changed definition of maxParamChange & tolx so that 
%                   the largest parameter doesn't dominate
% 2017-01-21 - Replace by the "reflection point" as long as 
%                   it is better than the worst point
% 2017-01-21 - Decrease relativeErrorTolerance && relativeParamTolerance from 0.05 to 0.01
% 2017-01-25 - Fixed outparams.neuronparams for the initial simplex ouput
% 2017-01-25 - Added initError, firstMaxErrorChange, firstMaxParamChange, 
%                   firsttolf, firsttolx to simplexOut
% 2017-01-25 - Fixed simplexOut.params to be re-transformed values and 
%                   change to simplexOut.paramsInUseValues
% 2017-01-25 - Absorbed outparams.neuronparams & optimizedparams into 
%                   simplexOut.neuronParamsTable
% 2017-01-25 - Fixed index of cm so that it is always within range
% 2017-05-02 - Added cprflag so that results can be saved separately
% 2017-05-02 - Changed outparams.runnum_auto -> outparams.simplexNum
% 2017-05-13 - Now gets outFolderName from outparams
% 2017-05-15 - Added simplexLogSuffix and changed disp() to fprintf()
% 2017-05-15 - Changed outparams.prefix so that '_cpr' is already incorporated
% 2017-05-15 - Now saves simplexParams as a .p file for each function evaluation
% 2017-05-16 - Changed simplexOut.maxErrorChange so that it is relative to fv(1)
% 2017-05-16 - Replaced TolX and TolFun with TolXRel and TolFunRel
% 2017-05-16 - Made optimization algorithm parameters part of simplexParams, 
%               options and simplexOut
% 2017-05-16 - Replaced optimget() with getfield() and optimset() with setfield()
% 2017-05-16 - Made 'iter' the default PrintType
% 2017-05-16 - firstMaxErrorChange, firstMaxParamChange now reflects 
%               first iteration instead of third
% 2017-05-16 - Now always updates simplex structure
% 2017-05-16 - Moved calculation of maxErrorChange and maxParamChange to 
%               within update_simplexOut()
% 2017-05-16 - Added simplexOut.totalError & simplexOut.err
% 2017-05-16 - Removed mat file saving; 
%                   simplexOut will be saved in m3ha_optimizer_4compgabab.m
% 2017-05-17 - Now uses the same prefix for all output files
% 2017-05-17 - Expanded update_errorhistory() so that it can plot lts errors
% 2017-05-18 - outparams must now have initial errors already set
% 2017-05-18 - Changed sin & arcsin to tan & arctan
% 2017-05-18 - Changed each point in initial simplex to reflect an 
%                   absolute change of UsualDelta * pi in each direction
% 2017-05-18 - Removed ZeroTermDelta and changed UsualDelta from 0.5 to 0.75
% 2017-05-18 - Now transforms log-scaled parameters to log-space first
% 2017-05-18 - Normalize parameter change by pi; no longer shift by 2*pi
% 2017-05-19 - Now plots sweeps for the first and last function evaluations 
%                   when outparams.plotsweepsflag is true 
% 2017-05-20 - Changed prefix to outparams.prefix
% 2017-05-22 - Changed line width and indentation
% 2017-05-23 - Added otherwise to all switch statements
% 2018-01-24 - Added isdeployed
% 2018-03-02 - Added onHpcFlag as an optional argument
% 2019-11-25 - Improved titel for simplex history figure
% TODO: Change all contents of outparams used to varargin
% TODO: Make neuronParamsTable a required argument
%

%% Hard-coded parameters
% Simplex log file
simplexLogSuffix = 'simplexlog.txt';

% For plotting
figNumberSimplexHistory = 301;

%% Default options
defaultopt = struct('PrintType', 'iter', ...        % what to put in log file
    'MaxIter', '200*numberOfVariables', ...         % maximum number of iterations
    'MaxFunEvals', '200*numberOfVariables', ...     % maximum number of function evaluations
    'TolFunRel', 0.01, ...  % relative error tolerance (w.r.t. smallest error)
    'TolXRel', 0.01, ...    % relative parameter change tolerance (w.r.t. best set of parameters)
    'UsualDelta', 0.75, ... % increment for transformed parameters multiple of pi (must be in (0 1))
    'Rho', 1, ...           % used in the computation of the "reflection point" and others
    'Chi', 2, ...           % used in the computation of the "expansion point"
    'Psi', 0.5, ...         % used in the computation of the "contraction points"
    'Sigma', 0.5);          % used in the performance of "shrink"

%% Default values for optional arguments
%{
onHpcFlagDefault = false;       % whether on a high performance computing
                                %   server by default
%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Check number of required arguments
if nargin < 1
    error(create_error_for_nargin(mfilename));
end

% Set up Input Parser Scheme
iP = inputParser;         
iP.FunctionName = mfilename;

%{
% Add required inputs to the Input Parser
%% TODO

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'OnHpcFlag', onHpcFlagDefault, ...   
    @(x) validateattributes(x, {'logical', 'numeric'}, {'binary'}));

% Read from the Input Parser
parse(iP, varargin{:});
onHpcFlag = iP.Results.OnHpcFlag;
%}

%% Extract from outparams
% TODO: Change to use varargin
if isfield(outparams, 'initialSimplexMode')
    initialSimplexMode = outparams.initialSimplexMode;
else
    initialSimplexMode = 0;
end
multipleRunsFlag = outparams.multipleRunsFlag;
ltsFeatureWeights = outparams.ltsFeatureWeights;
showSimplexHistoryFlag = outparams.showSimplexHistoryFlag;
autoLoggerFlag = outparams.autoLoggerFlag;

outFolderName = outparams.outFolder;
prefix = outparams.prefix;
simplexNum = outparams.simplexNum;
simMode = outparams.simMode;
simplexIterCount = outparams.simplexIterCount;
simplexParams = outparams.simplexParams;
simplexParamNames = outparams.simplexParamNames;
neuronParamsTableInit = outparams.neuronParamsTable;

% Temporary measure
% TODO: Change all instances of outFolderName to outFolder
outparams.outFolder = outFolderName;

%% Preparation
% Determine whether LTS errors are computed
computeLtsError = set_default_flag([], strcmp(simMode, 'active') && ...
                                ~all(ltsFeatureWeights == 0));

% Extract from neuronParamsTableInit
pNames = neuronParamsTableInit.Properties.RowNames;
initValues = neuronParamsTableInit.Value;
pIsLog = neuronParamsTableInit.IsLog;
pInUse = neuronParamsTableInit.InUse;
pLowerBound = neuronParamsTableInit.LowerBound;
pUpperBound = neuronParamsTableInit.UpperBound;

% Count the number of parameters
nParams = height(neuronParamsTableInit);

% Create new prefix
prefix = [prefix, '_simplexrun_', num2str(simplexNum)];

% Create figure to track simplex performance
fprintf('Plotting figure to track simplex performance for %s ...\n\n', prefix);
if showSimplexHistoryFlag
    visibleStatus = 'on';
else
    visibleStatus = 'off';
end
simplexfigure = set_figure_properties('ClearFigure', true, ...
                'Visible', visibleStatus, ...
                'FigNumber', figNumberSimplexHistory, ...
                'Name', 'Simplex performance tracker');

% Set color map with line colors
cm = colormap(lines);

% Open log file
fidLog = fopen(fullfile(outFolderName, [prefix, '_', simplexLogSuffix]), 'w');

%% Deal with simplex parameters
% Initialize options
options = defaultopt;

% Update options with simplexParamNames & simplexParams
for k = 1:numel(simplexParams)
    options = setfield(options, simplexParamNames{k}, simplexParams(k));
end

% Update simplex parameters from options
printType = getfield(options, 'PrintType');
maxFunctionEvaluations = getfield(options, 'MaxFunEvals');
maxIterations = getfield(options, 'MaxIter');
relativeErrorTolerance = getfield(options, 'TolFunRel');
relativeParamTolerance = getfield(options, 'TolXRel');
usualDelta = getfield(options, 'UsualDelta');
rho = getfield(options, 'Rho');
chi = getfield(options, 'Chi');
psi = getfield(options, 'Psi');
sigma = getfield(options, 'Sigma');

% Store simplex parameters in simplexOut
simplexOut.printType = printType;
simplexOut.maxFunctionEvaluations = maxFunctionEvaluations;
simplexOut.maxIterations = maxIterations;
simplexOut.relativeErrorTolerance = relativeErrorTolerance;
simplexOut.relativeParamTolerance = relativeParamTolerance;
simplexOut.usualDelta = usualDelta;
simplexOut.rho = rho;
simplexOut.chi = chi;
simplexOut.psi = psi;
simplexOut.sigma = sigma;

% Initialize the NEURON parameters table in simplexOut
simplexOut.neuronParamsTable = outparams.neuronParamsTable;

%% Determine parameter bounds
% Type of parameter bounding (cf. D'Errico)
% 0 --> unconstrained variable
% 1 --> lower bound only
% 2 --> upper bound only
% 3 --> dual finite bounds
% 4 --> fixed variable
% Current only supports case 3 or 4
bounds.BoundClass = zeros(nParams, 1);
bounds.LB = zeros(nParams, 1);
bounds.UB = zeros(nParams, 1);
for k = 1:nParams
    if pInUse(k)    % parameter should be used for fitting
        isfinitenum = @(x) isnumeric(x) && isfinite(x);
        % Check if parameter has dual finite bounds, otherwise exit fitting
        if isfinitenum(pLowerBound(k)) ...
            && isfinitenum(pUpperBound(k)) ...
            bounds.BoundClass(k) = 3;       % 'dual finite bounds'
            bounds.LB(k) = pLowerBound(k);  % lower bound
            bounds.UB(k) = pUpperBound(k);  % upper bound
        else
            fprintf('Parameter #%d must have dual finite bounds!', k);
            exitFlag = -1;
            return
        end
    else            % parameter should be fixed
        bounds.BoundClass(k) = 4;           % 'fixed'
        bounds.LB(k) = initValues(k);       % use initial value
        bounds.UB(k) = bounds.LB(k);        % use initial value
    end
end

%% Determine number of parameters in use
n = sum(bounds.BoundClass == 3);        % number of parameters in use
ncp = n;                                % number of points to compare against
two2np1 = 2:n+1;
one2n = 1:n;

% Compute maxFunctionEvaluations & maxIterations in case the defaults were gathered from calling optimget
if ischar(maxFunctionEvaluations)
    if isequal(lower(maxFunctionEvaluations), '200*numberofvariables')
        maxFunctionEvaluations = 200 * n;
    else
        error('MATLAB:fminsearch:OptMaxFunEvalsNotInteger', ...
            'Option ''MaxFunEvals'' must be an integer value if not the default.')
    end
end
if ischar(maxIterations)
    if isequal(lower(maxIterations), '200*numberofvariables')
        maxIterations = 200 * n;
    else
        error('MATLAB:fminsearch:OptMaxIterNotInteger', ...
            'Option ''MaxIter'' must be an integer value if not the default.')
    end
end

%% Initialize parameters
% For parameters used for fitting, 
% transform initial values into the range [3*pi/2 5*pi/2] nonlinearly using arctan
% Also check for out-of-range initial values (heavily borrowed from D'Errico's code)
p0Log = zeros(n, 1);        % stores parameter values transformed to log-values
p0N121 = zeros(n, 1);       % stores parameter values transformed to [-1 1]
p0N121Corr = zeros(n, 1);   % p0N121 corrected
p0 = zeros(n, 1);       % stores parameter values transformed to [3*pi/2 5*pi/2]
ct = 0;                 % counts unfixed parameters
for i = 1:nParams
    switch bounds.BoundClass(i)
    case 3        % 'dual finite bounds' parameter
        % Increment count
        ct = ct + 1;

        % Record parameter name
        simplexOut.paramsInUseNames{ct} = pNames{i};

        % First transform [LB UB] or [log(LB) log(UB)] linearly to [-1 1]
        if pIsLog(i)
            p0N121(ct) = 2 * (log(initValues(i)) - log(bounds.LB(i))) / ...
                        (log(bounds.UB(i)) - log(bounds.LB(i))) - 1;
        else
            p0N121(ct) = 2 * (initValues(i) - bounds.LB(i)) / ...
                        (bounds.UB(i) - bounds.LB(i)) - 1;
        end

        % Collapse to bounds if initial value out of bounds
        p0N121Corr(ct) = max(-1, min(1, p0N121(ct)));

        % Then transform to [-pi/2 pi/2] nonlinearly with arctan
        p0(ct) = atan(p0N121Corr(ct));
    case 4        % 'fixed' parameter
        % Drop before fminsearch sees it
        % count is not incremented for this parameter
    otherwise
        error('This BoundClass has not been implemented!\n');
    end 
end

% If all variables were fixed, there is no need to perform simplex
if ct == 0        % all variables were fixed
    simplexOut.lastHow = '';
    simplexOut.ctEvals = 0;
    simplexOut.ctIterations = 0;
    simplexOut.paramsInUseValues = [];
    simplexOut.totalError = NaN;
    fprintf(fidLog, 'All variables were fixed; end without simplex run\n\n');
    exitFlag = -1;
    return
end

%% Initialize simplex: a convex region in the n-dimensional space with n+1 vertices
ctIterations = 0;             % count of simplex iterations
ctEvals = 0;            % count of error evaluations
how = 'initial';        % stores the method of current iteration
v = zeros(n, n + 1);    % the simplex that stores parameters from all iterations
fv = zeros(1, n + 1);   % the total errors corresponding to each set of parameters of the simplex
errv = cell(1, n + 1);  % the error structures corresponding to each set of parameters of the simplex

% Compute initial total error and plot sweeps if user requests so
[err0, ctEvals] = eval_with_plots('_bef', p0, bounds, pIsLog, ctEvals, outparams);

% Update simplex performance figure
update_errorhistory(simplexfigure, ctIterations, err0, cm, ...
    multipleRunsFlag, computeLtsError, simplexIterCount);

% Place initial params in the simplex! (credit L.Pfeffer at Stanford)
v(:, 1) = p0;
fv(1) = err0.totalError;
errv{1} = err0;

% Initialize prnt
switch printType
case {'notify', 'notify-detailed'}
    prnt = 1;
case {'none', 'off'}
    prnt = 0;
case {'iter', 'iter-detailed'}
    prnt = 3;
case {'final', 'final-detailed'}
    prnt = 2;
case 'simplex'
    prnt = 4;
otherwise
    prnt = 3;
    fprintf('Warning: PrintType unrecognised, use default ''iter''!\n\n');
end

% Print out initial f(x) as 0th iteration
switch prnt
case 3
    fprintf(fidLog, ' \n');
    fprintf(fidLog, ' Iteration   Func-count     min f(x)         Procedure\n');
    fprintf(fidLog, sprintf(' %5.0f        %5.0f     %12.6g         %s\n', ctIterations, ctEvals, fv(1), how));
case 4
    clc
    formatsave = get(0, {'format', 'formatspacing'});    % save default format of objects
    format compact
    format short e
    fprintf(fidLog, ' \n');
    fprintf(fidLog, '%s\n', how);
    v
    fv
    ctEvals
otherwise
end

% Log initial errors and params
if autoLoggerFlag
    % Name of csv file to log errors and params
    logFileName = [prefix, '_errors_and_params_log.csv'];

    % Log initial errors and params in .csv file
    simplexOut.initError = fv(1);
    simplexOut = update_simplexOut(simplexOut, how, ctIterations, ctEvals, fv, errv, v, bounds, ncp, pIsLog);
    m3ha_log_errors_params(logFileName, outparams, errv{1}, simplexOut);
end

%% First iteration of simplex: set up simplex near the initial set of parameters
% Following improvement suggested by L. Pfeffer at Stanford
for j = 1:n    % for each parameter
    % Construct a set of parameters with this parameter altered only
    y = p0;            % start with initial set of parameters; all values are within [-pi/2 pi/2]

    % Decide on the change for the new vertex
    if initialSimplexMode == 0
        % Increment parameter j by usualDelta * pi
        y(j) = y(j) + usualDelta * pi;
    elseif initialSimplexMode == 1
        % Increment all parameters but j by usualDelta * pi
        for iParam = 1:n
            if iParam ~= j
                y(iParam) = y(iParam) + usualDelta * pi;
            end
        end
    else 
        error('Initial Simplex Mode %d has not been implemented!', ...
                initialSimplexMode);
    end

    % Place new set of parameters in simplex and evaluate its associated errors
    v(:, j+1) = y;
    [errv{j+1}, ctEvals] = eval_with_bounds(y, bounds, pIsLog, ctEvals, outparams);
    fv(j+1) = errv{j+1}.totalError;
end

% Sort simplex so that v(:, 1) has the lowest total error
[fv, I] = sort(fv);    % sort fv in ascending order
v = v(:, I);        % rearrange v in same order
errv = errv(I);        % rearrange errv in same order

% Update counts and record what was done
how = 'initial simplex';
ctIterations = ctIterations + 1;

% Update simplexOut structure and compute maximum error change, maximum parameter change
simplexOut = update_simplexOut(simplexOut, how, ctIterations, ctEvals, fv, errv, v, bounds, ncp, pIsLog); 

% Log results and update simplex performance figure
switch prnt
case 3
    fprintf(fidLog, sprintf(' %5.0f        %5.0f     %12.6g         %s\n', ctIterations, ctEvals, fv(1), how));
case 4
    fprintf(fidLog, ' \n');
    fprintf(fidLog, '%s\n', how);
    v
    fv
    ctEvals
otherwise
end

% Update simplex performance figure
update_errorhistory(simplexfigure, ctIterations, errv{1}, cm, ...
    multipleRunsFlag, computeLtsError, simplexIterCount);

% Log errors and params
if autoLoggerFlag
    m3ha_log_errors_params(logFileName, outparams, errv{1}, simplexOut);
end

%% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% next best ncp other points in the simplex is less than or equal to relativeParamTolerance. 
% Specifically, until max(|v2(i)-v1(i)|/|v1(i)|, |v2(i)-v1(i)|/|v1(i)|, 
%                ..., |v(ncp+1)(i)-v1(i)|/|v1(i)|) <= TolXRel for all i in 1:n
% where v1 holds the vertex with the current lowest function value; AND (Cannot use OR instead of AND.)
% (b) the corresponding relative differences in function values is less than or equal to TolFunRel; OR
% (c) the maximum number of iterations is exceeded; OR
% (d) the maximum number of function evaluations is exceeded

simplexOut.algorithm = 'Nelder-Mead simplex direct search';
while ctIterations < maxIterations && ctEvals < maxFunctionEvaluations ...
    && (simplexOut.maxErrorChange > relativeErrorTolerance || simplexOut.maxParamChange > relativeParamTolerance)

    % Find the worst point (pworst) and compute the average of the better n points (pbar)
    pworst = v(:, end);            % set of params that give the largest error in the last iteration
    pbar = sum(v(:, one2n), 2) / n;        % average of the n better sets of parameters

    % Compute the "reflection point" (pr): 
    %    the point rho*||pbar - pworst|| away from pbar
    %    in the opposite direction of pworst
    pr = (1 + rho) * pbar - rho * pworst;    

    % Compute the error (fpr) associated with the "reflection point"
    [errpr, ctEvals] = eval_with_bounds(pr, bounds, pIsLog, ctEvals, outparams);
    fpr = errpr.totalError;

    if fpr < fv(:, 1)    % if the "reflection point" is better than the previous best point
        % Compute the "expansion point" (pe):
        %    the point chi*rho*||pbar - pworst|| away from pbar
        %    in the opposite direction of pworst
        pe = (1 + chi * rho) * pbar - chi * rho * pworst;
        [errpe, ctEvals] = eval_with_bounds(pe, bounds, pIsLog, ctEvals, outparams);
        fpe = errpe.totalError;

        % Replace the worst point with the better of the "expansion point" and the "reflection point"
        if fpe < fpr    
            how = 'expand: < best point';
            v(:, end) = pe;
            fv(end) = fpe;
            errv{end} = errpe;
        else
            how = 'reflect: < best point';
            v(:, end) = pr;
            fv(end) = fpr;
            errv{end} = errpr;
        end
    else            % if the "reflection point" is not better than the previous best point
        if fpr < fv(end-1)    % if the "reflection point" is better than the second worst point
            % Replace the worst point with the "reflection point"
            how = 'reflect: < second worst point';
            v(:, end) = pr;
            fv(end) = fpr;
            errv{end} = errpr;
        else        % if the "reflection point" is not better than the second worst point
            % Perform a "contraction" to see if any better
            if fpr < fv(end)    % if the "reflection point" is better than the worst point
                % Compute the "outside contraction point" (pc):
                %    the point psi*rho*||pbar - pworst|| away from pbar
                %    in the opposite direction of pworst
                pc = (1 + psi * rho) * pbar - psi * rho * pworst;
                [errpc, ctEvals] = ...
                    eval_with_bounds(pc, bounds, pIsLog, ctEvals, outparams);
                fpc = errpc.totalError;

                % Replace the worst point with the better of 
                %    the "outside contraction point" and the "reflection point"
                if fpc <= fpr
                    how = 'contract outside: < worst point';
                    v(:, end) = pc;
                    fv(end) = fpc;
                    errv{end} = errpc;
                else
                    how = 'reflect: < second worst point';
                    v(:, end) = pr;
                    fv(end) = fpr;
                    errv{end} = errpr;
                end
            else            % if the "reflection point" is not better than the worst point
                % Compute the "inside contraction point" (pcc):
                %    the point psi*||pbar - pworst|| away from pbar
                %    in the SAME direction as pworst
                pcc = (1 - psi) * pbar + psi * pworst;
                [errpcc, ctEvals] = ...
                    eval_with_bounds(pcc, bounds, pIsLog, ctEvals, outparams);
                fpcc = errpcc.totalError;

                % Replace the worst point with the "inside contraction point" 
                %     if it's better than the worst point
                %     Otherwise, no direction of replacement is better; perform a "shrink"
                if fpcc < fv(end)
                    how = 'contract inside: < worst point';
                    v(:, end) = pcc;
                    fv(end) = fpcc;
                    errv{end} = errpcc;
                else
                    how = 'shrink';
                    for j = two2np1
                        % Replace all points p other than the best point
                        %    with the point sigma*||p - pbest|| away from pbest
                        %    in the direction of p
                        v(:, j) = v(:, 1) + sigma * (v(:, j) - v(:, 1));
                        psh = v(:, j); 
                        [errpsh, ctEvals] = ...
                            eval_with_bounds(psh, bounds, pIsLog, ctEvals, outparams);
                        fv(j) = errpsh.totalError;
                        errv{j} = errpsh;
                    end

                end
            end
            
        end
    end

    % Find best set of parameters in this iteration and update counts
    [fv, I] = sort(fv);         % sort errors in ascending order
    v = v(:, I);                % rearrange the parameters in same order
    errv = errv(I);             % rearrange errv in same order
    ctIterations = ctIterations + 1;  % increase iteration count


    % Update simplexOut structure and compute maximum error change, maximum parameter change
    simplexOut = update_simplexOut(simplexOut, how, ctIterations, ctEvals, fv, errv, v, bounds, ncp, pIsLog);

    % Log results of current iteration and update simplex performance figure
    switch prnt
    case 3
        fprintf(fidLog, sprintf(' %5.0f        %5.0f     %12.6g         %s\n', ctIterations, ctEvals, fv(1), how));
    case 4
        fprintf(fidLog, ' \n');
        fprintf(fidLog, '%s\n', how);
        v
        fv
        ctEvals
    otherwise
    end

    % Update simplex performance figure
    update_errorhistory(simplexfigure, ctIterations, errv{1}, cm, ...
                multipleRunsFlag, computeLtsError, simplexIterCount);

    % Log errors and params
    if autoLoggerFlag
        m3ha_log_errors_params(logFileName, outparams, errv{1}, simplexOut);
    end
end

%% Save optimization performance figure and close it
suptitle(['Simplex run #', num2str(simplexNum), ...
        ' for ', replace(prefix, '_', '\_')]);
saveas(simplexfigure, fullfile(outFolderName, [prefix, '_history.png']));
close(simplexfigure);

%% Save best parameters for this simplex run (always do this)
sheetPath = fullfile(outFolderName, ...
                    [prefix, '_best_params.csv']);
save_params(simplexOut.neuronParamsTable, 'FileName', sheetPath);

%% Generate exitFlag
if ctEvals >= maxFunctionEvaluations
    msg = sprintf(['Exiting: Maximum number of function evaluations has been exceeded\n' ...
            '         - increase MaxFunEvals option.\n' ...
            '         Current function value: %f \n'], fv(1));
    if prnt > 0
        fprintf(fidLog, ' \n');
        fprintf(fidLog, '%s\n', msg);
    end
    exitFlag = 0;
elseif ctIterations >= maxIterations
    msg = sprintf(['Exiting: Maximum number of iterations has been exceeded\n' ... 
            '         - increase MaxIter option.\n' ...
            '         Current function value: %f \n'], fv(1));
    if prnt > 0
        fprintf(fidLog, ' \n');
        fprintf(fidLog, '%s\n', msg);
    end
    exitFlag = 0;
else
    msg = sprintf(['Optimization terminated:\n', ...
            ' the current x satisfies the termination criteria using OPTIONS.TolXRel of %e \n' ...
            ' and F(X) satisfies the convergence criteria using OPTIONS.TolFunRel of %e \n'], ...
            relativeParamTolerance, relativeErrorTolerance);
    if prnt > 1
        fprintf(fidLog, ' \n');
        fprintf(fidLog, '%s\n', msg);
    end
    exitFlag = 1;
end
simplexOut.message = msg;

%% Reset default format of objects
if prnt == 4
    set(0, {'format', 'formatspacing'}, formatsave);
end

%% Close log file and figure
fclose(fidLog);

%% Evaluate again and plot if requested
eval_with_plots('_aft', v(:, 1), bounds, pIsLog, ctEvals, outparams)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function simplexOut = update_simplexOut(simplexOut, how, ctIterations, ctEvals, fv, errv, v, bounds, ncp, pIsLog)
%% Update simplexOut structure and compute maximum error change, maximum parameter change

% Update simplexOut structure
simplexOut.lastHow = how;
simplexOut.ctIterations = ctIterations;
simplexOut.ctEvals = ctEvals;
simplexOut.totalError = fv(:, 1);
simplexOut.err = errv{1};

% Compute maximum error change, maximum parameter change
if ctIterations == 0        % simplex has not been set up
    simplexOut.maxErrorChange = NaN;
    simplexOut.maxParamChange = NaN;
else
    simplexOut.maxErrorChange = max(abs(fv(1) - fv(2:ncp+1))) / fv(1);
    simplexOut.maxParamChange = max(max(abs(v(:, 2:ncp+1) - v(:, ones(1, ncp))) ./ pi));
end

% Save first maximum error change, maximum parameter change, used in m3ha_optimizer_4compgabab.m
if ctIterations == 1
    simplexOut.firstMaxErrorChange = simplexOut.maxErrorChange;
    simplexOut.firstMaxParamChange = simplexOut.maxParamChange;
end

% Read from simplexOut
neuronParamsTable = simplexOut.neuronParamsTable;

% Undo the variable transformations
[newParamValues, paramsInUseValues] = xtransform(v(:, 1), bounds, pIsLog);

% Place new values in the new neuronParamsTable
neuronParamsTable{:, 'Value'} = newParamValues;

% Update simplexOut
simplexOut.neuronParamsTable = neuronParamsTable;
simplexOut.paramsInUseValues = paramsInUseValues;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_errorhistory(fig, ctIterations, err, cm, ...
                    multipleRunsFlag, computeLtsError, simplexIterCount)
% Update simplex performance figure
% TODO: Merge with code in m3ha_optimizer_4compgabab.m

% Choose the iteration number
if multipleRunsFlag
    % Total number of times NEURON was run
    iter = simplexIterCount + ctIterations;
else
    % Iteration count of current simplex
    iter = ctIterations;
end

% Find a color for this iteration
color = cm(mod(iter, size(cm, 1)) + 1, :);        

% Choose the marker style
if ctIterations == 0
    % The initial error is a filled circle
    markerstyle = 'o';
else
    % The rest are dots
    markerstyle = '.';
end

set(0, 'CurrentFigure', fig);
if computeLtsError   % if lts error is computed
    % Plot the total error
    subplot(3, 2, 1);
    update_subplot(iter, err.totalError, [], 'total error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

    % Plot the average sweep error
    subplot(3, 2, 2);
    update_subplot(iter, err.avgSwpError, [], 'sweep error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

    % Plot the average LTS error
    subplot(3, 2, 3);
    update_subplot(iter, err.ltsMatchError, [], 'match error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

    % Plot the average LTS amp error
    subplot(3, 2, 4);
    update_subplot(iter, err.avgLtsAmpError, [], 'amp error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

    % Plot the average LTS time error
    subplot(3, 2, 5);
    update_subplot(iter, err.avgLtsDelayError, [], 'time error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

    % Plot the average LTS slope error
    subplot(3, 2, 6);
    update_subplot(iter, err.avgLtsSlopeError, [], 'slope error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);

else
    % If no lts error is computed, just plot the total error
    update_subplot(iter, err.totalError, 'Number of iterations', 'total error', ...
            markerstyle, color, ctIterations, multipleRunsFlag);
end

% Update plot as simplex proceeds
drawnow

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function update_subplot(x, y, xLabel, yLabel, markerstyle, color, ctIterations, multipleRunsFlag)
% Update subplot

% Set an appropriate upper y limit value
if y > 0
    yLimitsMax = y * 1.1;
else
    yLimitsMax = 1;
end

% Plot new marker
plot(x, y, 'Marker', markerstyle, 'Color', color, 'MarkerFaceColor', 'auto');

% Adjust y limits
if ctIterations == 0
    hold on;
    ylim([0, yLimitsMax]);
    if ~isempty(xLabel)
        xlabel(xLabel); 
    end
    if ~isempty(yLabel)
        ylabel(yLabel);
    end
else
    ylimits = get(gca, 'YLim');
    if ylimits(2) < y        % rescale axes if error is greater than axis limit
        ylim([0, yLimitsMax]);
    end
end
if multipleRunsFlag            % if simplex method is run multiple times
    if ylimits(2) > y * 10        % rescale axes if error is much smaller than initial error
        ylim([0, yLimitsMax]);
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err, ctEvals] = eval_with_plots(prefixMod, p, bounds, pIsLog, ctEvals, fixedParams)

% Restore plotIndividualFlag and plotOverlappedFlag
fixedParams = ...
    restore_fields(fixedParams, 'plotIndividualFlag', 'plotOverlappedFlag');

% Change prefix
fixedParams.prefix = [fixedParams.prefix, prefixMod];

% Run and plot sweeps
[err, ctEvals] = eval_with_bounds(p, bounds, pIsLog, ctEvals, fixedParams);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [err, ctEvals] = eval_with_bounds(p, bounds, pIsLog, ...
                                            funcEvalsOld, fixedParams)
%% Update neuronParamsTable and find associated errors

% Initialize a table from the previously stored neuronParamsTable
neuronParamsTable = fixedParams.neuronParamsTable;

% Transform variables to get new values
newParamValues = xtransform(p, bounds, pIsLog);

% Place new values in the new neuronParamsTable
neuronParamsTable{:, 'Value'} = newParamValues;

% Evaluate error with the new parameters table
err = m3ha_neuron_run_and_analyze(neuronParamsTable, fixedParams);

% Update evaluation count
ctEvals = funcEvalsOld + 1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [x, xpart] = xtransform(p, bounds, pIsLog)
%% Converts unconstrained parameters (an angle) into their original domains and replace fixed parameters
%  Note that parameters in p are unconstrained because sine is periodical

% Extract bounds
UB = bounds.UB;
LB = bounds.LB;

% Initialize variables
numinuse = numel(p);        % number of parameters in use
nParams = numel(LB);      % number of original parameters
p_01 = zeros(numinuse, 1);    % stores parameters transformed to [0 1]
xpart = zeros(numinuse, 1);    % stores transformed parameters
x = zeros(nParams, 1);    % stores transformed parameters and original parameters

ct = 0;                        % counts unfixed parameters that were not dropped from optimization
for i = 1:numel(bounds.LB)
    switch bounds.BoundClass(i)
    case 3        % 'dual finite bounds' parameter
        % Increment count
        ct = ct + 1;

        % First transform back to [0 1]
        p_01(ct) = (tan(p(ct)) + 1) / 2;

        % Then transform linearly to original interval [LB UB] or [log(LB) log(UB)]
        if pIsLog(i)
            x(i) = exp(p_01(ct) * (log(UB(i)) - log(LB(i))) + log(LB(i)));
        else
            x(i) = p_01(ct) * (UB(i) - LB(i)) + LB(i);
        end

        % Just in case of any floating point problems
        x(i) = max(LB(i), min(UB(i), x(i)));
        xpart(ct) = x(i);

    case 4        % 'fixed' parameter
        % Since bounds were set to initial value, set it at either bound
        x(i) = LB(i);
    otherwise
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
