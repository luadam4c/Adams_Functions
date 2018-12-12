function [simplexOut, exitFlag] = m3ha_fminsearch3(outparams)
%% Applies the Nelder-Mead simplex algorithm to optimize parameters (modified version of fminsearch for the m3ha project)
% Usage: [simplexOut, exitFlag] = m3ha_fminsearch3(outparams)
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
%       cd/m3ha_run_neuron_once.m
%       cd/save_params.m
%       cd/set_fields_zero.m
%       cd/restore_fields.m
%
% Used by:
%       ~/m3ha/optimizer4gabab/optimizer_4compgabab.m

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
% 2017-05-15 - Added simplexLogFile and changed disp() to fprintf()
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
%                   simplexOut will be saved in optimizer_4compgabab.m
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
% TODO: Change all contents of outparams used to varargin
% TODO: Make neuronParamsTable a required argument
%

% Simplex log file
simplexLogFile = 'simplexlog.txt';

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
    error(['Not enough input arguments, ', ...
            'type ''help %s'' for usage'], mfilename);
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
multipleRunsFlag = outparams.multipleRunsFlag;
ltsErrorFlag = outparams.ltsErrorFlag;
showSimplexHistoryFlag = outparams.showSimplexHistoryFlag;
autoLoggerFlag = outparams.autoLoggerFlag;

outFolderName = outparams.outFolderName;
prefix = outparams.prefix;
simplexNum = outparams.simplexNum;
simMode = outparams.simMode;
simplexIterCount = outparams.simplexIterCount;
simplexParams = outparams.simplexParams;
simplexParamNames = outparams.simplexParamNames;
neuronParamsTableInit = outparams.neuronParamsTable;

%% Preparation
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
    simplexfigure = figure(301);
else
    simplexfigure = figure('Visible', 'off');
end
set(simplexfigure, 'Name', 'Simplex performance tracker');
clf(simplexfigure);
cm = colormap(lines);   % set color map with line colors

% Open log file
fidLog = fopen(fullfile(outFolderName, [prefix, '_', simplexLogFile]), 'w');

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
            bounds.BoundClass(k) = 3;                   % 'dual finite bounds'
            bounds.LB(k) = pLowerBound(k);   % lower bound
            bounds.UB(k) = pUpperBound(k);   % upper bound
        else
            fprintf('Parameter #%d must have dual finite bounds!', k);
            exitFlag = -1;
            return
        end
    else                                % parameter should be fixed
        bounds.BoundClass(k) = 4;                       % 'fixed'
        bounds.LB(k) = initValues(k);           % use initial value
        bounds.UB(k) = bounds.LB(k);                        % use initial value
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
p0_log = zeros(n, 1);        % stores parameter values transformed to log-values
p0_n121 = zeros(n, 1);        % stores parameter values transformed to [-1 1]
p0_n121corr = zeros(n, 1);    % p0_n121 corrected
p0 = zeros(n, 1);        % stores parameter values transformed to [3*pi/2 5*pi/2]
ct = 0;                % counts unfixed parameters
for i = 1:nParams
    switch bounds.BoundClass(i)
    case 3        % 'dual finite bounds' parameter
        % Increment count
        ct = ct + 1;

        % Record parameter name
        simplexOut.paramsInUseNames{ct} = pNames{i};

        % First transform [LB UB] or [log(LB) log(UB)] linearly to [-1 1]
        if pIsLog(i)
            p0_n121(ct) = 2 * (log(initValues(i)) - log(bounds.LB(i))) / ...
                        (log(bounds.UB(i)) - log(bounds.LB(i))) - 1;
        else
            p0_n121(ct) = 2 * (initValues(i) - bounds.LB(i)) / ...
                        (bounds.UB(i) - bounds.LB(i)) - 1;
        end

        % Collapse to bounds if initial value out of bounds
        p0_n121corr(ct) = max(-1, min(1, p0_n121(ct)));

        % Then transform to [-pi/2 pi/2] nonlinearly with arctan
        p0(ct) = atan(p0_n121corr(ct));
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
    multipleRunsFlag, ltsErrorFlag, simplexIterCount, simMode);

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
for j = 1:n    % for each parame  ter
    % Construct a set of parameters with this parameter altered only
    y = p0;            % start with initial set of parameters; all values are within [-pi/2 pi/2]
    y(j) = y(j) + usualDelta * pi;    % increment parameter j by usualDelta * pi

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
    multipleRunsFlag, ltsErrorFlag, simplexIterCount, simMode);

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
                multipleRunsFlag, ltsErrorFlag, simplexIterCount, simMode);

    % Log errors and params
    if autoLoggerFlag
        m3ha_log_errors_params(logFileName, outparams, errv{1}, simplexOut);
    end
end

%% Save optimization performance figure and close it
suptitle(['Simplex run #', num2str(simplexNum)])
saveas(simplexfigure, fullfile(outFolderName, [prefix, '_history.png']));
close(simplexfigure);

%% Save best parameters for this simplex run (always do this)
sheetPath = fullfile(outFolderName, ...
                    [prefix, '_best_params.xlsx']);
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

% Save first maximum error change, maximum parameter change, used in optimizer_4compgabab.m
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
                    multipleRunsFlag, ltsErrorFlag, simplexIterCount, simMode)
% Update simplex performance figure

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
if strcmpi(simMode, 'active') && ltsErrorFlag   % if lts error is computed

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
    update_subplot(iter, err.avgLtsError, [], 'LTS error', ...
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

% Plot new marker
plot(x, y, 'Marker', markerstyle, 'Color', color, 'MarkerFaceColor', 'auto');

% Adjust y limits
if ctIterations == 0
    hold on;
    initymax = y * 1.1;        % initial ymax
    ylim([0 initymax]);
    if ~isempty(xLabel)
        xlabel(xLabel); 
    end
    if ~isempty(yLabel)
        ylabel(yLabel);
    end
else
    ylimits = get(gca, 'YLim');
    if ylimits(2) < y        % rescale axes if error is greater than axis limit
        ylim([0, y * 1.1]);
    end
end
if multipleRunsFlag            % if simplex method is run multiple times
    if ylimits(2) > y * 10        % rescale axes if error is much smaller than initial error
        ylim([0, y * 1.1]);
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

function [err, ctEvals] = eval_with_bounds(p, bounds, pIsLog, funcEvalsOld, fixedParams)
%% Update neuronParamsTable and find associated errors

% Initialize a table from the previously stored neuronParamsTable
neuronParamsTable = fixedParams.neuronParamsTable;

% Transform variables to get new values
newParamValues = xtransform(p, bounds, pIsLog);

% Place new values in the new neuronParamsTable
neuronParamsTable{:, 'Value'} = newParamValues;

% Evaluate error with the new parameters table
err = m3ha_run_neuron_once(neuronParamsTable, fixedParams);

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
%% OLD CODE

[outparams.neuronparams, ~ ] = xtransform(v(:, 1), bounds);    % undo the variable transformations
optimizedparams = outparams.neuronparams;
outparams.neuronparams = reshape(xtransform(v(:, 1), bounds), [], 1)';    % undo the variable transformations

            if strcmp(how, 'shrink')
            end

                % Replace the worst point with the "outside contraction point" 
                %     if it's better than the "reflection point"    
                %     Otherwise, prepare to perform a "shrink"
                    how = 'shrink';

%% Main algorithm: iterate until 
% (a) the maximum coordinate difference between the current best point and the 
% next best ncp other points in the simplex is less than or equal to simplexOut.tolx. 
% Specifically, until max(||v2-v1||, ||v2-v1||, ..., ||v(ncp+1)-v1||) <= TolX,
% where ||.|| is the infinity-norm, and v1 holds the 
% vertex with the current lowest function value; AND (Cannot use OR instead of AND.)
% (b) the corresponding difference in function values is less than or equal
% to TolFun. ; OR
% (c) the maximum number of iterations is exceeded; OR
% (d) the maximum number of function evaluations is exceeded

    simplexOut.tolx = max(v(:, 1)) * relativeParamTolerance;    % parameter change tolerance
    simplexOut.maxParamChange = max(max(abs(v(:, 2:ncp+1) - v(:, ones(1, ncp)))));

    ncp = min(2, n);                % number of points to compare against %%% Why min(2, n)?

        initymax = fv(:, 1) * 1.1;

    plot(outparams.simplexIterCount + 2, err0.totalError, ...
    plot(outparams.simplexIterCount + 1 + ctIterations, fv(:, 1), ...

%    if max(abs(fv(1)-fv(two2np1))) <= max(tolf,10*eps(fv(1))) && ...
%            max(max(abs(v(:,two2np1)-v(:,onesn)))) <= max(tolx,10*eps(max(v(:,1))))
%        break
%    end

    fprintf('ctIterations == %d\n', ctIterations);
    fprintf('how == %s\n\n', how);

% 
% % If just 'defaults' passed in, return the default options in X
% if nargin==1 && nargout <= 1 && isequal(funfcn,'defaults')
%     x = defaultopt;
%     return
% end
% 
% if nargin<3, options = []; end
% 
% % Detect problem structure input
% if nargin == 1
%     if isa(funfcn,'struct') 
%         [funfcn,x,options] = separateOptimStruct(funfcn);
%     else % Single input and non-structure
%         error('MATLAB:fminsearch:InputArg','The input to FMINSEARCH should be either a structure with valid fields or consist of at least two arguments.');
%     end
% end
% 
% if nargin == 0
%     error('MATLAB:fminsearch:NotEnoughInputs',...
%         'FMINSEARCH requires at least two input arguments');
% end
% 
% 
% % Check for non-double inputs
% if ~isa(x,'double')
%   error('MATLAB:fminsearch:NonDoubleInput', ...
%          'FMINSEARCH only accepts inputs of data type double.')
% end

% % Convert to function handle as needed.
% funfcn = fcnchk(funfcn,length(varargin));
% % Add a wrapper function to check for Inf/NaN/complex values
% if simplexOut.funValCheck
%     % Add a wrapper function, CHECKFUN, to check for NaN/complex values without
%     % having to change the calls that look like this:
%     % f = funfcn(x,varargin{:});
%     % x is the first argument to CHECKFUN, then the user's function,
%     % then the elements of varargin. To accomplish this we need to add the 
%     % user's function to the beginning of varargin, and change funfcn to be
%     % CHECKFUN.
%     varargin = {funfcn, varargin{:}};
%     funfcn = @checkfun;
% end

% Initialize the output and plot functions.
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'init',ctIterations, ...
        ctEvals, how, fv(:,1),varargin{:});
    if stop
        [x,fval,exitFlag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',ctIterations, ...
        ctEvals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitFlag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',ctIterations, ...
        ctEvals, how, fv(:,1),varargin{:});
    if stop  % Stop per user request.
        [x,fval,exitFlag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
        if  prnt > 0
            disp(output.message)
        end
        return;
    end
end
    % OutputFcn and PlotFcns call
    if haveoutputfcn || haveplotfcn
        [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',ctIterations, ...
            ctEvals, how, fv(:,1),varargin{:});
        if stop  % Stop per user request.
            [x,fval,exitFlag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
            if  prnt > 0
                disp(output.message)
            end
            return;
        end
    end
    %%% END ORIGINAL CODE

% OutputFcn and PlotFcns call
if haveoutputfcn || haveplotfcn
    callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,'done',ctIterations, ctEvals, how, fval, varargin{:});
end

%--------------------------------------------------------------------------
function [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,x,xOutputfcn,state,iter,...
    numf,how,f,varargin)
% CALLOUTPUTANDPLOTFCNS assigns values to the struct OptimValues and then calls the
% outputfcn/plotfcns.
%
% state - can have the values 'init','iter', or 'done'.

% For the 'done' state we do not check the value of 'stop' because the
% optimization is already done.
optimValues.iteration = iter;
optimValues.funccount = numf;
optimValues.fval = f;
optimValues.procedure = how;

xOutputfcn(:) = x;  % Set x to have user expected size
stop = false;
% Call output functions
if ~isempty(outputfcn)
    switch state
        case {'iter','init'}
            stop = callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimOutputFcns(outputfcn,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end
% Call plot functions
if ~isempty(plotfcns)
    switch state
        case {'iter','init'}
            stop = callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:}) || stop;
        case 'done'
            callAllOptimPlotFcns(plotfcns,xOutputfcn,optimValues,state,varargin{:});
        otherwise
            error('MATLAB:fminsearch:InvalidState', ...
                'Unknown state in CALLOUTPUTANDPLOTFCNS.')
    end
end

%--------------------------------------------------------------------------
function [x,FVAL,EXITFLAG,OUTPUT] = cleanUpInterrupt(xOutputfcn,optimValues)
% CLEANUPINTERRUPT updates or sets all the output arguments of FMINBND when the optimization
% is interrupted.

x = xOutputfcn;
FVAL = optimValues.fval;
EXITFLAG = -1;
OUTPUT.iterations = optimValues.iteration;
OUTPUT.funcCount = optimValues.funccount;
OUTPUT.algorithm = 'golden section search, parabolic interpolation';
OUTPUT.message = 'Optimization terminated prematurely by user.';

%--------------------------------------------------------------------------
function f = checkfun(x,userfcn,varargin)
% CHECKFUN checks for complex or NaN results from userfcn.

f = userfcn(x,varargin{:});
% Note: we do not check for Inf as FMINSEARCH handles it naturally.
if isnan(f)
    error('MATLAB:fminsearch:checkfun:NaNFval', ...
        'User function ''%s'' returned NaN when evaluated;\n FMINSEARCH cannot continue.', ...
        localChar(userfcn));  
elseif ~isreal(f)
    error('MATLAB:fminsearch:checkfun:ComplexFval', ...
        'User function ''%s'' returned a complex value when evaluated;\n FMINSEARCH cannot continue.', ...
        localChar(userfcn));  
end

%--------------------------------------------------------------------------
function strfcn = localChar(fcn)
% Convert the fcn to a string for printing

if ischar(fcn)
    strfcn = fcn;
elseif isa(fcn,'inline')
    strfcn = char(fcn);
elseif isa(fcn,'function_handle')
    strfcn = func2str(fcn);
else
    try
        strfcn = char(fcn);
    catch
        strfcn = '(name not printable)';
    end
end
%--------------------------------------------------------------------------


% Handle the output
outputfcn = optimget(options, 'OutputFcn', defaultopt, 'fast');
if isempty(outputfcn)
    haveoutputfcn = false;
else
    haveoutputfcn = true;
    xOutputfcn = x; % Last x passed to outputfcn; has the input x's shape
    % Parse OutputFcn which is needed to support cell array syntax for OutputFcn.
    outputfcn = createCellArrayOfFunctions(outputfcn,'OutputFcn');
end

% Handle the plot
%plotfcns = optimget(options,'PlotFcns',defaultopt,'fast');
if isempty(plotfcns)
    haveplotfcn = false;
else
    haveplotfcn = true;
    xOutputfcn = x; % Last x passed to plotfcns; has the input x's shape
    % Parse PlotFcns which is needed to support cell array syntax for PlotFcns.
    plotfcns = createCellArrayOfFunctions(plotfcns,'PlotFcns');
end

function [optimizedparams, simplexOut] = m3ha_fminsearch3(realData, outparams, outFolderName, tolfun)

%err0 = run_neuron_once_2comp(realData,outparams,outFolderName);

% outparams.neuronparams = reshape(xtransform(x,bounds),[],1);
paramnames = outparams.neuronparamnames;

        fprintf(fid, '%2.2g, ', numel(outparams.swpw)); %1
        fprintf(fid, repmat('%6.4e, ', 1, nParams), reshape(xtransform(bounds.p,bounds), [], 1)); %nParams

n = numel(x0);
nParams = n;
xin = x0(:);        % Force pin to be a column vector
x(:) = x0;        % Change x to the form expected by funfcn
%fv(:,1) = funfcn(x,varargin{:});  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % Log header
    fid = fopen(fullfile(outFolderName, logfilename), 'w');        % append
%    fprintf(fid, '%s, ', 'Experiment name'); %1
    fprintf(fid, '%s, ', 'Number of sweeps'); %1
    fprintf(fid, '%s, ', 'Total error'); %1
    fprintf(fid, '%s, ', 'LTS to sweep error ratio'); %1
    fprintf(fid, repmat('%s, ', 1, nParams), paramnames{:}); %nParams
%    fprintf(fid, repmat('%2.2g, ', 1, n), outparams.sortedswpnum); %n
    fprintf(fid, repmat('%s, ', 1, n), repmat({'sweep weight'}, 1, n)); %n
    fprintf(fid, repmat('%s, ', 1, n), repmat({'sweep left edge'}, 1, n)); %n
    fprintf(fid, repmat('%s, ', 1, n), repmat({'sweep right edge'}, 1, n)); %n
    fprintf(fid, repmat('%s, ', 1, n), repmat({'sweep error'}, 1, n)); %n
    fprintf(fid, repmat('%s, ', 1, 3), 'LTS amp weight', 'LTS time weight', 'LTS slope weight'); %3
%    fprintf(fid, '%s, %s, %s, ', 'LTS amp error', 'LTS time error', 'LTS slope error'); %3
    fprintf(fid, '%s, %s, %s, %s', 'ctEvals', 'maxFunctionEvaluations', 'ctIterations', 'maxIterations'); %4
    fprintf(fid, '\n');
    fclose(fid);

    % save initial params and error
    starterr = run_neuron_once_4compgabab(realData, outparams, outFolderName);

    fid = fopen(fullfile(outFolderName, logfilename), 'a');        % append
%    fprintf(fid, '%s, ',outparams.experimentname); %1
    fprintf(fid, '%2.2g, ', n); %1
    fprintf(fid, '%6.4e, ', starterr.totalError); %1
    fprintf(fid, '%6.4f, ', outparams.lts_to_swp_errratio); %1
    fprintf(fid, repmat('%6.4e, ', 1, nParams), outparams.neuronparams); %nParams
%    fprintf(fid, repmat('%2.2g, ', 1, n), outparams.sortedswpnum); %n
    fprintf(fid, repmat('%6.4f, ', 1, n), outparams.swpw); %n
    fprintf(fid, repmat('%6.4f, ', 1, n), outparams.swpedges(:,1)'); %n
    fprintf(fid, repmat('%6.4f, ', 1, n), outparams.swpedges(:,2)'); %n
    fprintf(fid, repmat('%6.4e, ', 1, n), starterr.swperr); %n
    fprintf(fid, repmat('%6.4f, ', 1, 3), outparams.ltsWeights); %3
%    fprintf(fid, '%6.4e, %6.4e, %6.4e, ', err.ltserr_maxamp_val,...
%            err.ltserr_maxamp_time, err.ltserr_maxdiffamp_val); %3
    fprintf(fid, '%6.4e, %6.4e, %6.4e, %6.4e', 0, maxFunctionEvaluations, 0, maxIterations);
    fprintf(fid, '\n');
    fclose(fid);

printType = optimget(options, 'Display', defaultopt, 'fast');
tolx = optimget(options, 'TolX', defaultopt, 'fast');
tolf = optimget(options, 'TolFun', defaultopt, 'fast');
maxFunctionEvaluations = optimget(options, 'MaxFunEvals', defaultopt, 'fast');
maxIterations = optimget(options, 'MaxIter', defaultopt, 'fast');
funValCheck = strcmp(optimget(options, 'FunValCheck', defaultopt, 'fast'), 'on');

%simplexOut.printType = getfield(options,'Display');
%simplexOut.tolx = getfield(options,'TolX');
%simplexOut.tolf = getfield(options,'TolFun');
%simplexOut.maxFunctionEvaluations = getfield(options,'MaxFunEvals');
%simplexOut.maxIterations = getfield(options,'MaxIter');
%simplexOut.funValCheck = strcmp(getfield(options,'FunValCheck'),'on');

% ctEvals = 1;
% ctIterations = 0;

xlswrite([outFolderName, '/simplexrun_', num2str(outparams.simplexNum), '.xls'], bestparams);

%options = defaultopt;
    %optimset not working due to library loading issue 2014-04-49
    %options=optimset(options,outparams.simplexParamNames{k},outparams.simplexParams(k));

%outputfcn = getfield(options,'OutputFcn');
%plotfcns = getfield(options,'PlotFcns');

%x = x0;

%     % OutputFcn and PlotFcns call
%     if haveoutputfcn || haveplotfcn
%         [xOutputfcn, optimValues, stop] = callOutputAndPlotFcns(outputfcn,plotfcns,v(:,1),xOutputfcn,'iter',ctIterations, ...
%             ctEvals, how, fv(:,1),varargin{:});
%         if stop  % Stop per user request.
%             [x,fval,exitFlag,output] = cleanUpInterrupt(xOutputfcn,optimValues);
%             if  prnt > 0
%                 disp(output.message)
%             end
%             return;
%         end
%     end

%outparams.neuronparams = x;
%outparams.neuronparams = x;
    %f = funfcn(x,varargin{:});

LB = bounds.LB;
UB = bounds.UB;

%p0 = x0;                        % initial parameter values unconstrained
%n = numel(x);
%xin = x(:);        % Force xin to be a column vector
    ctIterations = 0;
    ctEvals = 0;
    simplexOut.error = 0;


            if x0(i) <= LB(i)    % infeasible starting value
                % Make -90 degrees
                p0(k) = -pi/2;
            elseif x0(i) >= UB(i)    % infeasible starting value
                % Make +90 degrees
                p0(k) = pi/2;
            else
                % First transform [LB UB] linearly to [-1 1]
                p0_n121(k) = 2 * (x0(i) - LB(i))/(UB(i) - LB(i)) - 1;

                % Then transform to (-infty infty)
                % shift by 2*pi to avoid problems at zero in fminsearch
                % otherwise, the initial simplex is vanishingly small
                p0(k) = 2 * pi + asin(p0_n121(k));
%                p0(k) = 2 * pi+asin(max(-1, min(1, p0_n121(k))));    % not necessary 
            end

%% Transform starting values into unconstrained surrogates
% k allows some variables to be fixed, thus dropped from the optimization

    % Find initial errors and params
    starterr = run_neuron_once_4compgabab(realData, outparams, outFolderName);
    m3ha_log_errors_params(outFolderName, logfilename, outparams, starterr, simplexOut);

newopt = defaultopt;
options = newopt;

% If any parameters were fixed, shorten p0_n121, p0_n121corr & p0
if ct < nParams
    p0_n121((ct+1):nParams) = [];
    p0_n121corr((ct+1):nParams) = [];
    p0((ct+1):nParams) = [];
end
if isempty(p0)        % all variables were fixed

onesn = ones(1, n);
optimizedparams = outparams.neuronparams;


% Compute maximum error change, maximum parameter change and respective tolerances
'TolX', 1e-4, ...
'TolFun', 1e-4, ...
simplexOut.tolx = optimget(options, 'TolX', defaultopt, 'fast');
simplexOut.tolf = optimget(options, 'TolFun', defaultopt, 'fast');
simplexOut.maxErrorChange = max(abs(fv(1) - fv(2:ncp+1)));
simplexOut.tolf = fv(1) * simplexOut.relativeErrorTolerance;    % error tolerance 
simplexOut.tolx = simplexOut.relativeParamTolerance;        % parameter change tolerance (relative)
simplexOut.firsttolf = simplexOut.tolf;
simplexOut.firsttolx = simplexOut.tolx;
&& (simplexOut.maxErrorChange > simplexOut.tolf || simplexOut.maxParamChange > simplexOut.tolx)
simplexOut.maxErrorChange = simplexOut.tolf*2;    % initialize to greater than simplexOut.tolf so that simplex will start
simplexOut.maxParamChange = simplexOut.tolx*2;    % initialize to greater than simplexOut.tolx so that simplex will start

simplexOut.printType = optimget(options, 'Display', defaultopt, 'fast');
simplexOut.maxFunctionEvaluations = optimget(options, 'MaxFunEvals', defaultopt, 'fast');
simplexOut.maxIterations = optimget(options, 'MaxIter', defaultopt, 'fast');

options = optimset(defaultopt, 'Display', displayoption);
simplexOut.funValCheck = strcmp(getfield(options, 'FunValCheck'), 'on');

            'FunValCheck', 'off', ...    % TODO: Examine
            'OutputFcn', [], ...        % TODO: Examine
            'PlotFcns', []);        % TODO: Examine


%% Update simplexOut even if autoLoggerFlag is not true
simplexOut = update_simplexOut(simplexOut, how, ctIterations, ctEvals, fv, v, bounds);

simplexOut.maxErrorChange = simplexOut.relativeErrorTolerance * 2;    % initialize to greater than simplexOut.tolf
simplexOut.maxParamChange = simplexOut.relativeParamTolerance * 2;    % initialize to greater than simplexOut.tolx

save(fullfile(outparams.outFolderName, ...
        [outparams.prefix, '_simplexrun_', num2str(outparams.simplexNum), '.mat']), ...
        'bestparams', 'simplexOut', 'outparams', '-v7.3');

    msg = sprintf(['Optimization terminated:\n', ...
            ' the current x satisfies the termination criteria using OPTIONS.TolX of %e \n' ...
            ' and F(X) satisfies the convergence criteria using OPTIONS.TolFun of %e \n'], ...
            simplexOut.tolx, simplexOut.tolf);

plot(outparams.simplexIterCount + ctIterations, err0.totalError, ...
    'o', 'color', cm(mod(outparams.simplexNum, size(cm, 1)) + 1, :)); hold on;
    plot(iter, err.totalError, '.', 'color', cm(mod(outparams.simplexNum, size(cm, 1)) + 1, :));

    ylim auto;            % reset ylim to auto to enable further expansion

% transform initial values into the range [3*pi/2 5*pi/2] nonlinearly using arcsin
            % Then transform to [3*pi/2 5*pi/2] nonlinearly with arcsin
            p0(ct) = 2*pi + asin(p0_n121corr(ct));
            % First transform back to [0 1]
            p_01(ct) = (sin(p(ct)) + 1) / 2;


    if y(j) ~= 0        % if the parameter is not zero
        y(j) = (1 + simplexOut.usualDelta)*y(j);        % increment parameter by usualDelta (relative)
    else            % if the parameter is zero
        y(j) = simplexOut.zero_term_delta;            % increment parameter by zero_term_delta (absolute)
    end

            'ZeroTermDelta', 0.00025, ...    % absolute increment for zero parameters
            simplexOut.zero_term_delta = getfield(options, 'ZeroTermDelta');

            % Then transform to [3*pi/2 5*pi/2] nonlinearly with arctan
            % Shift by 2*pi to avoid problems at zero in fminsearch
            % Otherwise, the initial simplex is vanishingly small
            p0(ct) = 2*pi + atan(p0_n121corr(ct));

    y = p0;            % start with initial set of parameters; all values are within [3*pi/2 5*pi/2]

    simplexOut.maxParamChange = max(max(abs(v(:, 2:ncp+1) - v(:, ones(1, ncp))) ./ abs(v(:, ones(1, ncp)))));

%% Hard-coded parameters
saveParamsFlag = false;     % whether to save parameters as a .p file

figure(fig);

bestparams = [outparams.neuronparamnames', num2cell(simplexOut.neuronParamsTable'), ...
    num2cell(outparams.neuronparams_min'), num2cell(outparams.neuronparams_max')];
save(fullfile(outparams.outFolderName, [outparams.prefix, '_best.p']), 'bestparams');

%% Add directories to search path for required functions
if exist('/home/Matlab/', 'dir') == 7
    functionsdirectory = '/home/Matlab/';
elseif exist('/scratch/al4ng/Matlab/', 'dir') == 7
    functionsdirectory = '/scratch/al4ng/Matlab/';
else
    error('Valid functionsdirectory does not exist!');
end
if ~isdeployed
    addpath(fullfile(functionsdirectory, '/Adams_Functions/'));        
                                % for set_fields_zero.m, restore_fields.m,
end

%       cd/run_neuron_once_4compgabab.m

outparams.logconcisefilename = strrep(outparams.logfilename, '.csv', '_concise.csv');

save_params(sheetPath, ...
            outparams.neuronparamnames, simplexOut.neuronParamsTable, ...
            outparams.neuronparams_min, outparams.neuronparams_max);

pIsLog = outparams.neuronparamislog;

[err, ~, ~, ~] = m3ha_run_neuron_once(outparams, [], ...
                                        'RealData', realData, ...
                                        'SaveParamsFlag', false, ...
                                        'OnHpcFlag', onHpcFlag);
% Save original prefix
prefixOrig = outparams.prefix;

% Restore original prefix
prefix = prefixOrig;

% Do not plot sweeps for other function evaluations
outparams = ...
    set_fields_zero(outparams, 'plotIndividualFlag', 'plotOverlappedFlag');

% Plot sweeps if user requests so
if outparams.plotIndividualFlag || outparams.plotOverlappedFlag
end

[outparams] = restore_fields(outparams, 'plotIndividualFlag', 'plotOverlappedFlag');    % see if user wants to plot sweeps
prefixOrig = prefix;                        % save original prefix
prefix = [prefix, '_bef'];                % change prefix for "before fitting"
[err0, ctEvals] = eval_with_bounds(p0, bounds, pIsLog, ctEvals, outparams);
prefix = prefixOrig;                        % restore original prefix
[outparams] = set_fields_zero(outparams, 'plotIndividualFlag', 'plotOverlappedFlag');    % do not plot sweeps for other function evaluations

nParams = numel(outparams.neuronparams);

simplexOut.paramsInUseNames{ct} = outparams.neuronparamnames{i};

% Determine indices of parameters in use
indInUse = neuronParamsTable{:, 'InUse'};

p_01 = zeros(1, numinuse);    % stores parameters transformed to [0 1]
xpart = zeros(1, numinuse);    % stores transformed parameters
x = zeros(1, nParams);    % stores transformed parameters and original parameters

[simplexOut.neuronParamsTable, simplexOut.paramsInUseValues] = ...
    xtransform(v(:, 1), bounds, pIsLog);

%}