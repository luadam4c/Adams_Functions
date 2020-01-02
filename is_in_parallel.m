function answer = is_in_parallel ()
%% Checks whether in a parfor loop
% Usage: answer = is_in_parallel ()
% Explanation:
%       % TODO
% Example(s):
%       is_in_parallel
%       parfor i = 1:10
%            test(i) = is_in_parallel;
%       end
% Outputs:
%       answer      - whether in a parallel loop
%                   specified as a logical scalar
% Arguments:
%
% Used by:
%       cd/array_fun.m
%       cd/load_neuron_outputs.m

% File History:
% 2018-11-01 Copied from https://www.mathworks.com/matlabcentral/answers/
%               58228-am-i-running-in-parallel-best-way-to-check
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
try
    answer = ~isempty(getCurrentTask());
catch err
    if ~strcmp(err.identifier, 'MATLAB:UndefinedFunction')
        rethrow(err);
    end
    answer = false;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%