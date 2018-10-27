function variableName = get_variable_name (variable)
%% Returns a variable's name as a string
% Usage: variableName = get_variable_name (variable)
% Explanation:
%       get_variable_name (myVar) == 'myVar'
% Example(s):
%       get_variable_name (myVar)
% Outputs:
%       variableName    - name of variable
%                       specified as a character vector
% Arguments:    
%       variable        - variable in workspace
%
% Used by:    
%       cd/compute_weighted_average.m

% File History:
% 2018-10-26 Adapted from https://www.mathworks.com/matlabcentral/answers/
%               382503-how-can-i-get-the-name-of-a-matlab-variable-as-a-string
% 

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Do the job
variableName = inputname(1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%