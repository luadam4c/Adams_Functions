function myTable = updatevars (myTable, varValue)
%% Replace a variable in a table or add it if it doesn't exist
% Usage: myTable = updatevars (myTable, varValue)
% TODO: Add input parser
% Used by:
%       cd/m3ha_simulate_population.m
%
% File History:
%   2020-08-04 Moved from m3ha_simulate_population.m

varName = inputname(2);
if is_field(myTable, varName)
    myTable.(varName) = varValue;
else
    myTable = addvars(myTable, varValue, 'NewVariableNames', varName);
end