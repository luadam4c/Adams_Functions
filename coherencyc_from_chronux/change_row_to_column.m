function outputData = change_row_to_column (inputData)
%% Transforms 1D arrays into column vectors needed by Chronux routines
% Usage: outputData = change_row_to_column (inputData)
% Explanation:
%       Checks if the input is a vector or a structure array. 
%       If it is a vector, it forces it into a column (Nx1) orientation.
%       If it is a structure, it extracts the first field and linearizes it 
%       if the structure array itself is 1D.
%
% Example(s):
%       out = change_row_to_column([1 2 3]) -> [1; 2; 3]
%
% Outputs:
%       outputData  - The data transformed into form samples x channels/trials
%
% Arguments:
%       inputData   - Data to transform
%                   must be a matrix, vector, or structure
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from change_row_to_column.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 1
    error('Need inputData');
end

%% Do the job
dataTemp = [];

% Check if input is a structure
if isstruct(inputData)
    numStructElements = length(inputData);
    
    if numStructElements == 1
        % Extract field names
        fieldNames = fieldnames(inputData);
        
        % extract the first field's data
        % eval(['dtmp=data.' fnames{1} ';']) -> Replaced with dynamic field access
        dataTemp = inputData.(fieldNames{1});
        
        % Linearize to column
        outputData = dataTemp(:);
    else
        % If structure array has multiple elements, return as is
        outputData = inputData;
    end
else
    % Input is a matrix or vector
    [nRows, nCols] = size(inputData);
    
    if nRows == 1 || nCols == 1
        % Linearize to column
        outputData = inputData(:);
    else
        % Return matrix as is
        outputData = inputData;
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%