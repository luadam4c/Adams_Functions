function [nRows, nCols] = check_consistency (data1, data2, isSpikeData)
%% Checks consistency of data dimensions between two datasets
% Usage: [nRows, nCols] = check_consistency (data1, data2, isSpikeData)
% Explanation:
%       Verifies that data1 and data2 have compatible dimensions (channels/trials).
%       Stops with an error if dimensions do not match.
%
% Outputs:
%       nRows   - Number of samples (left empty if struct)
%       nCols   - Number of channels/trials
%
% Arguments:
%       data1       - First dataset (matrix or struct)
%       data2       - Second dataset (matrix or struct)
%       isSpikeData - (opt) 1 if one dataset is spike times (1d array)
%                   default == 0
%
% Requires:
%
% Used by:
%       coherencyc.m
%       coherencycpt.m
%
% File History:
% 2026-01-13 Reorganized and reannotated by Gemini from check_consistency.m
%   in Chronux version 2.12 v03 package Latest release (2018-Oct-25)
%   \\moorelaboratory.dts.usc.edu\Shared\Code\Downloaded_Functions\chronux_2_12

%% Deal with arguments
if nargin < 3 || isempty(isSpikeData)
    isSpikeData = 0;
end

nRows1 = []; 
nRows2 = [];

%% Do the job

% Get dimensions for Data 1
if isstruct(data1)
    nCols1 = length(data1);
else
    [nRows1, nCols1] = size(data1);
end

% Get dimensions for Data 2
if isstruct(data2)
    nCols2 = length(data2);
else
    [nRows2, nCols2] = size(data2);
end

% Check Channel/Trial consistency
if nCols1 ~= nCols2
    error('inconsistent dimensions: Number of channels/trials do not match');
end

% Check Sample consistency (if not spike data and not structs)
if isSpikeData == 0
    if ~isstruct(data1) && ~isstruct(data2)
        if nRows1 ~= nRows2
            error('inconsistent dimensions: Number of time points do not match');
        end
    end
end

%% Output results
nRows = nRows1; 
nCols = nCols1;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%