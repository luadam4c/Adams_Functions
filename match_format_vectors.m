function [output1] = match_format_vectors (varargin)
%% Match the format of individual vectors by making them all column vectors with the same length
% Usage: [output1] = match_format_vectors (varargin)
% Explanation:
%       TODO
% Example(s):
%       TODO
% Outputs:
%       output1     - TODO: Description of output1
%                   specified as a TODO
% Arguments:
%       vec1     - TODO: Description of vec1
%                   must be a TODO
%       varargin    - 'param1': TODO: Description of param1
%                   must be a TODO
%                   default == TODO
%
% Requires:
%       /TODO:dir/TODO:file
%
% Used by:
%       /TODO:dir/TODO:file

% File History:
% 2018-12-16 Created by Adam Lu

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% If there are more outputs requested than inputs, return error
if nargout > nargin
    error(['Cannot request more outputs than provided inputs, ', ...
            'type ''help %s'' for usage'], mfilename);
end


%% Preparation


%% Do the job
% Place all arguments in an input structure
%   Note: varargin is a row cell array
myStructInputs = cell2struct(varargin, myFieldNames, 2);

%% Output results
% 


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
