function [xData, yData] = extract_data_from_lines(varargin)
%% Extracts data from line objects on an axes, with optional data override
% Usage: [xData, yData] = extract_data_from_lines(varargin)
% Explanation:
%       This function finds all objects of type 'Line' on a given axes,
%       extracts their XData and YData properties, and concatenates them
%       into single column vectors. The user can also provide their own
%       XData or YData to replace the extracted values.
%
% Example(s):
%       x1 = 1:10; y1 = x1 + randn(size(x1));
%       x2 = 1:12; y2 = x2 * 1.5 + randn(size(x2));
%       y2(3) = NaN;
%       figure;
%       plot(x1, y1, 'o');
%       hold on;
%       plot(x2, y2, 's');
%       [x, y] = extract_data_from_lines('RemoveNaN', true);
%       disp(['Number of data points with NaNs removed: ', num2str(numel(x))]);
%
%       % Standard extraction
%       figure; plot(1:10, (1:10) + randn(1, 10), 'o');
%       [x, y] = extract_data_from_lines;
%
%       % Override YData with a known vector
%       yNew = linspace(0, 1, numel(x));
%       [xFinal, yFinal] = extract_data_from_lines('YData', yNew);
%       figure; plot(xFinal, yFinal, 'rs-'); title('Data with Y-values replaced');
%
% Outputs:
%       xData           - all x data from the axes' line objects
%                       specified as a numeric column vector
%       yData           - all y data from the axes' line objects
%                       specified as a numeric column vector
%
% Arguments:
%       varargin    - 'AxesHandle': axes handle to extract data from
%                   must be a empty or an axes object handle
%                   default == gca
%                   - 'XData': user-provided x-data to replace extracted data
%                   must be a numeric vector
%                   default == []
%                   - 'YData': user-provided y-data to replace extracted data
%                   must be a numeric vector
%                   default == []
%                   - 'RemoveNaN': whether to remove NaN values
%                   must be a logical scalar
%                   default == false
%
% Requires:
%       cd/apply_over_cells.m
%       cd/extract_fields.m
%       cd/force_column_vector.m
%       cd/set_axes_properties.m
%
% Used by:
%       cd/plot_correlation_coefficient.m
%       cd/plot_regression_line.m

% File History:
% 2025-09-17 - Extracted from plot_regression_line.m by Gemini
% 2025-09-17 - Added 'XData' and 'YData' optional arguments to override extraction

%% Hard-coded parameters
% None

%% Default values for optional arguments
axHandleDefault = [];           % gca by default
xDataDefault = [];
yDataDefault = [];
removeNaNDefault = false;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Deal with arguments
% Set up Input Parser Scheme
iP = inputParser;
iP.FunctionName = mfilename;

% Add parameter-value pairs to the Input Parser
addParameter(iP, 'AxesHandle', axHandleDefault);
addParameter(iP, 'XData', xDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'YData', yDataDefault, ...
    @(x) validateattributes(x, {'numeric'}, {'2d'}));
addParameter(iP, 'RemoveNaN', removeNaNDefault, ...
    @(x) validateattributes(x, {'logical', 'numeric'}, {'scalar'}));

% Read from the Input Parser
parse(iP, varargin{:});
axHandle = iP.Results.AxesHandle;
xDataUser = iP.Results.XData;
yDataUser = iP.Results.YData;
removeNaN = iP.Results.RemoveNaN;

%% Preparation
% Decide on the axes to use
axHandle = set_axes_properties('AxesHandle', axHandle);

% Find all line objects on the axes
% Note: plot, scatter, etc. create 'Line' objects for data points
lines = findobj(axHandle, 'Type', 'Line');

%% Extraction
% Extract x data values unless provided by the user
if ~isempty(xDataUser)
    xData = force_column_vector(xDataUser);
else
    if isempty(lines)
        xData = [];
    else
        xData = extract_fields(lines, 'XData');
        xData = force_column_vector(xData, 'ToLinearize', true);
        xData = apply_over_cells(@vertcat, xData);
    end
end

% Extract y data values unless provided by the user
if ~isempty(yDataUser)
    yData = force_column_vector(yDataUser);
else
    if isempty(lines)
        yData = [];
    else
        yData = extract_fields(lines, 'YData');
        yData = force_column_vector(yData, 'ToLinearize', true);
        yData = apply_over_cells(@vertcat, yData);
    end
end

%% Validation
% Ensure that xData and yData have the same number of elements
if ~isempty(xData) && ~isempty(yData) && numel(xData) ~= numel(yData)
    error('XData and YData must have the same number of elements!');
end

%% Post-processing
% Return empty if no data
if isempty(xData) && isempty(yData)
    return;
end

% Remove NaN values if requested
if removeNaN && (any(isnan(xData)) || any(isnan(yData)))
    % Find rows where either x or y is NaN
    toRemove = isnan(xData) | isnan(yData);
    
    % Remove those rows from both vectors
    xData(toRemove) = [];
    yData(toRemove) = [];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%