%% Loads example data structures for testing
%
% Requires:
%       cd/distribute_balls_into_boxes.m
%
% Used by:    
%       cd/test_passive_fit.m

% File History:
% 2018-XX-XX Created by Adam LU
% 2018-06-12 Moved xlsData from Test_xlWrite.m

%% Hard-coded parameters
nRows = 3;
nCols = 3;
nSamples = 300;
siMs = 0.1;
tau1 = 10;
tau2 = 1;
amp1 = 2;
amp2 = 5;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Logical arrays
myLogicalScalar = true;
myLogicalRow = [true, false, true, true, false];
myLogicalCol = [false; false; true; false; false];
myLogical2D = [true, false; false, true];

% Character arrays
myString = 'I love you';

% Numeric arrays
myNumericScalar = nCols;
myNumericCol = ones(nRows, 1);
myNumericRow = ones(1, nCols);
myNumeric2D = magic(nRows);
myNumeric3D = floor(rand(3, 4, 2) * 100);

% Cell arrays of character arrays
names = {'Adam', 'Ashley', 'Mark', 'Peter', 'Katie'};
myCellStrScalar = {'I love you'};
myCellStrRow = {'I love you', 'You love me', 'Blab hooray!', 'Why?'};
myCellStrCol = {'I love you'; 'You love me'; 'Blab hooray!'; 'Why?'};
myCellStr2D = {'I love you', 'You love me'; 'Blab hooray!', 'Why?'};

% Cell arrays of numeric arrays
myCellNumeric2D = {myNumeric2D; myNumeric2D * 2; myNumeric2D * 0.5};
myCellNumeric3D = {myNumeric3D; myNumeric3D * 2; myNumeric3D * 0.5};
myCellColumnVecs = {myNumericCol; myNumericCol * 2; myNumericCol * 0.5};
myCellRowVecs = {myNumericRow; myNumericRow * 2; myNumericRow * 0.5};

% Cell arrays of cell arrays
myCellCell = {myCellStrRow, myCellStrRow, myCellStrRow};

% Structures
blab = struct;
blab.students = names;
blab.isMarried = myLogicalCol;
myStruct = struct;
myStruct.myLogicalScalar = myLogicalScalar;
myStruct.myLogicalCol = myLogicalCol;
myStruct.myLogicalRow = myLogicalRow;
myStruct.myLogical2D = myLogical2D;
myStruct.myNumericScalar = myNumericScalar;
myStruct.myNumericCol = myNumericCol;
myStruct.myNumericRow = myNumericRow;
myStruct.myNumeric2D = myNumeric2D;
myStruct.names = names;
myStruct.myCellStrScalar = myCellStrScalar;
myStruct.myCellStrRow = myCellStrRow;
myStruct.myCellStrCol = myCellStrCol;
myStruct.myCellStr2D = myCellStr2D;
myStruct.blab = blab;
myScalarStruct.a = floor(rand(1) * 10);
myScalarStruct.b = floor(rand(1) * 10);

% Structure arrays
myStructArray = [myStruct, myStruct, myStruct];
for i = 1:10
    myScalarStructArray(i).a = floor(rand(1) * 10);
    myScalarStructArray(i).b = floor(rand(1) * 10);
end

% Sheet Data
xlsData = {'A Number' 'Boolean Data' 'Empty Cells' 'Strings';...
    1 true [] 'String Text';...
    5 false [] 'Another very descriptive text';...
    -6.26 false 'This should have been an empty cell but I made an error' 'This is text';...
    1e8 true [] 'Last cell with text';...
    1e3 false NaN NaN;...
    1e2 true [] 'test'};

% Normally distributed data
normData = randn(nSamples, 1);

% Normally distributed data with outliers
normDataWOutliers = [normData; -100; 200; 1000];

% Time vector
myTimeVec = transpose(1:nSamples) * siMs;

% A random time series
myRandomSignal = rand(nSamples, 1);
myRandomSignals1 = rand(nSamples, 3) + repmat(1:3, [nSamples, 1]);
myRandomSignals2 = (-1 + rand(nSamples, 20) * 2) + ...
                    repmat(sin(myTimeVec), [1, 20]);
myRandomSignals3 = force_column_numeric(myRandomSignals2);
myRandomTruth2 = repmat(sin(myTimeVec), [1, 20]);
myRandomTruth3 = force_column_numeric(myRandomTruth2);

% Exponentials
myRisingExponential = 5*(1-exp(-myTimeVec/2));
myFallingExponential = 5*(exp(-myTimeVec/2));

% Pulses
nSamplesEachBin1 = distribute_balls_into_boxes(nSamples, 3, 'evenly', true);
nBin1 = nSamplesEachBin1(1);
nShift = floor(nBin1 / 2);
nBin2 = nSamplesEachBin1(2);
nBin3 = nSamplesEachBin1(3);
nSamplesEachBin2 = [nBin1 - nShift, nBin2 - nShift, nBin3 + 2*nShift];
myPulse = [zeros(nBin1, 1); ones(nBin2, 1); zeros(nBin3, 1)];
myPulse1 = [zeros(nBin1, 1); ones(nBin2, 1); zeros(nBin3, 1)];
myPulse2 = [zeros(nSamplesEachBin2(1), 1); ...
            (-1) * ones(nSamplesEachBin2(2), 1); ...
            zeros(nSamplesEachBin2(3), 1)];
myPulse3 = [ones(nBin1 + nShift, 1); -0.5*ones(nBin2, 1); ones(nBin3 - nShift, 1)];
myPulse4 = [zeros(nBin1 - nShift, 1); 2*ones(nBin2, 1); zeros(nBin3 + nShift, 1)];
myPulses = [myPulse1, myPulse2, myPulse3, myPulse4];
myPulsesCell = {myPulse1, myPulse2, myPulse3, myPulse4};

% Pulse responses
myPulseResponse1a = zeros(nSamples, 1);
myPulseResponse1b = zeros(nSamples, 1);
myPulseResponse2a = zeros(nSamples, 1);
myPulseResponse2b = zeros(nSamples, 1);

cumNSamplesEachBin1 = cumsum(nSamplesEachBin1);
cumNSamplesEachBin2 = cumsum(nSamplesEachBin2);

indPulse1 = cumNSamplesEachBin1(1):cumNSamplesEachBin1(2);
indPulse2 = cumNSamplesEachBin2(1):cumNSamplesEachBin2(2);

indAfter1 = cumNSamplesEachBin1(2):nSamples;
indAfter2 = cumNSamplesEachBin2(2):nSamples;

shiftedTimeVec1 = siMs * ((1:length(indPulse1)) - 1);
myPulseResponse1a(indPulse1) = (amp1 + amp2)*(1 - exp(-(shiftedTimeVec1)/tau1));
myPulseResponse1b(indPulse1) = amp1*(1 - exp(-(shiftedTimeVec1)/tau1)) + ...
                                amp2*(1 - exp(-(shiftedTimeVec1)/tau2));
shiftedTimeVec2 = siMs * ((1:length(indPulse2)) - 1);
myPulseResponse2a(indPulse2) = (-amp1 + -amp2)*(1 - exp(-(shiftedTimeVec2)/tau1));
myPulseResponse2b(indPulse2) = -amp1*(1 - exp(-(shiftedTimeVec2)/tau1)) + ...
                                -amp2*(1 - exp(-(shiftedTimeVec2)/tau2));

responseAmp1a = (amp1 + amp2)*(1 - exp(-(shiftedTimeVec1(end))/tau1));
responseAmp1b = [amp1*(1 - exp(-shiftedTimeVec1(end)/tau1)), ...
                    amp2*(1 - exp(-shiftedTimeVec1(end)/tau2))];
responseAmp2a = (-amp1 + -amp2)*(1 - exp(-(shiftedTimeVec2(end))/tau1));
responseAmp2b = [-amp1*(1 - exp(-shiftedTimeVec2(end)/tau1)), ...
                    -amp2*(1 - exp(-shiftedTimeVec2(end)/tau2))];

shiftedTimeVec1 = siMs * ((1:length(indAfter1)) - 1);
myPulseResponse1a(indAfter1) = responseAmp1a*(exp(-(shiftedTimeVec1)/tau1));
myPulseResponse1b(indAfter1) = responseAmp1b(1)*(exp(-(shiftedTimeVec1)/tau1)) + ...
                                responseAmp1b(2)*(exp(-(shiftedTimeVec1)/tau2));
shiftedTimeVec2 = siMs * ((1:length(indAfter2)) - 1);
myPulseResponse2a(indAfter2) = responseAmp2a*(exp(-(shiftedTimeVec2)/tau1));
myPulseResponse2b(indAfter2) = responseAmp2b(1)*(exp(-(shiftedTimeVec2)/tau1)) + ...
                                responseAmp2b(2)*(exp(-(shiftedTimeVec2)/tau2));

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%{
OLD CODE:

%}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%