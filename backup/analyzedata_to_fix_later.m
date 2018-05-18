function [data] = analyzedata(fn)

[data, si] = abf2load(fn); % si is the sampling interval
