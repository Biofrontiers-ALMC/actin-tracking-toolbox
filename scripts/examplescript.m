clearvars
clc

%Create an ActinTracker object
AT = ActinTracker;

%Process the file
process(AT, ...
    'C:\Users\Jian Tay\Dropbox\Actin_Tracking_PSU\Leica export Tiff\Project001_001_ch00.tif', ...
    'D:\Projects\Toolboxes\actin-tracking-toolbox\data\Jinghua\results')

%%

%After processing is complete, create an ActinData object to analyze
AD = ActinData;

%Import the resulting data
AD = importdata(AD, 'D:\Projects\Toolboxes\actin-tracking-toolbox\data\Jinghua\results\Project001_001_ch00.mat');

%If data is from a TIF-file, it will be missing important metadata. You can
%add these in manually. pxSize is the image pixel resolution (i.e., how
%many microns is each pixel in the image). meanDeltaT is the time between
%frames. Right now, we assume that this is constant.
AD = setFileMetadata(AD, 'pxSize', 2.119e-4, 'meanDeltaT', 1.19);

%Analyze the resulting data
AD = analyze(AD);

%Export the data to CSV for further analysis as required or we can develop
%a MATLAB script to work with the ActinData object.
export(AD, 'test.csv')


