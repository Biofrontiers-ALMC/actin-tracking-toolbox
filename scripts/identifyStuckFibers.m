%Import data
% AD = ActinData;
% AD = importdata(AD, 'RQPeri_2spf_1.mat');
%AT = setFileMetadata(AT, 'meanDeltaT', 2);

% AD = analyze(AD);

WT = ActinData;
WT = importdata(WT, 'WTbeta_0,5spf_1.mat');

%Have to correct the timestamp metadata
WT = setFileMetadata(WT, 'meanDeltaT', 5);
WT = analyze(WT);

%Collect th
stuckFibers = [100, 4, 8, 81, 87, 89, 332, 310, 272, 247];
notStuckFibers = [440, 124, 315, 535, 156, 92];

stuckFibers_avgSpeed = mean([WT.Tracks(stuckFibers).AverageSpeed]);
notStuckFibers_avgSpeed = mean([WT.Tracks(notStuckFibers).AverageSpeed]);

%Plot histogram of the WT data

histogram([WT.Tracks.AverageSpeed], 100)