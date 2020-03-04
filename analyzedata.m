clearvars

%Example of importing the data to the analyzer object
AD = ActinData;
AD = importdata(AD, 'D:\Projects\2020Feb Leinwand Actin\data\trackdata_20200302.mat');

AD = setFileMetadata(AD, 'PxSize', 0.13, ...
    'DeltaT', 2);

AD = analyze(AD);

%Example of filtering
minFilamentLength = 3;
minNumFrames = 3;

tracks_pass = [AD.Tracks.NumFrames] > minNumFrames & [AD.Tracks.MeanFilamentLength] > minFilamentLength;

%Compute average instant speeds of all tracks that pass
meanSpeed = mean([AD.Tracks(tracks_pass).MeanInstantSpeed]);

%Histogram
figure;
histogram([AD.Tracks(tracks_pass).MeanInstantSpeed]);

%Example of plotting - track 3
track = getTrack(AD, 3);
tt = track.Frames * AD.FileMetadata.DeltaT;

%TODO: How many stuck

figure;
plot(tt, track.InstantSpeed)
xlabel('Time (s)');
ylabel('Instantaneous Speed (\mum/s)');


export(AD, 'test.csv');