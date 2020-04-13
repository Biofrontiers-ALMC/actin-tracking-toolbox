clearvars

%Example of importing the data to the analyzer object
AD = ActinData;
%AD = importdata(AD, 'D:\Projects\2020Feb Leinwand Actin\data\20200330_fibermetric\RQPeri_2spf_1.mat');
AD = importdata(AD, 'D:\Projects\2020Feb Leinwand Actin\data\20200330_fibermetric\WTbeta_0,5spf_1.mat');

%The file has the issue where timestamps are not read correctly
%AD = setFileMetadata(AD, 'meanDeltaT', 2);  %RQPeri
AD = setFileMetadata(AD, 'meanDeltaT', 0.5);  %WTbeta0

AD = analyze(AD);

export(AD, 'WTbeta_0,5spf_1.csv');
return

%%
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
track = getTrack(AD, 5);
tt = track.Frames * AD.FileMetadata.meanDeltaT;

figure;
plot(tt, track.InstantSpeed)
xlabel('Time (s)');
ylabel('Instantaneous Speed (\mum/s)');

%Export data to a CSV file
export(AD, 'test.csv');

%TODO: How many stuck
tracks_pass_id = find(tracks_pass);
numNotMoving = 0;
notMovingIds = [];

for iT =  1%:numel(tracks_pass)
    
    track = getTrack(AD, iT);
    
    %Compute the distance travelled between subsequent frames
    distTravelled = sqrt(sum((diff(track.Centroid)).^2, 2));
    
    %Determine if there was a period of time in which the fiber was "stuck"
    %or not moving much
    isNotMoving = nnz(distTravelled < 1);

    if isNotMoving > 5
        %Let's just say if the fiber doesn't move much in 5 frames, call it
        %not moving
        numNotMoving = numNotMoving + 1;
        notMovingIds(end + 1) = iT;        
    end
end

notMovingFrac = numNotMoving / numel(tracks_pass_id)


























