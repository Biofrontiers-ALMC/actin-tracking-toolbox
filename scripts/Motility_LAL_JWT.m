% %% Motility analysis
% % This script outlines how to analyze motility videos
% 
% %% Process Movie - Find Centroids
% AT = ActinTracker;
% AT.ChannelToSegment = 'Texas Red';
% AT.LinkScoreRange = [0 30];
% process(AT)
% 
% %% Process All Movies In a Folder
% 
% moviesLocation = '/Volumes/LINDSEY LEE/20210727_motility'; %copy file to this line
% movies = dir(moviesLocation);
% 
% for i = 3:length(movies)
%     
%     AT = ActinTracker;
%     AT.ChannelToSegment = 'Texas Red';
%     AT.LinkScoreRange = [0 30];
%     process(AT,[moviesLocation '/' movies(i).name],[moviesLocation '/Results'])
%     
% end
%% Review Video
% Review .avi movie file created by process(AT). Identify example stuck and
% non-stuck fibers.

stuckFibersID = [91, 77, 311, 3, 84];
movingFibersID = [40, 57, 34, 95, 101];

%% Data analysis
% Only run after stuck and non-stuck fiber examples have been defined above

dataPath = ''; % if this is empty, the script will open a dialog box
if ~isempty(dataPath)
    data = ActinData_LAL2_JWT;
    data = importdata(data,dataPath);
else
    data = ActinData_LAL2_JWT;
    data = importdata(data);
end

% Set proper timestamp by changing "T" to frame rate

%data = setFileMetadata(data, 'meanDeltaT', 1); % The third input to this function is the average time between frames in seconds

%HACK to grab and use timestamps
nd2Filepath = 'D:\Projects\2020Feb Leinwand Actin\data\WT_1_3fps_1001.nd2';
nd2 = ND2reader(nd2Filepath);

%Get the timestamps
data.timestamps = getTimestamps(nd2, 1, 1);

%Made an additional filter for filament length and number of frames
minFilamentLength = 1; %number = anything less than will be filtered out
numFrames = 15; %filters out <= number provided, use 15 to select filaments in at least half the frames
data = analyze(data,minFilamentLength,numFrames);

% Calculate average stuck and non-stuck velocities
stuckFibers_avgSpeed = mean([data.Tracks(stuckFibersID).AverageSpeed], 'all', 'omitnan');
notStuckFibers_avgSpeed = mean([data.Tracks(movingFibersID).AverageSpeed], 'all', 'omitnan');

data.maxStuckSpeed = max([data.Tracks(stuckFibersID).AverageSpeed], [], 'all', 'omitnan')+.001; % The +.001 is to get around a greater than or equal comparision.

data = analyze(data,minFilamentLength,numFrames);
% Print results - represents just filaments selected above
fprintf(['Stuck fibers average speed is: ' num2str(stuckFibers_avgSpeed) ' um/s\n'])
fprintf(['Non-stuck fibers average speed is: ' num2str(notStuckFibers_avgSpeed) ' um/s\n'])


%% Export Data
newFormat = true; % true = Lindsey's format, false = original format + stuck
[FILEPATH,NAME,EXT] = fileparts (data.FileMetadata.filename);
resultsLocation = [FILEPATH '\Results'];
export(data,[resultsLocation '\' NAME '.csv'],newFormat)

%% Show stuck and non stuck frames
frame = 1;
showlabels(data, frame, 79) % If last input is not zero, then the corresponding filament label will be green

%% Check single filament
track = getTrack(data,8); % second input is filament number

tt = track.Frames * data.FileMetadata.meanDeltaT;

% Plots single filament's velocity (um/s) vs. time (sec)
figure;
plot(tt, track.InstantaneousSpeed)
xlabel('Time (s)');
ylabel('Instantaneous Speed (\mum/s)');

%% More plots
% All filaments
for i = 1:length(data.Tracks)
    filamentLength(i) = data.Tracks(i).AverageFilamentLength;
    averageSpeed(i) = data.Tracks(i).AverageSpeed;
end

figure
hold on
scatter(filamentLength,averageSpeed)
xlabel('Average Filament Length (um)')
ylabel('Average Velocity (um/s)')
title('All Filament Length vs. Velocity')
grid on
xlim([0 15]) % [lowerLenghtBound upperLengthBound
ylim([0 7]) % [lowerVelocityBound upperVelocityBound]
hold off

% Stuck filaments
count = 1;
filamentLength = [];
averageSpeed = [];
for i = 1:length(data.Tracks)
    if data.Tracks(i).isStuck
    filamentLength(i) = data.Tracks(i).AverageFilamentLength;
    averageSpeed(i) = data.Tracks(i).AverageSpeed;
    count = count + 1;
    end
end

figure
hold on
scatter(filamentLength,averageSpeed)
xlabel('Average Filament Length (um)')
ylabel('Average Velocity (um/s)')
title('Stuck Filament Length vs. Velocity')
grid on
xlim([0 15]) % [lowerLenghtBound upperLengthBound
ylim([0 .5]) % [lowerVelocityBound upperVelocityBound]
hold off

% Not Stuck filaments
count = 1;
filamentLength = [];
averageSpeed = [];
for i = 1:length(data.Tracks)
    if ~data.Tracks(i).isStuck
    filamentLength(i) = data.Tracks(i).AverageFilamentLength;
    averageSpeed(i) = data.Tracks(i).AverageSpeed;
    count = count + 1;
    end
end

figure
hold on
scatter(filamentLength,averageSpeed)
xlabel('Average Filament Length (um)')
ylabel('Average Velocity (um/s)')
title('Non-Stuck Filament Length vs. Velocity')
grid on
xlim([0 15]) % [lowerLenghtBound upperLengthBound
ylim([0 7]) % [lowerVelocityBound upperVelocityBound]
hold off


%% Histogram

% All filaments histogram
avgSpeed = zeros([1 length(data.Tracks)]);
for i = 1:length(data.Tracks)
    avgSpeed(i) = data.Tracks(i).AverageSpeed;
end
figure
hold on
grid on
title('Average Speed Histogram of All Filaments (um/s)')
xlabel('Speed (um/s)')
ylabel('Count')
histogram(avgSpeed,30','BinLimits',[0 1]);

% Stuck filaments histogram
count = 1;
for i = 1:length(data.Tracks)
    if data.Tracks(i).isStuck
        avgSpeedStuck(count) = data.Tracks(i).AverageSpeed;
        count = count + 1;
    end
end
figure
hold on
grid on
title('Average Speed Histogram of All Stuck Filaments (um/s)')
xlabel('Speed (um/s)')
ylabel('Count')
histogram(avgSpeedStuck,30','BinLimits',[0 1]);

% Not Stuck filaments histogram
count = 1;
for i = 1:length(data.Tracks)
    if ~data.Tracks(i).isStuck
        avgSpeedNotStuck(count) = data.Tracks(i).AverageSpeed;
        count = count + 1;
    end
end
figure
hold on
grid on
title('Average Speed Histogram of All Not Stuck Filaments (um/s)')
xlabel('Speed (um/s)')
ylabel('Count')
histogram(avgSpeedNotStuck,30','BinLimits',[0 1]);

