clearvars
load tracks.mat

pxSize = 0.13;
tUnits = 2;

%Compute instantaneous speed - histogram?
%length
for ii = 181%:numel(TrackArrayData.Tracks)
    
    %Concatenate the Track Position data
    posData = cat(1, TrackArrayData.Tracks(ii).Data.Centroid{:});
    
    %Displacement
    distTravelled = [0; sqrt(sum((diff(posData, 1)).^2, 2))];
    
    instantSpeed = (distTravelled * pxSize)/tUnits;
    
end

tt = (1:numel(instantSpeed)) * tUnits;
plot(tt, instantSpeed)
xlabel('Time (s)');
ylabel('Instantaneous Speed (\mum/s)');
