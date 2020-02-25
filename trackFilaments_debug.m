%Instatnenous speed
%Length
%pos from leading edge
clearvars
bfr = BioformatsImage('D:\Projects\2020Feb Leinwand Mitochondria\data\RQPeri_2spf_1.nd2');
Linker = LAPLinker;
Linker.LinkScoreRange = [0, 15];

for iT = 1:4

    I = double(getPlane(bfr, 1, 1, iT));
    
    %%
    %Segmentation
    
    %Background subtract - there is bleedthrough of intensity
    bgImg = imopen(I, strel('disk', 8));
    
    Isub = I - bgImg;
    Isub(Isub < 0) = 0;
    Isub = Isub ./ max(Isub(:));
    
    B = fibermetric(Isub, [4, 6, 10]);
    
    fiberMask = B > 0.7;
    fiberMask = bwareaopen(fiberMask, 8);
    fiberMask = imdilate(fiberMask, 1);
    
    fiberMask = bwskel(fiberMask);
%     showoverlay(Isub, bwperim(fiberMask));
    
    data = regionprops(fiberMask, 'Centroid', 'PixelIdxList');
    
%     showoverlay(Isub, bwperim(fiberMask));
%     centroid = cat(1, data.Centroid);
%     hold on
%     plot(centroid(:, 1), centroid(:, 2), 'o')
%     hold off

    if iT == 4
        keyboard
    end
    Linker = assignToTrack(Linker, iT, data);

    disp(Linker.NumTracks);
    disp(numel(data));
    
    %Generate output image
    Iout = showoverlay(Isub, bwperim(fiberMask));
    
    for iL = Linker.activeTrackIDs
        
        currTrack = getTrack(Linker, iL);
        
        Iout = insertText(Iout, currTrack.Centroid(end, :), iL);        
    end
    
    Iout = Iout./ max(Iout(:));
    
    storeData(iT).Centroid = cat(1, data.Centroid);
    
end

%%

figure(1)
showoverlay(Isub, bwperim(fiberMask));
hold on
plot(storeData(1).Centroid(:, 1), storeData(1).Centroid(:, 2), 'or');
plot(storeData(2).Centroid(:, 1), storeData(2).Centroid(:, 2), 'ob');
plot(storeData(3).Centroid(:, 1), storeData(3).Centroid(:, 2), 'og');
plot(storeData(4).Centroid(:, 1), storeData(4).Centroid(:, 2), 'oy');
for ii = 1:Linker.NumTracks
    
    track = getTrack(Linker, ii);
    plot(track.Centroid(:, 1), track.Centroid(:, 2));
    
end
hold off
% 
% 
% figure(2);
% imshow(Iout)
