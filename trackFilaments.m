%Instatnenous speed
%Length
%pos from leading edge
clearvars

bfr = BioformatsImage('D:\Projects\2020Feb Leinwand Mitochondria\data\RQPeri_2spf_1.nd2');
Linker = LAPLinker;
Linker.LinkScoreRange = [0, 30];

vid = VideoWriter('test.avi');
open(vid)

for iT = 1:bfr.sizeT

    I = double(getPlane(bfr, 1, 1, iT));
    
    %%
    %Segmentation
    
    %Background subtract - there is bleedthrough of intensity
    bgImg = imerode(I, strel('disk', 10));
    
    Isub = I - bgImg;
    Isub(Isub < 0) = 0;
    Isub = Isub ./ max(Isub(:));
    
    B = fibermetric(Isub, [3, 4, 5]);
    
    fiberMask = B > 0.6;
    fiberMask = bwareaopen(fiberMask, 7);
    fiberMask = imdilate(fiberMask, 1);
    
%     fiberMask = bwmorph(fiberMask, 'skel', Inf);
%     fiberMask_ep = bwmorph(fiberMask, 'endpoints');
%     fiberMask_bp = bwmorph(fiberMask, 'branchpoints');
%     
%     fiberMask = fiberMask - fiberMask_bp - fiberMask_ep;
%     fiberMask = fiberMask > 0.5;
%     
%     fiberMask = bwareaopen(fiberMask, 3);
    
    %showoverlay(Isub, bwperim(fiberMask > 0.5), 'Opacity', 30);
    
        
    data = regionprops(fiberMask, 'Centroid', 'PixelIdxList', 'MajorAxisLength');
    
%     showoverlay(Isub, bwperim(fiberMask));
%     centroid = cat(1, data.Centroid);
%     hold on
%     plot(centroid(:, 1), centroid(:, 2), 'o')
%     hold off

    
    Linker = assignToTrack(Linker, iT, data);

    %Generate output video
    Iout = showoverlay(Isub, bwperim(fiberMask));
    
    for iL = 1:numel(Linker.activeTrackIDs)
        
        currTrack = getTrack(Linker, Linker.activeTrackIDs(iL));
        
        Iout = insertText(Iout, currTrack.Centroid(end, :), Linker.activeTrackIDs(iL), ...
            'BoxOpacity', 0, 'textcolor', 'w');        
    end
    
    Iout = Iout./ max(Iout(:));
    
    writeVideo(vid, Iout);
    
end

close(vid)

fnOut = fileparts(bfr.filename);
save([fnOut, '.mat'], 'Linker')

fh = figure;
Iout = getPlane(bfr, 1, 1, iT);
imshow(Iout)
hold on
for ii = 1:Linker.NumTracks
    
    track = getTrack(Linker, ii);
    plot(track.Centroid(:, 1), track.Centroid(:, 2));
    
end
hold off

saveas(fh, [fnOut, '.png']);