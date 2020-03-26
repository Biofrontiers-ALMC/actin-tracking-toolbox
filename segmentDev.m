%Read and pre-process the image
I = imread('fullframe.tiff');

Irescale = imresize(I, 3, 'nearest');
Irescale = imgaussfilt(Irescale, 1);

bg = imerode(Irescale, strel('disk', 15));

Isub = Irescale - bg;
Isub(Isub == 0) = 1;
%Try to normalize intensity of image

%Take a rough mask to identify objects
mask = imbinarize(Isub, 'adaptive', 'Sensitivity', 0.001);
mask = imopen(mask, strel('disk', 4));
mask = bwareaopen(mask, 10);
showoverlay(Isub, mask)

%Measure intensity of detected objects
rp = regionprops(mask, Isub, 'Centroid', 'MeanIntensity');

%Fit measured data to a 2D Gaussian
gauss2D =  fittype('A * exp(-((xx - B).^2 + (yy - C).^2) / (2*D.^2))',...
    'independent', {'xx', 'yy'});

xyFit = cat(1, rp.Centroid);
intFit = cat(1, rp.MeanIntensity);

[maxInt, maxIntInd] = max(intFit);

initGuess = [maxInt, xyFit(maxIntInd, 1), xyFit(maxIntInd, 2), 100];

fitObj = fit(xyFit, intFit, gauss2D, 'StartPoint', initGuess);

figure;
plot(fitObj, xyFit, intFit)

%Reconstruct the intensity profile
xx = 1:size(Isub, 2);
yy = 1:size(Isub, 1);
[xx, yy] = meshgrid(xx, yy);

reconInt = fitObj.A * exp(-((xx - fitObj.B).^2 + (yy - fitObj.C).^2) / (2*fitObj.D.^2));

figure;
imshow(reconInt, [])

Icorr = double(Isub) ./ reconInt;
imshow(Icorr, [])

%% Segment
Icorr = imgaussfilt(Icorr, 0.5);

%Assume bg is bottom 5th percentile
bgVal = prctile(Icorr(:), 10);

mask = Icorr > (bgVal * 20);

%Skeletonize - seems to work better if we keep the image larger
mask_skel = bwskel(mask, 'MinBranchLength', 5);
mask_skel = bwareaopen(mask_skel, 5);

%Detect branching lines
bp = bwmorph(mask_skel, 'branchpoints');
bpInd = find(bp);

figure;
Iout = showoverlay(Isub, mask_skel);
showoverlay(Iout, bp, 'Color', [1 0 0]);

%% Estimate arc length by using straight-line segments
%An easier/likely good estimate is simply the number of pixels

cc = bwconncomp(mask_skel, 8);

%Find middle point by finding segment that is closest to middle?
centroid = zeros(cc.NumObjects, 2);
removeObj = false(1, cc.NumObjects);

for ii = 1:cc.NumObjects
    
    %Check if the object has a branch point. If yes, then remove it
    if any(ismember(bpInd, cc.PixelIdxList{ii}))
        removeObj(ii) = true;
        
        %Edit the mask
        mask_skel(cc.PixelIdxList{ii}) = false;
    else
        
        %--Find the mid point--
        
        %Convert pixel indices to row, column coordinates
        [II, JJ] = ind2sub(size(Irescale), cc.PixelIdxList{ii});
        cc.PixelCoords{ii} = [II, JJ];
        
        %Find points that are furthest apart from each other
        dist = pdist(cc.PixelCoords{ii});
        dist = squareform(dist);
        
        [~, maxDistInd] = max(dist, [], 'all', 'linear');
        [subI, subJ] = ind2sub(size(dist), maxDistInd);
        
        %Choose one end
        cc.endPtLoc{ii} = cc.PixelCoords{ii}(subJ, :);
        
        %Sort the points from distance to end point by connecting the
        %points
        sortedPts = zeros(size(cc.PixelCoords{ii}, 1), 2);
        sortedPts(1, :) = cc.endPtLoc{ii};  %Initialize with end point location
        
        unsortedPts = cc.PixelCoords{ii};
        
        ctr = 1;
        while ctr <= size(sortedPts, 1)
        
            %Find next nearest-neighbor
            distToEndPt = sum((unsortedPts - cc.endPtLoc{ii}).^2, 2);
            [~, minDistToEnd] = min(distToEndPt);
            
            sortedPts(ctr, :) = unsortedPts(minDistToEnd, :);
            unsortedPts(minDistToEnd, :) = [];
        
            ctr = ctr + 1;
        end
                
        cc.sortedPts{ii} = sortedPts;
        
        %Calculate length of line using the sorted coordinates
        cc.totalLen{ii} = sum(sqrt((diff(sortedPts(:, 1))).^2 + (diff(sortedPts(:, 2))).^2));
        
        cc.midPtCoord{ii} = cc.sortedPts{ii}(round(size(cc.sortedPts{ii}, 1)/2), :);
        
    end
end


%Plot estimated centroids
figure;
showoverlay(Isub, bwperim(mask_skel));
hold on
centroid = cat(1,  cc.midPtCoord{:});
plot(centroid(:, 2), centroid(:, 1), 'rx');

endPt = cat(1,  cc.endPtLoc{:});
plot(endPt(:, 2), endPt(:, 1), 'yo');

hold off













