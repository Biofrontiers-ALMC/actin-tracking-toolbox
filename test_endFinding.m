clearvars
clc

mask = rgb2gray(imread('curvedLine.png'));
mask = mask > 0;

data = regionprops(mask, 'PixelIdxList');

[yy, xx] = ind2sub(size(mask), data.PixelIdxList);

%Start at the first point
unsortedPts = [xx, yy];

%Pick a point - we'll just start at the beginning
sortedPts = unsortedPts(1, :);
unsortedPts(1, :) = [];

%Go one way
while ~isempty(unsortedPts)
    diffXY = sortedPts(end, :) - unsortedPts;
    dist = sum(diffXY.^2, 2);
    
    [minDist, minInd] = min(dist);
    
    if minDist <= 2
        sortedPts = [sortedPts; unsortedPts(minInd, :)];
        unsortedPts(minInd, :) = [];
    else
        break
    end
end

%Go the other
while ~isempty(unsortedPts)
    
    diffXY = sortedPts(1, :) - unsortedPts;
    dist = sum(diffXY.^2, 2);
    
    [minDist, minInd] = min(dist);
    
    if minDist <= 2
        sortedPts = [unsortedPts(minInd, :); sortedPts];
        unsortedPts(minInd, :) = [];
    else
        break
    end   
end

%Find the mid-point
distFromEnd = [0; cumsum(sqrt(sum((diff(sortedPts)).^2, 2)))];

lineLength = distFromEnd(end);

[~, midPtLoc] = min(abs(distFromEnd - lineLength/2));
midPtCoord = sortedPts(midPtLoc, :);

imshow(mask)
hold on
plot(xx(1), yy(1), 'x')
plot(sortedPts(end, 1), sortedPts(end, 2), 'ro')
plot(sortedPts(1, 1), sortedPts(1, 2), 'ro')
plot(midPtCoord(1), midPtCoord(2), 'ks');
hold off