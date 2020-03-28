classdef ActinTracker
    
    properties
        
        %Segmentation options
        RescaleFactor = 3;  %Factor to increase image by
        BackgroundPercentile = 10;  %After background subtraction
        ThresholdFactor = 20;  %Fiber = background * threshold factor
        NormalizeIntensity = true;  %Should intensity be normalized
        
        ChannelToSegment = 'GFP';
        
    end
    
    
    methods
        
        function process(obj, varargin)
            %PROCESS  Set up processing
            %
            %  PROCESS(OBJ) will open a dialog box allowing you to select
            %  files to process.
            
            
            %  This method checks that files exist, then calls the static
            %  method processFile on each file. This allows the whole
            %  process to be carried out in a parfor loop if needed.
            
            
        end
        

        
    end
    
    methods (Static)
        
        function intProfile = generateIntensityProfile(filename, channel)
            %GENERATEINTENSITYPROFILE  Generates an intensity profile
            %
            %  P = GENERATEINTENSITYPROFILE(filename, channel) estimates an
            %  intensity profile P.
            %
            %  Examples:
            %             figure;
%             imshow(intProfile, [])
%             
%             Icorr = double(I) ./ intProfile;
%             imshow(Icorr, [])
%             

            
            %Open the file
            reader = BioformatsImage(filename);
            
            %Compute the mean intensity of all frames
            I = zeros(reader.height, reader.width);
            
            for iT = 1:reader.sizeT
                I = mean(cat(3, I, double(getPlane(reader, 1, channel, iT))), 3);
            end
            
            %Take a rough mask to identify objects
            mask = imbinarize(I, 'adaptive');
            
            %Fit measured data to a 2D Gaussian
            gauss2D =  fittype('A * exp(-((xx - B).^2 + (yy - C).^2) / (2*D.^2))',...
                'independent', {'xx', 'yy'});
                        
            %Generate the fitting coordinates
            [yyFit, xxFit] = find(mask);
            
            %Generate an initial guess
            [maxInt, maxIntInd] = max(I(mask));
            initGuess = [maxInt, xyFit(maxIntInd, 1), xyFit(maxIntInd, 2), 100];
            
            %Fit the detected pixel intensities
            fitObj = fit([xxFit, yyFit], I(mask), gauss2D, 'StartPoint', initGuess);
            
            %figure;
            %plot(fitObj, xyFit, intFit)
            
            %Generate the output intensity profile
            xx = 1:size(I, 2);
            yy = 1:size(I, 1);
            [xx, yy] = meshgrid(xx, yy);
            
            intProfile = fitObj.A * exp(-((xx - fitObj.B).^2 + (yy - fitObj.C).^2) / (2*fitObj.D.^2));
           
        end
        
    end
    
    methods (Static, Hidden)
        
        function processFile(obj, file, settings)
            %PROCESSFILE  Segment and track fibers in a file
            %
            %  ActinTracker.PROCESSFILE(FN, SETTINGS) will process the
            %  file specified by the filename FN. SETTINGS must be a struct
            %  that contains the same fieldnames as the public properties
            %  of this object.
            %
            %  This is a static method.
            
            %Normalize intensity if set
            if settings.NormalizeIntensity
                
                intProfile = generateIntensityProfile(filename, channel);
                
            end
            
            reader = BioformatsImage(file);
                        
            %Read and pre-process the image
            I = imread('fullframe.tiff');
            
            %Increasing the resolution of the image seems to help with
            %segmentation            
            Iseg = imresize(I, settings.RescaleFactor, 'nearest');
            Iseg = imgaussfilt(Iseg, 1);
            
            bg = imerode(Iseg, strel('disk', 15));
            
            IbgSub = Iseg - bg;
            IbgSub(IbgSub == 0) = 1;
            
            %TODO: Normalize the intensity of the image
            
            
            %-- Segment ---%
            
            %Determine background value by smoothing the image a bit, then
            %taking the darkest pixels
            Ismooth = imgaussfilt(Iseg, 0.5);
            
            %Assume bg is bottom 5th percentile
            bgVal = prctile(Ismooth(:), 10);
            clearvars('Ismooth');
            
            %Make a mask of the fibers
            mask = Iseg > (bgVal * 20);
            
            %Skeletonize the fibers to get the center line
            mask = bwskel(mask, 'MinBranchLength', 5);  %Remove small branches
            mask = bwareaopen(mask, 5); %Clean up small objects
            
            %Detect branching lines - objects with branches will be removed
            %in the final mask
            bpMask = bwmorph(mask, 'branchpoints');
            bpInd = find(bpMask);
            
            %--DEBUG
            %figure;
            %Iout = showoverlay(IbgSub, mask);
            %showoverlay(Iout, bpMask, 'Color', [1 0 0]);
            
            clearvars bpMask;
                        
            %-- Measure --
            
            %Find connected components
            cc = bwconncomp(mask, 8);
            
            %Declare structure to store data
            data = struct('PixelIdxList', {});
           
            for ii = 1:cc.NumObjects
                
                %Check if the object has a branch point - these are most
                %likely intersecting fibers                
                if any(ismember(bpInd, cc.PixelIdxList{ii}))
                    
                    %Remove branching objects from the mask
                    mask(cc.PixelIdxList{ii}) = false;
                else
                    
                    %Otherwise, compute and store data
                    idx = numel(data) + 1;
                    data(idx).PixelIdxList = cc.PixelIdxList{ii};
                    
                    %Procedure to compute the actual length of the object
                    %and the location of the mid point:
                    %
                    %  * First find an pole of the object
                    %  * Order the pixels by distance to the pole
                    %  * Length is then the sum of consecutive sequences
                    %  * The midpoint is the coordinate that is closest to
                    %    half the length of the object
                    
                    %For each fiber, convert the list of pixel indices to
                    %coordinates
                    [ycoord, xcoord] = ind2sub(size(Iseg), cc.PixelIdxList{ii});
                    
                    %Compute the pairwise distance of each point with every
                    %other point
                    dist = pdist([xcoord, ycoord]);
                    dist = squareform(dist);
                    
                    %Find the two points that are furthest apart from each
                    %other. ACTUALLY THIS WOULDN'T WORK FOR LINES THAT
                    %CURVE ON THEMSELVES
                    
                    %Start with a random point, then travel one way, and
                    %travel the other. Say that for a line to be connected,
                    %the points have to be within <2 pixels
                    
                    
                    
                    
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
            showoverlay(IbgSub, bwperim(mask));
            hold on
            centroid = cat(1,  cc.midPtCoord{:});
            plot(centroid(:, 2), centroid(:, 1), 'rx');
            
            endPt = cat(1,  cc.endPtLoc{:});
            plot(endPt(:, 2), endPt(:, 1), 'yo');
            
            hold off
            




        end

    end
    
    
    
end