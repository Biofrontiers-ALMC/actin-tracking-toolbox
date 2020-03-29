classdef ActinTracker
    
    properties
        
        %Segmentation options
        RescaleFactor = 4;  %Factor to increase image by - likely should be even
        BackgroundPercentile = 10;  %After background subtraction
        ThresholdFactor = 15;  %Fiber = background * threshold factor
        NormalizeIntensity = true;  %Should intensity be normalized
        
        ChannelToSegment = 1;
        
        %Relevant LAPLinker settings
        LinkScoreRange = [0 30];
        MaxTrackAge = 1;        
        TrackDivision = false;
        
    end
        
    methods
        
        function process(obj, file)
            %PROCESS  Set up processing
            %
            %  PROCESS(OBJ) will open a dialog box allowing you to select
            %  files to process.
            
            
            %  This method checks that files exist, then calls the static
            %  method processFile on each file. This allows the whole
            %  process to be carried out in a parfor loop if needed.
            
            C = metaclass(obj);
            P = C.Properties;
            for k = 1:length(P)
                if ~P{k}.Dependent
                    settings.(P{k}.Name) = obj.(P{k}.Name);                    
                end
            end
            
            ActinTracker.processFile(file, 'D:\Projects\2020Feb Leinwand Actin\data\test', settings)
            
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
        
        function processFile(file, outputDir, settings)
            %PROCESSFILE  Segment and track fibers in a file
            %
            %  ActinTracker.PROCESSFILE(FN, OUTDIR, SETTINGS) will process
            %  the file specified by the filename FN. SETTINGS must be a
            %  struct that contains the same fieldnames as the public
            %  properties of this object. OUTDIR should be a pathname to
            %  the output folder, which must already exist.
            %
            %  This is a static method.
            
            %TODO: Normalize intensity if set
            if settings.NormalizeIntensity
                disp('Noramlize intensity - to be done')
                %intProfile = generateIntensityProfile(filename, channel);
                
            end
            
            reader = BioformatsImage(file);

            %Set up the linker object
            linker = LAPLinker(settings);

            %Set up output file details
            [~, fn] = fileparts(file);
                        
            vid = VideoWriter(fullfile(outputDir, [fn, '.avi']));
            vid.FrameRate = 10;
            open(vid);
            
            for iT = 1:reader.sizeT
                
                %Read and pre-process the image
                I = getPlane(reader, 1, settings.ChannelToSegment, iT);
                
                %-- Pre-processing --$
                %Increasing the resolution of the image to help with
                %segmentation
                Iseg = imresize(I, settings.RescaleFactor, 'nearest');
                Iseg = imgaussfilt(Iseg, 1);
                
                bg = imerode(Iseg, strel('disk', 15));
                
                IbgSub = Iseg - bg;
                IbgSub(IbgSub == 0) = 1;
                
                %-- Segment ---%
                
                %Determine background value by smoothing the image a bit, then
                %taking the darkest pixels
                Ismooth = imgaussfilt(IbgSub, 0.5);
                
                %Assume bg is bottom 5th percentile
                bgVal = prctile(Ismooth(:), settings.BackgroundPercentile);
                clearvars('Ismooth');
                
                %Make a mask of the fibers
                mask = IbgSub > (bgVal * settings.ThresholdFactor);
                
                %Reduce the mask back to the original size
                mask = imresize(mask, 1/settings.RescaleFactor, 'nearest');
                
                %Skeletonize the fibers to get the center line
                mask = bwskel(mask);  %Remove small branches
                %mask = bwareaopen(mask, 5); %Clean up small objects
                
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
                        
                        %Measure the centroid and the length of the line
                        [yy, xx] = ind2sub(size(mask), data(idx).PixelIdxList);
                        lineData = ActinTracker.getLineData([xx, yy]);
                        
                        data(idx).Centroid = lineData.Centroid;
                        data(idx).Length = lineData.Length;
                        
                    end
                end
                
                %Track the fibers
                linker = assignToTrack(linker, iT, data);
                
                %Make a movie - TODO
                Iout = showoverlay(I, mask);
                for iTrack = 1:numel(linker.activeTrackIDs)
                    
                    track = getTrack(linker, linker.activeTrackIDs(iTrack));
                    Iout = insertText(Iout, track.Centroid(end, :), linker.activeTrackIDs(iTrack), ...
                        'BoxOpacity', 0, 'TextColor', 'yellow');
                    %TODO OPACITY AND BOX COLOR                    
                    
                    Iout = double(Iout);
                    Iout = Iout./max(Iout(:));
                end
                writeVideo(vid, Iout);
                clearvars Iout
            
            end
            
            close(vid)
            tracks = linker.tracks;
            save(fullfile(outputDir, [fn, '.mat']), 'tracks');
            %Save data - TODO METADATA
            
        end
        
        function lineData = getLineData(coords)
            %GETLINEDATA  Measures data of a line
            %
            %  S = GETLINEDATA(C) returns a struct S containing data about
            %  the line specified by the coordinate vector C. C must be an
            %  Nx2 list of coordinates where C(:, 1) is x and C(:, 2) is y.
            %  It is expected that the points in C is connected by no more
            %  than 1.4 pixels (i.e. the object mask is 8-connected).
            %
            %  S has the following fields:
            %     Centroid = Center coordinate of the line
            %     Length = Length of the line
            %     SortedCoords = Sorted coordinates of the line
            %     Excluded = Any coordinates that were excluded
            %
            %  The SortedCoords will be an Nx2 matrix specifying
            %  coordinates going from one end of the line to the other.
            %
            %  Excluded contains a list of coordinates that were excluded
            %  from the line. These are points that did not connect to the
            %  traced line, e.g. if there was a branch.
            %
            %  The algorithm works by first tracing the points to find an
            %  end point. The cumulative distance to the end point is then
            %  computed and used to order the points as well as getting the
            %  length of the line. The middle point is then the point that
            %  is closest to the length/2.
            %
            %  Example:
            %  %Assume after segmentation and skeletonization you have a
            %  %mask. Note the use of 'PixelList' instead of
            %  %'PixelIdxList'.
            %  data = regionprops(mask, 'PixelList');
            %
            %  lineData = ActinTracker.findLineCenter(data.PixelList);
            
            
            %Initialize variables            
            isSorted = [true; false(size(coords, 1) - 1, 1)];  %Flag if point has been sorted
            sortIdx = [1; zeros(size(coords, 1) - 1, 1)];  %To store sorted indices
            ptrLastIndex = 1;  %Position of sort index to add to
            
            %Find an end point by travelling in one direction from the
            %first pixel.
            minInd = 1;  %Start at the first coordinate
            while ~all(isSorted)
                
                %Find the next nearest pixel
                sqDistToPtr = sum((coords - coords(minInd, :)).^2, 2);
                sqDistToPtr(isSorted) = Inf;
                
                [minDist, minInd] = min(sqDistToPtr);
                
                if minDist <= 2
                    isSorted(minInd) = true;
                    
                    %Append to sorted indices
                    ptrLastIndex = ptrLastIndex + 1;
                    sortIdx(ptrLastIndex) = minInd;
                    
                else
                    break;
                end
            end

            %Shift the indices to the end of the array
            sortIdx = circshift(sortIdx, nnz(~isSorted));
            
            %Pointer will now count upwards so update the value
            ptrLastIndex = nnz(~isSorted) + 1;
            
            minInd = 1; %Reset the point back to the first coordinate
            while ~all(isSorted)
                
                %Find the next nearest pixel
                sqDistToPtr = sum((coords - coords(minInd, :)).^2, 2);
                sqDistToPtr(isSorted) = Inf;
                
                [minDist, minInd] = min(sqDistToPtr);
                
                if minDist <= 2
                    isSorted(minInd) = true;
                    
                    %Add to sorted indices going upwards
                    ptrLastIndex = ptrLastIndex - 1;
                    sortIdx(ptrLastIndex) = minInd;
                    
                else
                    break;
                end
                                
            end
            
            %Sort the array
            lineData.SortedCoords = coords(sortIdx, :);
            
            %Compute the line length
            distFromEnd = [0; cumsum(sqrt(sum((diff(lineData.SortedCoords)).^2, 2)))];
            lineData.Length = distFromEnd(end);
            
            %Find the center coordinate of the line
            [~, midPtLoc] = min(abs(distFromEnd - lineData.Length/2));
            lineData.Centroid = lineData.SortedCoords(midPtLoc, :);
            
            %Report any excluded data points
            lineData.Excluded = coords(~isSorted, :);
            
            
        end

    end
    
    
    
end














