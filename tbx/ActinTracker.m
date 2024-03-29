classdef ActinTracker
    
    properties
        
        %Segmentation options
        ChannelToSegment = 1;
        ThresholdFactor = 12;  %Fiber = background * threshold factor
        MinBranchLength = 3;  %Minimum length of branches
        
        %Relevant LAPLinker settings
        LinkScoreRange = [0 30];
        MaxTrackAge = 0;        
        TrackDivision = false;
        
    end
        
    methods
        
        function process(obj, files, outputDir)
            %PROCESS  Start processing files
            %
            %  PROCESS(OBJ) will open a dialog box allowing you to select
            %  ND2 files to process, as well as the output directory. The
            %  fibers in the selected files will be segmented and tracked.
            %
            %  The code outputs a video (*.avi) file to allow you to check
            %  segmentation and tracking. The code will also generate a
            %  MAT-file (*.mat) that contains the tracked the data. This
            %  MAT-file can be opened with the ActinAnalyzer object for
            %  further analysis. Both files will have the same name as the
            %  original image file.
            %
            %  During the segmentation step, each frame is resized by the
            %  RescaleFactor. This seems to improve detection of the
            %  fibers. The background of the image is then estimated and
            %  subtracted before the fibers are identified by intensity
            %  thresholding to generate a binary mask.
            %
            %  The detected fiber mask is then resized to its original
            %  size, then skeletonized to get the center line of each
            %  object. The length and center point of each line is then
            %  measured. Any lines with branches (usually indicating
            %  intersecting fibers) are removed.
            %
            %  The fibers are then tracked using a nearest-neighbor
            %  algorithm based on the center point.
            %
            %  There are several settings that you can adjust to improve
            %  segmentation and tracking:
            %
            %  Segmentation
            %     ChannelToSegment - Channel to segment
            %     RescaleFactor - Factor to resize images (must be even)
            %     BackgroundPercentile - Level of background
            %     ThresholdFactor - Theshold level above background
            %  
            %  Tracking
            %     LinkScoreRange - Range of pixels to link objects
            %     MaxTrackAge - Whether missing fibers should be retracked
            %     TrackDivision - Whether division should be tracked
            %                     (should be false)
            %
            %  PROCESS(OBJ, FILES, OUTPUTDIR) can be used to specify
            %  programatically the files and output directory. If more than
            %  one file is to be processed, FILES should be a cell array.
            %  Files is different directories can be specified by providing
            %  their full path. OUTPUTDIR should be the path to the output
            %  directory.
            
            
            %  This method checks that files exist, then calls the static
            %  method processFile on each file. This allows the whole
            %  process to be carried out in a parfor loop if needed.

            %Validate inputs
            if nargin == 1
                
                %Get files to process
                [files, fpath] = uigetfile({'*.nd2', 'ND2 files (*.nd2)'; ...
                                            '*.tif; *.tiff', 'TIFF files (*.tif; *.tiff)'; ...
                                            '*.*', 'All files (*.*)'}, ...
                                            'Select file(s) to process', ...
                                            'MultiSelect', 'on');
                %Cancel button pressed
                if isequal(files, 0)
                    return;                    
                end
                
                %Convert to cell if only a single file was selected
                if ~iscell(files)
                    files = {files};                    
                end
                
                %Append the full path name to each file
                for iF = 1:numel(files)
                    files{iF} = fullfile(fpath, files{iF});                    
                end

                %Get output directory
                outputDir = uigetdir(fpath, 'Select output directory');
                
                %Cancel button pressed
                if isequal(outputDir, 0)
                    return;                    
                end
                
            elseif nargin == 3
                
                if ~iscell(files)
                    files = {files};
                end
                
                %Check if files exist
                for ii = 1:numel(files)
                    if ~exist(files{ii}, 'file')
                        error('ActinTracker:process:FileNotFound', ...
                            'File %s was not found.', files{ii});
                    end
                end
                
                %Check if output directory exists. Otherwise, create it.
                if ~exist(outputDir, 'dir')
                    mkdir(outputDir);                    
                end
                
            else
                
                error('ActinTracker:process:InvalidNumberOfArguments', ...
                    'Invalid number of input arguments. Expected one or three.')
                
            end
            
            %Pack the settings into a struct            
            C = metaclass(obj);
            P = C.Properties;
            for k = 1:length(P)
                if ~P{k}.Dependent
                    settings.(P{k}.Name) = obj.(P{k}.Name);                    
                end
            end
            
            %Export settings
            exportSettings(obj, fullfile(outputDir, 'settings.txt'));
            
            %Start processing
            for iF = 1:numel(files)
                %Print a progress statement
                fprintf('%s: Started %s\n', datestr(now), files{iF});
                
                %Process file but catch errors that occur
                try
                    
                    ActinTracker.processFile(files{iF}, outputDir, settings);
                    
                catch ME
                    
                    %Print an error statement
                    fprintf('%s: Error processing %s\n', datestr(now), files{iF});
                    
                    msgText = getReport(ME);
                    fprintf('%s\n', msgText);
                
                end
                fprintf('%s: Completed %s\n', datestr(now), files{iF});
            end
            
        end
        
        function varargout = testsegment(obj, file, frame)
            %TESTSEGMENT  Test the segmentation settings
            %
            %  TESTSEGMENT(OBJ, FILE, FRAME) runs the segmentation routine
            %  on the file and frame specified.
            %
            
            reader = BioformatsImage(file);
            
            I = getPlane(reader, 1, obj.ChannelToSegment, frame);
            
            %Pack the settings into a struct
            C = metaclass(obj);
            P = C.Properties;
            for k = 1:length(P)
                if ~P{k}.Dependent
                    settings.(P{k}.Name) = obj.(P{k}.Name);
                end
            end
           
            [mask, bpInd, Isub] = ActinTracker.segment(I, settings);
            
            %Convert branch points indices to mask
            bpMask = false(size(mask));
            bpMask(bpInd) = true;
            
            Iout = ActinTracker.showoverlay(Isub, mask);
            ActinTracker.showoverlay(Iout, bpMask, 'color', [1 0 0]);
            
            if nargout >= 1
                varargout{1} = mask;
            end
            if nargout >= 2 
                varargout{2} = bpInd;
            end

            
        end
        
        function obj = importSettings(obj, varargin)
            %IMPORTSETTINGS  Import settings from file
            %
            %  IMPORTSETTINGS(OBJ, FILENAME) will load the options from the
            %  file specified.
            %
            %  IMPORTSETTINGS(OBJ) will open a dialog box for the user to
            %  select the option file.
            
            %Get the options file
            if isempty(varargin)
                [fname, fpath] = uigetfile({'*.txt','Text file (*.txt)';...
                    '*.*','All files (*.*)'},...
                    'Select settings file');
                
                if fname == 0
                    return;
                end
                
                optionsFile = fullfile(fpath,fname);
                
            elseif numel(varargin) == 1
                optionsFile = varargin{1};
                
            else
                error('ActinTracker:importSettings:TooManyInputs', 'Too many input arguments.');
                
            end
            
            fid = fopen(optionsFile,'r');
            
            if isequal(fid,-1)
                error('ActinTracker:importSettings:ErrorReadingFile',...
                    'Could not open file %s for reading.',optionsFile);
            end
            
            ctrLine = 0;
            while ~feof(fid)
                currLine = strtrim(fgetl(fid));
                ctrLine = ctrLine + 1;
                
                if isempty(currLine) || strcmpi(currLine(1),'%') || strcmpi(currLine(1),'#')
                    %Empty lines should be skipped
                    %Lines starting with '%' or '#' are comments, so ignore
                    %those
                    
                else
                    
                    parsedLine = strsplit(currLine,'=');
                    
                    %Check for errors in the options file
                    if numel(parsedLine) < 2
                        error('ActinTracker:importSettings:ErrorReadingOption',...
                            'Error reading <a href="matlab:opentoline(''%s'', %d)">file %s (line %d)</a>',...
                            optionsFile, ctrLine, optionsFile, ctrLine);
                    elseif isempty(parsedLine{2})
                        error('ActinTracker:importSettings:ErrorReadingOption',...
                            'Missing value in <a href="matlab:opentoline(''%s'', %d)">file %s (line %d)</a>',...
                            optionsFile, ctrLine, optionsFile, ctrLine);
                    end
                    
                    %Get parameter name (removing spaces)
                    parameterName = strtrim(parsedLine{1});
                    
                    %Get value name (removing spaces)
                    value = strtrim(parsedLine{2});
                    
                    if isempty(value)
                        %If value is empty, just use the default
                    else
                        obj.(parameterName) = eval(value);
                    end
                    
                end
                
            end
            
            fclose(fid);
        end
        
        function exportSettings(obj, exportFilename)
            %EXPORTOPTIONS  Export tracking options to a file
            %
            %  L.EXPORTOPTIONS(filename) will write the currently set
            %  options to the file specified. The options are written in
            %  plaintext, no matter what the extension of the file is.
            %
            %  L.EXPORTOPTIONS if the filename is not provided, a dialog
            %  box will pop-up asking the user to select a location to save
            %  the file.
            
            if ~exist('exportFilename','var')
                
                [filename, pathname] = uiputfile({'*.txt','Text file (*.txt)'},...
                    'Select output file location');
                
                exportFilename = fullfile(pathname,filename);
                
            end
            
            fid = fopen(exportFilename,'w');
            
            if fid == -1
                error('ActinTracker:exportSettings:CouldNotOpenFile',...
                    'Could not open file %s to write.', exportFilename)
            end
            
            %First, write the date
            fprintf(fid,'%%%s \r\n\r\n', datestr(now));
            
            propertyList = properties(obj);
            
            %Write output data depending on the datatype of the value
            for ii = 1:numel(propertyList)
                
                if ischar(obj.(propertyList{ii}))
                    fprintf(fid,'%s = ''%s'' \r\n',propertyList{ii}, ...
                        obj.(propertyList{ii}));
                    
                elseif isnumeric(obj.(propertyList{ii}))
                    fprintf(fid,'%s = %s \r\n',propertyList{ii}, ...
                        mat2str(obj.(propertyList{ii})));
                    
                elseif islogical(obj.(propertyList{ii}))
                    
                    if obj.(propertyList{ii})
                        fprintf(fid,'%s = true \r\n',propertyList{ii});
                    else
                        fprintf(fid,'%s = false \r\n',propertyList{ii});
                    end
                    
                end
                
            end
            
            fclose(fid);
            
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
            
%             %TODO: Normalize intensity if set
%             if settings.NormalizeIntensity
%                 disp('Noramlize intensity - to be done')
%                 %intProfile = generateIntensityProfile(filename, channel);
%                 
%             end
%             
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
                
                [mask, bpInd, Isub] = ActinTracker.segment(I, settings);
                 
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
                        
                        %Account for the resizing of the mask
                        data(idx).LineCoords = [xx, yy];% ./ settings.RescaleFactor;
                        
                        lineData = ActinTracker.getLineData(data(idx).LineCoords);
                        
                        data(idx).Centroid = lineData.Centroid;
                        data(idx).Length = lineData.Length;
                        
                    end
                end
                
                %Track the fibers
                linker = assignToTrack(linker, iT, data);
                
                %Make a movie
                Iout = double(Isub);
                Iout = Iout ./ max(Iout(:));
                for iTrack = 1:numel(linker.activeTrackIDs)
                    
                    track = getTrack(linker, linker.activeTrackIDs(iTrack));
                    
                    %Mark the centroid of each line
                    Iout = insertShape(Iout, 'FilledCircle', [track.Centroid(end, :), 3], 'Color', 'red');
                    
                    Iout = insertText(Iout, track.Centroid(end, :), linker.activeTrackIDs(iTrack), ...
                        'BoxOpacity', 0, 'TextColor', 'yellow');

                    Iout = Iout./max(Iout(:));
                end
                writeVideo(vid, Iout);
                clearvars Iout
            
            end
            
            close(vid)
            
            %Set file metadata
            linker = updateMetadata(linker, ...
                'filename', file,...
                'imgSize', [reader.height, reader.width]);
            
            %Update metadata with image information
            [~, ~, fext] = fileparts(reader.filename);

            if ismember(fext, {'.tif', '.tiff'})

                %TIFF files do not have metadata
                linker = updateMetadata(linker, ...
                    'pxSize', NaN, ...
                    'pxSizeUnits', NaN, ...
                    'timestamps', NaN, ...
                    'timestampUnits', NaN, ...
                    'meanDeltaT', NaN);

            else

                [ts, tsunits] = getTimestamps(reader, 1, settings.ChannelToSegment);

                linker = updateMetadata(linker, ...
                    'pxSize', reader.pxSize, ...
                    'pxSizeUnits', reader.pxUnits, ...
                    'timestamps', ts, ...
                    'timestampUnits', tsunits, ...
                    'meanDeltaT', mean(diff(ts)));

            end

            tracks = linker.tracks;            
            save(fullfile(outputDir, [fn, '.mat']), 'tracks');
            
        end
        
        function [mask, bpInd, Isub] = segment(I, settings)
            %SEGMENT  Segment an image and return mask and branchpoints
            %
            %  [MASK, BP] = SEGMENT(I, SETTINGS) segments the image I with
            %  the SETTINGS in the struct provided. The outputs are a MASK
            %  and the branchpoint indices BP.
            
            Isub = imtophat(I, ones(10));
            
            FM = fibermetric(Isub, [2, 4, 6]);
            mask = FM > (max(FM(:)) / settings.ThresholdFactor);
            
            %Skeletonize the fibers to get the center line
            mask = bwskel(mask, 'MinBranchLength', settings.MinBranchLength);  %Remove small branches
            
            %Detect branching lines - objects with branches will be removed
            %in the final mask
            bpMask = bwmorph(mask, 'branchpoints');
            bpInd = find(bpMask);
   
            %--DEBUG
            %figure;
            %Iout = showoverlay(Isub, mask);
            %showoverlay(Iout, bpMask, 'Color', [1 0 0]);
            
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
        
    methods (Static)
        
        function varargout = showoverlay(img, mask, varargin)
            %SHOWOVERLAY  Overlays a mask on to a base image
            %
            %  SHOWOVERLAY(I, M) will overlay mask M over the image I, displaying it in
            %  a figure window.
            %
            %  C = SHOWOVERLAY(I, M) will return the composited image as a matrix
            %  C. This allows multiple masks to be composited over the same image. C
            %  should be of the same class as the input image I. However, if the input
            %  image I is a double, the output image C will be normalized to between 0
            %  and 1.
            %
            %  Optional parameters can be supplied to the function to modify both the
            %  color and the transparency of the masks:
            %
            %     'Color' - 1x3 vector specifying the color of the overlay in
            %               normalized RGB coordinates (e.g. [0 0 1] = blue)
            %
            %     'Transparency' - Value between 0 - 100 specifying the alpha level of
            %                      the overlay
            %
            %  Examples:
            %
            %    %Load a test image
            %    testImg = imread('cameraman.tif');
            %
            %    %Generate a masked region
            %    maskIn = false(size(testImg));
            %    maskIn(50:70,50:200) = true;
            %
            %    %Store the image to a new variable
            %    imgOut = SHOWOVERLAY(testImg, maskIn);
            %
            %    %Generate a second mask
            %    maskIn2 = false(size(testImg));
            %    maskIn2(100:180, 50:100) = true;
            %
            %    %Composite and display the second mask onto the same image as a
            %    %magenta layer with 50% transparency
            %    SHOWOVERLAY(imgOut, maskIn2, 'Color', [1 0 1], 'Transparency', 50);
            
            % Author: Jian Wei Tay (jian.tay@colorado.edu)
            % Version 2018-Feb-01
            
            ip = inputParser;
            ip.addParameter('Color',[0 1 0]);
            ip.addParameter('Opacity',100);
            ip.addParameter('Normalize', true);
            ip.parse(varargin{:});
            
            alpha = ip.Results.Opacity / 100;
            
            %Get the original image class
            imageClass = class(img);
            imageIsInteger = isinteger(img);
            
            % %Process the input image
            if ip.Results.Normalize
                img = double(img);
                img = img ./ max(img(:));
            end
            
            if size(img,3) == 1
                %Convert into an RGB image
                img = repmat(img, 1, 1, 3);
            elseif size(img,3) == 3
                %Do nothing
            else
                error('showoverlay:InvalidInputImage',...
                    'Expected input to be either a grayscale or RGB image.');
            end
            
            %Process the mask
            mask = double(mask);
            mask = mask ./ max(mask(:));
            
            if size(mask,3) == 1
                %Convert mask into an RGB image
                %mask = repmat(mask, 1, 1, 3);
                
                replacePx = mask ~= 0;
                
                for iC = 1:3
                    %mask(:,:,iC) = mask(:,:,iC) .* ip.Results.Color(iC);
                    
                    currC = img(:,:,iC);
                    currC(replacePx) = (currC(replacePx) .* (1 - alpha)) + (mask(replacePx) .* alpha .* ip.Results.Color(iC));
                    
                    img(:,:,iC) = currC;
                end
                
                
            elseif size(mask,3) == 3
                
                %Make the composite image
                replacePx = mask ~= 0;
                
                img(replacePx) = img(replacePx) .* (1 - alpha) + mask(replacePx) .* alpha;
                
            else
                error('showoverlay:InvalidMask',...
                    'Expected mask to be either a logical or RGB image.');
            end
            
            
            %Recast the image into the original image class
            if imageIsInteger
                multFactor = double(intmax(imageClass));
                img = img .* multFactor;
                img = cast(img, imageClass);
            else
                %multFactor = 1;
            end
            
            
            
            %Produce the desired outputs
            if nargout == 0
                imshow(img,[])
            else
                varargout = {img};
            end
            
        end
        
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
    
    
end














