%NOTE: This script has additional code to filter filament data by length
%and number of frames tracked


classdef ActinData_LAL2_JWT < TrackArray
    %ACTINDATA  Analyze actin data
    %
    %  OBJ = ACTINDATA creates a new ActinData object to analyze data from
    %  the tracked movies. Data should be imported into the project using
    %  the method importdata.
    %
    
    properties
        
        maxStuckSpeed = 0.06;
        
        timestamps = [];
        
    end
    
    methods
        
        function obj = importdata(obj, filename)
            
            %IMPORTDATA  Import data into the analyzer object
            %
            %  OBJ = IMPORTDATA(OBJ, FILENAME) imports data from the
            %  MAT-file and carries out additional data cleaning and
            %  processing steps (e.g. assigning colony and generation
            %  numbers, computing growth rate)
            
            if ~exist('filename', 'var')
                
                [filename, pathname] = uigetfile({'*.mat', 'MAT-file'}, ...
                    'Select file(s) to process');
                
                if isequal(filename, 0) || isequal(pathname, 0)
                    %Cancel
                    return;
                end
                
                tmp = load(fullfile(pathname, filename));
                
            else
                
                tmp = load(filename);
            end
            
            data = fieldnames(tmp);
            
            if numel(data) == 1 && isa(tmp.(data{1}), 'TrackArray')
                
                obj = importobj(obj, tmp.(data{1}));
                
            else
                error('Expected data to be a TrackArray object. Other formats not currently supported.')
                
            end
            
            
            
        end
        
        function obj = analyze(obj,minFilamentLength,numFrames)
            %ANALYZE  Analyze the data
            %
            %  OBJ = ANALYZE(OBJ) updates the track array with additional
            %  data fields:
            %    * Instantaneous speed per frame (um/s)
            %    * Number of frames
            %    * Average speed (um/s)
            %    * Standard deviation of speed (ums)
            %    * Average filament length (um)
            %    * Standard deviation of length (um)
            %    * Total distance traveled (um)
            %    * Total displacement (distance between first and last
            %      frames) (um)
            
            for iT = 1:numel(obj.Tracks)
                
                %-- Compute instantaneous speed --
                %Get position of cells and concatenate into a matrix
                pos = cell2mat(obj.Tracks(iT).Data.Centroid') .* obj.FileMetadata.pxSize(1);
                
                %Compute the distance traveled between frames
                distBetweenFrames = [0; sqrt(sum((diff(pos, 1)).^2, 2))];
                
                %Compute the instaneous speed
                if size(pos, 1) > 1
                    for jj = 1:numel(distBetweenFrames)
                        if jj == 1
                            obj.Tracks(iT).Data.InstantaneousSpeed{jj} = NaN;
                        else
                            
                            if ~isempty(obj.timestamps)
                                %Compute dTime
                                try
                                    
                                    dTime = obj.timestamps(obj.Tracks(iT).Frames(jj)) - obj.timestamps(obj.Tracks(iT).Frames(jj - 1));
                                    
                                catch
                                    keyboard
                                end
                                
                                obj.Tracks(iT).Data.InstantaneousSpeed{jj} = distBetweenFrames(jj) / dTime ;
                            else
                                obj.Tracks(iT).Data.InstantaneousSpeed{jj} = distBetweenFrames(jj) / obj.FileMetadata.meanDeltaT;
                            end
                        end
                    end
                else
                    obj.Tracks(iT).Data.InstantaneousSpeed = {NaN};
                end
                
                %-- Compute Statistics --
                %Number of frames
                obj.Tracks(iT).NumFrames = numel(obj.Tracks(iT).Frames);
                
                %Average speed
                obj.Tracks(iT).AverageSpeed = mean([obj.Tracks(iT).Data.InstantaneousSpeed{:}], 'all', 'omitnan');
                
                %Standard deviation of speed
                obj.Tracks(iT).StdDevSpeed = std([obj.Tracks(iT).Data.InstantaneousSpeed{:}], 0, 'all', 'omitnan');
                
                %Average filament length
                len = cell2mat(obj.Tracks(iT).Data.Length') * obj.FileMetadata.pxSize(1);
                obj.Tracks(iT).AverageFilamentLength = mean(len);
                
                %Standard deviation of filament length
                obj.Tracks(iT).StdDevFilamentLength = std(len);
                
                %Total distance traveled
                obj.Tracks(iT).DistanceTraveled = sum(distBetweenFrames);
                
                %Total displacement
                obj.Tracks(iT).Displacement = sqrt(sum((obj.Tracks(iT).Data.Centroid{1} - obj.Tracks(iT).Data.Centroid{end}).^2, 2));
                
                if obj.Tracks(iT).AverageSpeed < obj.maxStuckSpeed
                    obj.Tracks(iT).isStuck = true;
                else
                    obj.Tracks(iT).isStuck = false;
                end
                
                if ~obj.Tracks(iT).isStuck && length(obj.Tracks(iT).Frames) >= numFrames && obj.Tracks(iT).AverageFilamentLength > minFilamentLength
                    obj.Tracks(iT).filtered = true;
                else
                    obj.Tracks(iT).filtered = false;
                end
            end
        end
        
        function Iout = showlabels(obj, frame ,trackId)
            singleTrack = trackId > 0;
            bfr = BioformatsImage(obj.FileMetadata.filename);
            
            I = getPlane(bfr, 1, 1, frame);
            
            maskStuck = false(size(I));
            maskNotStuck = false(size(I));
            
            for iT = 1:numel(obj)
                
                idx = find(obj.Tracks(iT).Frames == frame);
                
                if ~isempty(idx)
                    
                    if obj.Tracks(iT).isStuck
                        color = 'y';
                        if singleTrack
                            if iT == trackId
                                color = 'g';
                            end
                        end
                        maskStuck(obj.Tracks(iT).Data.PixelIdxList{idx}) = true;
                        I = insertText(I, obj.Tracks(iT).Data.Centroid{idx}, iT, ...
                            'BoxOpacity', 0, 'TextColor', color);
                        
                    else
                        color = 'w';
                        if singleTrack
                            if iT == trackId
                                color = 'g';
                            end
                        end
                        maskNotStuck(obj.Tracks(iT).Data.PixelIdxList{idx}) = true;
                        I = insertText(I, obj.Tracks(iT).Data.Centroid{idx}, iT, ...
                            'BoxOpacity', 0, 'TextColor', color);
                        
                    end
                end
            end
            
            if any(maskStuck, 'all')
                I = ActinTracker.showoverlay(I, bwperim(maskStuck), 'color', [1 0 0]);
            end
            
            if any(maskNotStuck, 'all')
                I = ActinTracker.showoverlay(I, bwperim(maskNotStuck), 'color', [0 0 1]);
            end
            
            if nargout == 1
                Iout = I;
            else
                imshow(I)
            end
            
        end
        
    end
        methods (Access = protected)
            
            function exportToCSV(obj, fn, new)
                %EXPORTTOCSV  Export track and filemetadata as CSV files
                %
                %  EXPORTTOCSV(OBJ, FN) exports the data as a CSV file.
                
                fid = fopen(fn, 'w');
                
                if fid < 0
                    error('Error opening file %s for writing.', fn);
                end
                
                %Print filemetadata if not empty
                if ~isempty(obj.FileMetadata)
                    metadataFields = fieldnames(obj.FileMetadata);
                    
                    for iMD = 1:numel(metadataFields)
                        
                        switch class(obj.FileMetadata.(metadataFields{iMD}))
                            
                            case 'char'
                                fprintf(fid, '%s, %s\n', metadataFields{iMD}, obj.FileMetadata.(metadataFields{iMD}));
                                
                            case 'double'
                                fprintf(fid, '%s, %s\n', metadataFields{iMD}, mat2str(obj.FileMetadata.(metadataFields{iMD})));
                                
                        end
                    end
                    fprintf(fid, '\n');
                end
                
                %Print column headers
                datafields = fieldnames(obj.Tracks(1).Data);
                
                fprintf(fid, ['Track ID, Num Frames, Average Speed, StdDev Speed,'...
                    'Average Filament Length, StdDev Filament Length, Distance Traveled, Displacement, Stuck, Filtered, Frame']);
                
                fprintf(fid, ', %s', datafields{:});
                fprintf(fid, '\n');
                
                for iTrack = 1:numel(obj.Tracks)
                    
                    %Print track metadata
                    fprintf(fid, '%.0f, %.0f, %.3f, %.3f, %.3f, %.3f, %.3f, %.3f,%.3f,%.3f', ...
                        obj.Tracks(iTrack).ID, ...
                        obj.Tracks(iTrack).NumFrames, ...
                        obj.Tracks(iTrack).AverageSpeed, ...
                        obj.Tracks(iTrack).StdDevSpeed, ...
                        obj.Tracks(iTrack).AverageFilamentLength, ...
                        obj.Tracks(iTrack).StdDevFilamentLength, ...
                        obj.Tracks(iTrack).DistanceTraveled, ...
                        obj.Tracks(iTrack).Displacement,...
                        obj.Tracks(iTrack).isStuck,...
                        obj.Tracks(iTrack).filtered);
                    if ~new
                        %Print track data
                        for iF = 1:numel(obj.Tracks(iTrack).Frames)
                            
                            if iF > 1
                                fprintf(fid, ', , , , , , , ,');
                            end
                            
                            %Print frame index
                            fprintf(fid, ', %.0f', obj.Tracks(iTrack).Frames(iF));
                            
                            %Print data fields
                            for iP = 1:numel(datafields)
                                
                                switch class(obj.Tracks(iTrack).Data.(datafields{iP}){iF})
                                    
                                    case 'char'
                                        fprintf(fid, ', %s', obj.Tracks(iTrack).Data.(datafields{iP}){iF});
                                        
                                    case 'double'
                                        fprintf(fid, ', %s', mat2str(obj.Tracks(iTrack).Data.(datafields{iP}){iF}));
                                        
                                end
                                
                            end
                            fprintf(fid, '\n');
                        end
                    end
                    fprintf(fid, '\n');
                    
                end
                
                fclose(fid);
                
            end
            
            
        end
        
        
    end