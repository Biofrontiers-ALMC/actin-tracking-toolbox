classdef ActinData < TrackArray
    
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
        
        function obj = analyze(obj)
            %ANALYZE
            
            %Compute the:
            % 1. Instantaneous velocity
            % 2. Average Instantaneous velocity
            % 3. Filament length (avg?)
            % 3. Total number of frames
            % 4. If fiber is stuck
            
            for iT = 1:numel(obj.Tracks)
                
                %AD.Tracks(1).Data.Centroid{1}
                
                %Compute instantaneous speed
                pos = cell2mat(obj.Tracks(iT).Data.Centroid') .* obj.FileMetadata.PxSize;
                
                instantSpd = [0; sqrt(sum((diff(pos, 1)).^2, 2))];
                
                for jj = 1:numel(instantSpd)
                    obj.Tracks(iT).Data.InstantSpeed{jj} = instantSpd(jj) / obj.FileMetadata.DeltaT;
                end
                
                obj.Tracks(iT).MeanInstantSpeed = mean(instantSpd(2:end));
                
                len = cell2mat(obj.Tracks(iT).Data.MajorAxisLength') * obj.FileMetadata.PxSize;
                obj.Tracks(iT).MeanFilamentLength = mean(len);
                
                obj.Tracks(iT).NumFrames = numel(obj.Tracks(iT).Frames);
                
            end
            
            
        end
        
    end
    
    methods (Access = protected)
       
        function exportToCSV(obj, fn)
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
            
            fprintf(fid, 'Track ID, Num Frames, Mean Instantaneous Speed, Mean Filament Length, Frame');
            
            fprintf(fid, ', %s', datafields{:});
            fprintf(fid, '\n');
            
            for iTrack = 1:numel(obj.Tracks)
                
                %Print track metadata
                fprintf(fid, '%.0f, %.0f, %.3f, %.3f', ...
                    obj.Tracks(iTrack).ID, ...
                    obj.Tracks(iTrack).NumFrames, ...
                    obj.Tracks(iTrack).MeanInstantSpeed, ...
                    obj.Tracks(iTrack).MeanFilamentLength);
                    
                %Print track data
                for iF = 1:numel(obj.Tracks(iTrack).Frames)
                    
                    if iF > 1
                        fprintf(fid, ', , ,');
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
                fprintf(fid, '\n');
                
            end
            
            fclose(fid);
                        
        end
        
        
    end
    
    
end