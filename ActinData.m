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
            
            if numel(data) == 1 && isstruct(tmp.(data{1}))
                
                obj.LastID = tmp.(data{1}).LastID;
                obj.Tracks = tmp.(data{1}).Tracks;
                obj.FileMetadata = tmp.(data{1}).FileMetadata;
                obj.CreatedOn = tmp.(data{1}).CreatedOn;
                
            else
                error('Expected data to be a struct. Other formats not currently supported.')
                
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
end