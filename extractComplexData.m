function [rawResults, normalizedMeanResults, normalizedMedianResults, nPuffs] = extractComplexData(dataStruct) 
    rawResults = [];
    normalizedMeanResults = [];
    normalizedMedianResults = [];
    nPuffs = [];
    
    for i=1:numel(dataStruct)
        load([dataStruct(i).source filesep 'Tracking/ProcessedTracks.mat'], 'tracks');
        load([dataStruct(i).source filesep 'Classification/RFResults.mat'], 'puffs');
        rR = [];
        nmeanR = [];
        nmedianR = [];
        for j=puffs
            
            [peak, time_to_peak] = max(tracks(j).Ac);
            
            if ~isempty(tracks(j).endBuffer)
                Ac_with_buffer = [tracks(j).Ac (tracks(j).endBuffer.A+tracks(j).endBuffer.c)];
            else
                Ac_with_buffer = tracks(j).Ac;
            end
            
            tau_one_half = NaN;
            for x=time_to_peak:numel(Ac_with_buffer)
                if Ac_with_buffer(x) <= 0.5*peak
                    tau_one_half = x - time_to_peak;
                end
            end
            
            plateau = 1;
            for y=time_to_peak:numel(Ac_with_buffer)
                if Ac_with_buffer(x) >= 0.9*peak
                    plateau = plateau +1;
                end
            end
            
            % cell_num - lifetime - int_density - start frame -
            %   time_to_peak tau1/2 plateau
            datum = [i tracks(j).lifetime_s trapz(tracks(j).Ac) tracks(j).start ((time_to_peak*0.1) - 0.1) tau_one_half*0.1 plateau*0.1 mean(tracks(j).x)/dataStruct(i).imagesize(2) mean(tracks(j).y)/dataStruct(i).imagesize(1)];
            % - 1, cell number
            % - 2, lifetime
            % - 3, int density
            % - 4, track start frame
            % - 5, time to peak
            % - 6, t1/2
            % - 7, plateau
            rR = [rR; datum];
            nmeanR = [nmeanR; datum];
            nmedianR = [nmedianR; datum];
        end
        for k = 1:size(nmedianR,2)
            nmedianR(:,k) = nmedianR(:,k)./median(nmedianR(:,k));
            nmeanR(:,k) = nmeanR(:,k)./mean(nmeanR(:,k));
        end
        rawResults = [rawResults; rR];
        normalizedMeanResults = [normalizedMeanResults; nmeanR];
        normalizedMedianResults = [normalizedMedianResults; nmedianR];
        nPuffs = [nPuffs; numel(puffs)];
    end

end