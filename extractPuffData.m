function [rawResults, normalizedMeanResults, normalizedMedianResults] = extractPuffData(dataStruct, datumString) 
    rawResults = [];
    normalizedMeanResults = [];
    normalizedMedianResults = [];
    
    for i=1:numel(dataStruct)
        load([dataStruct(i).source filesep 'Tracking/ProcessedTracks.mat'], 'tracks');
        load([dataStruct(i).source filesep 'Classification/RFResults.mat'], 'puffs');
        rR = [];
        nmeanR = [];
        nmedianR = [];
        for j=puffs
            eval(['datum = [' datumString '];']);
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
    end

end