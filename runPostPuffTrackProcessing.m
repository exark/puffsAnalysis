function runPostPuffTrackProcessing (data, varargin)

ip.addParamValue('Overwrite', false, @islogical);
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat', @ischar);

tPath = [data.source 'Tracking' filesep opts.ProcessedTracks];
if exist(tPath, 'file')==2
    tracks = load([data.source 'Tracking' filesep 'ProcessedTracks.mat']);
else
    fprintf('runPostPuffTrackProcessing: no processed tracks data found for %s\n', getShortPath(data));
    return;
end



