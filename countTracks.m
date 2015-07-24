%runSlaveChannelClassification(data, varargin) identifies trajectories with significant slave channel fluorescence
%
% Input:
%          data : structure returned by loadConditionData()
%
% Options ('specifier', value):
%          'np' : number of points to use for randomized detections
%    'Cutoff_f' : minimum track length to consider for classification (in frames)
%
% Notes: This function modifies the output of runTrackProcessing(),
%        by default saved in Tracking/ProcessedTracks.mat

% Francois Aguet, October 2010 (last modified: 10/09/2012)

function countTracks(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addParamValue('ChannelNames', []);
ip.parse(varargin{:});

filename = 'ProcessedTracks.mat';
nd = numel(data);

for i=1:nd
  % load tracks (all)
  ts = load([data(i).source 'Tracking' filesep filename]);
  tracks = [ts.tracks];
  catFilter = find([ts.tracks.catIdx] == 1 | [ts.tracks.catIdx] == 5);

  total = numel(tracks(catFilter));
  significantMasters = [ts.tracks.significantMaster];
  significantSlaves = [ts.tracks.significantSlave];

  significantMasters = significantMasters(:,catFilter);
  significantSlaves = significantSlaves(:,catFilter);

  outputStr = '\nTotal tracks detected in cell %d: %d\n';

  %(TP)load cellmask
  maskPath = [data(i).source 'Detection' filesep 'cellmask.tif'];
  mask = double(imread(maskPath));
  mSize = sqrt(numel(mask(find(mask))));
  mSizeStr = 'Cell size is %d pixels^2\n';

  %(TP)normalize total pitcount to cellmask size
  normTracks = total/mSize;
  normStr = 'Normalized total track count: %d tracks/pixels^2\n';

  for c=2:numel(data(i).channels)
      counts(c) = numel(find(significantMasters(c,:,:) == 1));
      normCounts(c) = counts(c)/mSize;
      if isempty(ip.Results.ChannelNames)
         outputStr = [outputStr 'Significant tracks in channel ' c ': %d\n'];
         normStr = [normStr 'Normalized track count in channel' c ': %d tracks/pixels^2\n'];
      else
         outputStr = [outputStr ip.Results.ChannelNames{c} '-positive tracks: %d\n'];
         normStr = [normStr ip.Results.ChannelNames{c} '-normalized positive track count: %d tracks/pixels^2\n'];
      end
  end

  fprintf([outputStr mSizeStr normStr], [i total counts(2:end)], mSize, [normTracks normCounts(2:end)])
end
