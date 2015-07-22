% intensityAnalysis(data, varargin) outputs graphs looking at pit intensities

% Zach Weinberg 7/22/2015

function intensityAnalysis(data, varargin)

% load in lifetimeData.mat for data. Most of this is from runSlaveClassification

% is a list of movies using data as returned by loadConditionData()
parfor i = 1:length(data)
  %load([data(i).source 'Tracking' filesep 'ProcessedTracks.mat']);
  %load([data(i).source 'Analysis' filesep 'LifetimeData.mat']);
  LifeTimeData = load('LifeTimeData.mat')

  for t = 1:length(lifetime_s)
    %if tracks(i).catIdx = 1 | tracks(i).catIdx = 2 | tracks(i).catIdx = 5
      intensity = nanmean(tracks(t).A)
      category = track(t).catIdx
  end
end



