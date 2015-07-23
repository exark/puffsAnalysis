% intensityAnalysis(data, varargin) outputs graphs looking at pit intensities

% Zach Weinberg 7/22/2015

function intensityAnalysis(data, varargin)

% load in lifetimeData.mat for data. Most of this is from runSlaveClassification
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
%ip.addOptional('xl', []);
%ip.addOptional('xa', []);
ip.addParamValue('ExcludeVisitors', false, @islogical);
ip.addParamValue('Cutoff_f', 1, @isscalar);
ip.addParamValue('FirstNFrames', [], @isvector);
ip.addParamValue('DisplayFunction', @sqrt);
ip.addParamValue('Channel', 1, @isposint);
ip.addParamValue('Legend', []);
%ip.addParamValue('Parent', []);
ip.addParamValue('LifetimeData', 'LifetimeData.mat');
%ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat');
ip.addParamValue('PlotIndividual', false, @islogical);
ip.addParamValue('NormX', true, @islogical);
ip.addParamValue('FontSize', 10);
ip.addParamValue('Width', 4, @isposint);
ip.addParamValue('AmplitudeCorrection', []);
ip.addParamValue('allInt',[]); 
ip.addParamValue('allCat',[]);
ip.addParamValue('allLT',[]);
ip.parse(data, varargin{:});

% (TP)  taken from plotMaxIntVsLifetime
if ip.Results.PlotIndividual
    lftData = getLifetimeData(data, 'Overwrite', false, 'Mask', true,...
        'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
        'ReturnValidOnly', false, 'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
    data = arrayfun(@(i) i, data, 'unif', 0);
    lftData = arrayfun(@(i) i, lftData, 'unif', 0);
    nd = numel(data); %(TP) with data as an array, nd = number of cells
    
else
    if ~iscell(data)
        data = {data};
    end
    nd = numel(data); %(TP) with data as a cell, all the cells' info is put into 1 cell, -> nd always = 1
    lftData = cell(1,nd);
    for i = 1:nd
        lftData{i} = getLifetimeData(data{i}, 'Overwrite', false, 'Mask', true,...
            'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
            'ReturnValidOnly', false, 'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
    end
end
% taken from plotMaxVsLifetime

parfor i = 1:length(data) %numcells
  %load([data(i).source 'Tracking' filesep 'ProcessedTracks.mat']);
  %load([data(i).source 'Analysis' filesep 'LifetimeData.mat']);
  LifetimeData = load('LifetimeData.mat')

  for t = 1:numel(data(i)) %numtracks
    int = LifetimeData(t).A_all
    cat = LifetimeData(t).catIdx
    lt = LifetimeData(t).lifetime_s
    meanIntensity = nanmean(int)
    allInt(i,t) = meanIntensity 
    allCat(i,t) = cat 
    allLT(i,t) = lt
  end
end
scatter(



