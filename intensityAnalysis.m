% intensityAnalysis(data, varargin) outputs graphs looking at pit intensities

% Zach Weinberg, Tiffany Phan 7/22/2015

function intensityAnalysis(data, varargin)

% load in lifetimeData.mat for data. Most of this is from runSlaveClassification
ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data');
ip.addParamValue('Cutoff_f', 1, @isscalar);
ip.addParamValue('FirstNFrames', [], @isvector);
ip.addParamValue('LifetimeData', 'LifetimeData.mat');
ip.addParamValue('ProcessedTracks', 'ProcessedTracks.mat');
ip.addParamValue('AmplitudeCorrection', []);
ip.addParamValue('PlotAllCat',false, @islogical); %(TP) new addition
ip.parse(data, varargin{:});

% (TP) for a single cell at a time
load([data.source 'Analysis' filesep 'LifetimeData.mat']);%opts.LifetimeData]);
nt = numel(A_all); %(TP) number of tracks

for t = 1:nt %for every track
    if ~isempty(ip.Results.FirstNFrames)
        tmp = nanmean(A_all{t}(ip.Results.Cutoff_f:ip.Results.FirstNFrames));
    else
        tmp = nanmean(A_all{t}(ip.Results.Cutoff_f:end)); 
    meanA_all(t) = tmp;
    end
end

if ip.Results.PlotAllCat
    categories = []
    for c = 1:8
        if numel(find(catIdx == c)) ~= 0
            categories = [categories,c]
        end 
    end 
else 
   categories = [1,2,5];
end 

legendStr = {}
for c = 1:numel(categories)
% if numel(find(catIdx == categories(c))) ~= 0
 hold on;
 scatter(lifetime_s(catIdx == categories(c)),meanA_all(catIdx == categories(c)));
 legendStr = [legendStr, {num2str(categories(c))}]
end 
legend(legendStr);

% (TP) Everything from HERE DOWN is for numerous cells at a time. In
% progress. Need to figure out why all catIdx is 1 in lftData.

% nd = numel(data); %(TP) nd = number of cells in data
% %lftData = cell(1,nd);
% if ~iscell(data)
%     data = {data};
% end 
% 
% %(TP) each cell in lftData contains all the lifetime data for one cell
% for i = 1:nd
%     [lftData{i},~] = getLifetimeData(data{i}, 'Overwrite', false, 'Mask', true,...
%         'ProcessedTracks', ip.Results.ProcessedTracks, 'LifetimeData', ip.Results.LifetimeData,...
%         'ReturnValidOnly', false, 'AmplitudeCorrectionFactor', ip.Results.AmplitudeCorrection);
% end
% 
% meanA_all = cell(nd,1);
% 
% for i = 1:nd %(TP)for every cell
%     nt = numel(lftData{i}.A_all); %(TP) number of track
%     for t = 1:nt %for every track
%         if ~isempty(ip.Results.FirstNFrames)
%             tmp = nanmean([lftData{i}.A_all{t,ip.Results.Cutoff_f:ip.Results.FirstNFrames,1}]);
%         else
%             tmp = nanmean([lftData{i}.A_all{t,ip.Results.Cutoff_f:end,1}]); 
%         meanA_all{i}{t} = tmp;
%         end
%     end
% end
% 
% for i = 1:nd 
%     if ip.Results.PlotAllCat
%         catFilter = [find(lftData{i}.catIdx_all == 1)' find(lftData{i}.catIdx_all == 2)'...
%             find(lftData{i}.catIdx_all == 3)' find(lftData{i}.catIdx_all == 4)'...
%             find(lftData{i}.catIdx_all == 5)' find(lftData{i}.catIdx_all == 6)'...
%             find(lftData{i}.catIdx_all == 7)' find(lftData{i}.catIdx_all == 8)']';
%         categories = [1:8];
%     else 
%         catFilter = [find(lftData{i}.catIdx_all == 1)' find(lftData{i}.catIdx_all == 2)'...
%             find(lftData{i}.catIdx_all == 5)']';
%         categories = [1,2,5];
%     legendStr = (num2str(categories','Category %-d'));
%     scatter(lftData{i}.lifetime_s(catFilter),(meanA_all(1))(catFilter));
%     legend(legendStr);
%     end 
% end 


   
