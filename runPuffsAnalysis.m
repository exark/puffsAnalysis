%[puffsRes] = runPuffsAnalysis(data, varargin) computes stats for puffs
%
% Inputs:
%               data : structure returned by loadConditionData()
%
%
% Outputs:
%             puffsRes : structure containing number of puffs, nonpuffs and cell areas.

% Zara Weinberg

function [puffsRes] = runPuffsAnalysis(data, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @(x) isstruct(x) && numel(unique([data.framerate]))==1);
ip.addParamValue('RFResults', 'RFResults.mat', @ischar);

ip.parse(data, varargin{:});

mCh = find(strcmp(data(1).source, data(1).channels));
res = struct([]);

nd = numel(data);

puffsRes.cellArea = zeros(nd,1);
puffsRes.numPuffs = zeros(nd,1);
puffsRes.normalizedPuffs = zeros(nd,1);
for i = 1:nd

    px = data(i).pixelSize / data(i).M; % pixels size in object space
    mask = logical(getCellMask(data(i)));
    puffsRes.cellArea(i) = sum(mask(:)) * px^2 * 1e12; % in µm^2

    % Load and count number of puffs
    rpath = [data(i).source 'Classification' filesep ip.Results.RFResults];
    load(rpath);

    puffsRes.numPuffs(i) = numel(puffs);
    puffsRes.normalizedPuffs(i) = puffsRes.numPuffs(i)/puffsRes.cellArea(i); % in puffs/µm^2

end

figure();
bar(puffsRes.normalizedPuffs);
