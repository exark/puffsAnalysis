% puffsModelChecker(data, puffs, nonpuffs, varargin) displays movies with associated detection and tracking results.
%
% Inputs:
%             data : single movie structure returned by loadConditionData.m
%
% Options (specifier/value format):
%     'NumSamples' : Number of tracks to show from each condition. Default: 8
%     'LoadTracks' : {true}|false specifies whether tracking data is loaded
%     'LoadFrames' : {true}|false specifies whether the raw movie data is loaded
%                    This is required for visualization of tracks overlaid on data
%       'LoadMask' : {true}| false specifies whether the cell outline mask is loaded
%       'Cutoff_f' : Minimum length of tracks to load, in frames. Default: 3
%       'Subset'   : matrix of subset track indices to load
%

% Zara Weinberg, 2016

function puffsModelChecker(data, puffs, nonpuffs, varargin)

ip = inputParser;
ip.CaseSensitive = false;
ip.addRequired('data', @isstruct);
% ip.addRequired('puffs', @ismatrix);
% ip.addRequired('nonpuffs', @ismatrix);
ip.addParamValue('NumSamples', 12, @isscalar);
ip.addParamValue('LoadTracks', true, @islogical);
ip.addParamValue('LoadFrames', true, @islogical);
ip.addParamValue('LoadMask', true, @islogical);
ip.addParamValue('RelativePath', 'Tracking', @ischar);
ip.addParamValue('Cutoff_f', 3);
ip.addParamValue('Subset', []);
ip.parse(data, varargin{:});

% Handles/settings are stored in 'appdata' of the figure handle
handles.data = data;

% detect number of channels (up to 4)
nCh = length(data.channels);
if nCh>4
    error('Max. 4 channels supported.');
end

handles.nCh = nCh;
% master channel index
handles.mCh = find(strcmp(data.source, data.channels));


fidx = 1;
tcur = 1;
xs = [];
ys = [];

nx = data.imagesize(2);
ny = data.imagesize(1);
nf = data.movieLength;

lcolor = hsv2rgb([0.55 0.5 0.8]);

%===============================================================================
% Load movie and associated analysis results
%===============================================================================
% readfct = @(path, i) imread(path, i);

stack = cell(1,nCh);
if ip.Results.LoadFrames
    fprintf('Loading frames ... ');
    if ~iscell(data.framePaths{1})
        for c = 1:nCh
            %stack{c} = readtiff(data.framePaths{c});
            stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
            for i = 1:data.movieLength
                stack{c}(:,:,i) = imread(data.framePaths{c}, i);
            end
        end
    else
        for c = 1:nCh
            stack{c} = zeros([data.imagesize data.movieLength], 'uint16');
            for i = 1:data.movieLength
                stack{c}(:,:,i) = imread(data.framePaths{c}{i});
            end
        end
    end
    fprintf('done.\n');
end

%-------------------------------------------------------------------------------
% Load detection masks
%-------------------------------------------------------------------------------
dpath = [data.source 'Detection' filesep 'Masks'];
if exist(dpath, 'dir')==7 && ip.Results.LoadFrames
    fprintf('Loading detection masks ... ');
    dmask = zeros(ny,nx,nf, 'uint8');
    if ~iscell(data.framePaths{1})
        for i = 1:nf
            dmask(:,:,i) = imread(data.maskPaths, i);
        end
    else
        for i = 1:data.movieLength
            dmask(:,:,i) = imread(data.maskPaths{i});
        end
    end
    fprintf('done.\n');
else
    dmask = [];
end

if exist([data.source 'Detection' filesep 'cellmask.tif'], 'file')==2 && ip.Results.LoadMask
    cellMask = imread([data.source 'Detection' filesep 'cellmask.tif']);
else
    cellMask = [];
end

%-------------------------------------------------------------------------------
% Load detection files
%-------------------------------------------------------------------------------
if ip.Results.LoadFrames
    % for c = 1:nCh
    detectionFile = [data.channels{1} 'Detection' filesep 'detection_v2.mat'];
    if (exist(detectionFile, 'file')==2)
        frameInfo = load(detectionFile);
        frameInfo = frameInfo.frameInfo;
    else
        frameInfo = [];
    end
    % end
end

%-------------------------------------------------------------------------------
% Load tracks
%-------------------------------------------------------------------------------
if ip.Results.LoadTracks
    fprintf('Loading tracks ... ');
    tracks = [];
    bgA = [];
    maxA = [];
    % identify track file
    fileList = dir([data.source ip.Results.RelativePath filesep 'ProcessedTracks*.mat']);
    fileList = {fileList.name};
    if numel(fileList)>1
        idx = 0;
        while ~(idx>=1 && idx<=numel(fileList) && round(idx)==idx)
            fprintf('Tracking results found for this data set:\n');
            for i = 1:numel(fileList)
                fprintf('[%d] %s\n', i, fileList{i});
            end
            idx = str2double(input('Please enter the number of the set to load: ', 's'));
        end
        fileName = fileList{idx};
    elseif numel(fileList)==1
        fileName = fileList{1};
    else
        fileName = [];
    end
end

if exist([data.source ip.Results.RelativePath filesep fileName], 'file')==2 && ip.Results.LoadTracks
    tmp = load([data.source ip.Results.RelativePath filesep fileName]);
    if isempty(ip.Results.Subset)
        tracks = tmp.tracks;
    else
        tracks=tmp.tracks(ip.Results.Subset);
    end
    if isfield(tmp, 'bgA')
        bgA = cellfun(@(i) prctile(i, 95, 2), tmp.bgA, 'unif', 0);
        bgA = [bgA{:}];
    end
    clear tmp;
   % (TP) tracks = tracks([tracks.lifetime_s] >= data.framerate*ip.Results.Cutoff_f);
    %[~, sortIdx] = sort([tracks.lifetime_s], 'descend');
    %tracks = tracks(sortIdx);

    % apply cell mask
    if ~isempty(cellMask)
        nt = numel(tracks);
        x = NaN(1,nt);
        y = NaN(1,nt);
        for t = 1:nt
            x(t) = round(nanmean(tracks(t).x(1,:)));
            y(t) = round(nanmean(tracks(t).y(1,:)));
        end
        idx = sub2ind([ny nx], y, x);
        tracks = tracks(cellMask(idx)==1);
    end
    nt = numel(tracks);
    selIndex = true(1,nt);

    nseg = [tracks.nSeg];

    np = sum(nseg);
    X = NaN(nf, np);
    Y = NaN(nf, np);
    G = false(nf, np);
    % for significance values, store vectors
    mvec = [tracks.hval_Ar];
    if isfield(tracks, 'significantVsBackground')
        svec = [tracks.significantVsBackground];
    else
        svec = [];
    end
    fvec = [tracks.f];
    xvec = [tracks.x];
    yvec = [tracks.y];

    % vector of start indexes since multiple segments/track
    tidx = cumsum([1 nseg(1:end-1)]);

    trackStarts = [tracks.start];
    trackEnds = [tracks.end];
    mu_x = NaN(1,nt);
    mu_y = NaN(1,nt);

    for t = 1:nt
        if nseg(t)==1
            X(tracks(t).f, tidx(t)) = tracks(t).x(1,:);
            Y(tracks(t).f, tidx(t)) = tracks(t).y(1,:);
            G(tracks(t).f, tidx(t)) = tracks(t).gapVect;
            mu_x(t) = nanmean(X(:,tidx(t)));
            mu_y(t) = nanmean(Y(:,tidx(t)));
        else
            sep = find(isnan(tracks(t).t));
            sep = [0 sep numel(tracks(t).f)+1]; %#ok<AGROW>
            for s = 1:tracks(t).nSeg
                sidx = sep(s)+1:sep(s+1)-1;
                X(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).x(1,sidx);
                Y(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).y(1,sidx);
                G(tracks(t).f(sidx), tidx(t)+s-1) = tracks(t).gapVect(sidx);
            end
            mu_x(t) = nanmean(nanmean(X(:,tidx(t):tidx(t)+s-1)));
            mu_y(t) = nanmean(nanmean(Y(:,tidx(t):tidx(t)+s-1)));
        end
    end

    % segment index
    % Example: [1 1 2 3 4 4 ... ] first two cols are from same track
    idx = diff([tidx size(X,2)+1]);
    idx = arrayfun(@(i) i+zeros(1, idx(i)), 1:numel(idx), 'unif', 0);
    tstruct.idx = [idx{:}];
    tstruct.n = numel(tracks);

    % min/max track intensities
    maxA = arrayfun(@(t) max(t.A, [], 2), tracks, 'unif', 0);
    maxA = [maxA{:}];
    maxInt = prctile(maxA, 99, 2);
    da = floor(log10(maxInt));
    % y-axis unit
    yunit = round(maxInt ./ 10.^da) .* 10.^(da-1);
    maxInt = ceil(maxInt./yunit) .* yunit;
end
fprintf('done.\n');

%===============================================================================
% Setup main GUI window/figure
%===============================================================================
hfig = figure('Units', 'normalized', 'Position', [0 0 1 1],...
    'PaperPositionMode', 'auto', 'Toolbar', 'none',...
    'Color', [1 1 1], 'Resize', 'on',...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', 'ModelChecker'); % switch back to [getDirFromPath(getExpDir(data)) filesep getCellDir(data)]

pos = get(hfig, 'Position'); % [pixels]

%===============================================================================
% Setup panels and axes for displaying tracks.
%===============================================================================

% Setup tracks in 4 columns, two columns of each condition.
cols = 6;
rows = ip.Results.NumSamples/(cols/2);

% increments for x and y positioning
xinc = pos(3)/cols;
yinc = (pos(4)-75)/rows;

% pick random tracks from each of puffs and nonpuffs
puffsToDisplay = reshape(datasample(puffs, 12, 'Replace', false), rows, (cols/2));
nonpuffsToDisplay = reshape(datasample(nonpuffs, 12, 'Replace', false), rows, (cols/2));

for i=1:cols
  for j=1:rows
    if i <= (cols/2)
      tcur = puffsToDisplay(j,i);
      bgcolor = [0 1 0];
    else
      tcur = nonpuffsToDisplay(j,i-(cols/2));
      bgcolor = [1 0 0];
    end

    ph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Track: ' num2str(tcur)], 'BackgroundColor', bgcolor,...
                 'Position', [(pos(1)+((i-1)*xinc)+5) (pos(1)+((j-1)*yinc)+5) (xinc-10) (yinc-1)]);
    ax = axes('Parent', ph, 'Box', 'on', 'Units', 'pixels');

    plotTrack(handles.data, tracks(tcur), 'Handle', ax);
    tPanels(i,j) = ph;
    tAxes(i,j) = ax;
  end
end

handles.tPanels = tPanels;
handles.tAxes = tAxes;

setappdata(hfig, 'handles', handles); % write 'handles' to hfig

end

%======================================
% Helper functions:
%======================================

function plotTrackMontageLocal(track, trackStack, xa, ya, ax, width)

[nc, nf] = size(trackStack);

fontName = 'Helvetica'
fontSize = 14;
ip.addParamValue('FramesPerRow', 20, @isscalar);
ip.addParamValue('ShowDetection', false, @islogical);
ip.addParamValue('ShowMarkers', false, @islogical);
ip.addParamValue('DynamicRange', []);
ip.parse(track, trackStack, varargin{:});
width = ip.Results.Width;
labels = ip.Results.Labels;

xa = ip.Results.xa;
ya = ip.Results.ya;
if isempty(xa)
    w = (size(trackStack{1,1},2)-1)/2;
    xa = -w:w;
end
if isempty(ya)
    w = (size(trackStack{1,1},1)-1)/2;
    ya = -w:w;
end

if ~isempty(labels)
    rgbColors = arrayfun(@(x) hsv2rgb([x 1 0.9]), getFluorophoreHues(labels), 'UniformOutput', false);
else
    rgbColors = [];
end

% dynamic range for each channel
if ~isempty(track.startBuffer)
    sb = numel(track.startBuffer.t);
else
    sb = 0;
end

% dynamic range of the entire stack
dRange = ip.Results.DynamicRange;
if isempty(dRange)
    dRange = zeros(nc,2);
    for c = 1:nc
        cstack = cat(3, trackStack{c,:});
        dRange(c,:) = [min(cstack(:)) max(cstack(:))];
    end
end

% nx is now the number of frames displayed per row
nx = ip.Results.FramesPerRow;


% display with fixed #frames/width
[wxi, dxi, dci, nr, height] = getProportions(width, nx, nf, nc);


hf = figure('Visible', 'off', 'PaperPositionMode', 'auto', 'Position', [50, 200, width, height], 'Color', 'w');
ha = axes('Position', [0 0 1 1], 'XLim', [0 width], 'YLim', [0 height]);
set(hf, 'Units', 'pixels');

% compute max. text length
textWidth = zeros(1,nc);
if ~isempty(labels)
    for c = 1:nc
        ht = text(0, 0, labels{c}, 'FontName', ip.Results.FontName, 'FontSize', ip.Results.FontSize);
        extent = get(ht, 'extent');
        textWidth(c) = extent(3);
        delete(ht);
    end
    maxTextWidth = max(textWidth);
    offset = maxTextWidth + 2*dci;
else
    offset = 0;
end
delete(ha);

% dt = track.lifetime_s/(track.end-track.start+1);
% gapIdx = round((track.t+dt)/dt);
% gapIdx = gapIdx(track.gapVect==1) - track.start + 1;

ha = zeros(nc,nf);
set(hf, 'Position', [50, 100, width+offset, height], 'Visible', ip.Results.Visible);%, 'ResizeFcn', {@resizeCallback});

% stack index si: x + (rowi-1)*nc*nx + (c-1)*nx
for si = 1:nf
    x = rem(si-1,nx)+1;
    rowi = ceil(si/nx); % each row contains all channels

    for c = 1:nc
        ha(c,si) = axes('Units', 'pixels',...
            'Position', [offset+(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
            'XLim', [0 wxi], 'YLim', [0 wxi]);
        imagesc(xa{si}, ya{si}, trackStack{c, si}); axis image off; caxis(dRange(c,:));
        hold on;

        % frame index, relative to movie
        fi = si-sb+track.start-1;

        if ip.Results.ShowDetection && track.start<=fi && fi<=track.end% && fi>sb && fi<=(track.end-track.start+1)+sb && c==1
            idx = find(track.f==fi);
            for i = idx
                %if track.gapVect(i)==1
                    %plot(track.x(c,i), track.y(c,i), 'o', 'Color', [1 1 1], 'MarkerSize', 12);
                %else
                if ~isnan(track.isPSF(c,i)) && track.isPSF(c,i)
                    plot(track.x(c,i), track.y(c,i), 'o', 'Color', [0 1 0], 'MarkerSize', 12);
                else
                    plot(track.x(c,i), track.y(c,i), 'o', 'Color', [1 0 0], 'MarkerSize', 12);
                end
            end
        end

        if ip.Results.ShowMarkers && track.start<=fi && fi<=track.end
            % start marker
            if fi==track.start
                if c==1
                    plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
                if c==nc && nc>1
                    plot(mean(xa{si}([1 end])), ya{si}(end), '^', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
            end

            % gaps
            if any(track.gapVect(track.f==fi)==1) && c==1
                plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 10, 'MarkerFaceColor', [0.8 0 0]);
            end

            % end marker
            if fi==track.end
                if c==1
                    plot(mean(xa{si}([1 end])), ya{si}(1), 'v', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
                if c==nc && nc>1
                    plot(mean(xa{si}([1 end])), ya{si}(end), '^', 'MarkerEdgeColor', 'none', 'Markersize', 8, 'MarkerFaceColor', [0 0 0]);
                end
            end
        end

        % Channel labels
        if x==1 && rowi==1 && ~isempty(labels)
            ht(c) = text(-dci, wxi/2, labels{c}, 'Units', 'pixels',...
                'HorizontalAlignment', 'right', 'VerticalAlignment', 'baseline',...
                'FontUnits', 'pixels', 'FontName', ip.Results.FontName, 'FontSize', wxi/2.5, 'Color', rgbColors{c});
        end
    end
end
colormap(gray(256));
end



function [wxi, dxi, dci, nr, height] = getProportions(width, nx, nf, nc)

% fixed proportions
wx = 1;
dx = 1/15; % gap between frames, relative to frame width
dc = 1/5; % gap between channels

% number of rows
nr = ceil(nf/nx);

iwidth = nx*wx + (nx-1)*dx;
iheight = nr*(nc*wx+(nc-1)*dx) + (nr-1)*dc;

height = width*iheight/iwidth;

wxi = width / ((nx/dx + nx-1)*dx);
dxi = wxi*dx;
dci = wxi*dc;

end


function [tstack, xa, ya] = getTrackStack(t, w, reference)

    sigma = frameInfo(1).s;
    w = ceil(w*sigma);

    % coordinate matrices
    x0 = tracks(t).x;
    y0 = tracks(t).y;

    % start and end buffer sizes
    if ~isempty(tracks(t).startBuffer)
        sb = numel(tracks(t).startBuffer.t);
        x0 = [tracks(t).startBuffer.x x0];
        y0 = [tracks(t).startBuffer.y y0];
    else
        sb = 0;
    end
    if ~isempty(tracks(t).endBuffer)
        eb = numel(tracks(t).endBuffer.t);
        x0 = [x0 tracks(t).endBuffer.x];
        y0 = [y0 tracks(t).endBuffer.y];
    else
        eb = 0;
    end

    % frame index
    tfi = tracks(t).start-sb:tracks(t).end+eb;
    tnf = length(tfi);


    if tracks(t).nSeg==1 && strcmpi(reference, 'track') % align frames to track
        xi = round(x0(handles.mCh,:));
        yi = round(y0(handles.mCh,:));
        % ensure that window falls within frame bounds
        x0 = xi - min([xi-1 w]);
        x1 = xi + min([nx-xi w]);
        y0 = yi - min([yi-1 w]);
        y1 = yi + min([ny-yi w]);
        % axes for each frame
        xa = arrayfun(@(i) x0(i):x1(i), 1:tnf, 'unif', 0);
        ya = arrayfun(@(i) y0(i):y1(i), 1:tnf, 'unif', 0);
    else
        % window around track mean
        mu_x = round(nanmean(x0,2));
        mu_y = round(nanmean(y0,2));
        x0 = max(1, min(mu_x)-w);
        x1 = min(data.imagesize(2), max(mu_x)+w);
        y0 = max(1, min(mu_y)-w);
        y1 = min(data.imagesize(1), max(mu_y)+w);
        xa = repmat({x0:x1}, [tnf 1]);
        ya = repmat({y0:y1}, [tnf 1]);
    end

    tstack = cell(nCh,tnf);
    for ci = 1:nCh
        for k = 1:tnf
            tstack{ci,k} = stack{ci}(ya{k}, xa{k}, tfi(k));
        end
    end
end
