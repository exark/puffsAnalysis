function isPuff = puffScorer(data, oldIsPuff, varargin)

if ~isempty(varargin)
    cmp = varargin{1};
    cmp = strcmp('compare',cmp);
    if cmp
        disp('You said compare, but that aint implemented yet');
    end
end

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


%-------------------------------------------------------------------------------
% Load detection masks
%-------------------------------------------------------------------------------
dpath = [data.source 'Detection' filesep 'Masks'];
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


if exist([data.source 'Detection' filesep 'cellmask.tif'], 'file')==2
    cellMask = imread([data.source 'Detection' filesep 'cellmask.tif']);
else
    cellMask = [];
end

%-------------------------------------------------------------------------------
% Load detection files
%-------------------------------------------------------------------------------

detectionFile = [data.channels{1} 'Detection' filesep 'detection_v2.mat'];
if (exist(detectionFile, 'file')==2)
    frameInfo = load(detectionFile);
    frameInfo = frameInfo.frameInfo;
else
    frameInfo = [];
end

%-------------------------------------------------------------------------------
% Load tracks
%-------------------------------------------------------------------------------

fprintf('Loading tracks ... ');
tracks = [];
bgA = [];
maxA = [];
% identify track file
fileList = dir([data.source 'Tracking' filesep 'ProcessedTracks*.mat']);
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


if exist([data.source 'Tracking' filesep fileName], 'file')==2
    tmp = load([data.source 'Tracking' filesep fileName]);
    tracks = tmp.tracks;
    if isfield(tmp, 'bgA')
        bgA = cellfun(@(i) prctile(i, 95, 2), tmp.bgA, 'unif', 0);
        bgA = [bgA{:}];
    end
    clear tmp;

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
        %tracks = tracks(cellMask(idx)==1);
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
assignin('base','newtracks',tracks);
fprintf('done.\n');

hfig = figure('Units', 'pixels', 'Position', [250 250 1200 800],...
    'PaperPositionMode', 'auto', 'Toolbar', 'none', 'Resize', 'on',...
    'CloseRequestFcn', @quitCallback,...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', 'PuffScorer');
pos = get(hfig, 'Position');

gcph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Scoring:',...
               'Position', [0 100 pos(3)/2 100]);
yesButton = uicontrol(gcph, 'Style', 'pushbutton', 'String', 'PUFF!',...
                      'Position', [50 25 100 50], 'FontSize', 26,...
                      'HorizontalAlignment','center',...
                      'Callback', @yesCallback);
noButton = uicontrol(gcph, 'Style', 'pushbutton', 'String', 'NONPUFF!',...
                      'Position', [220 25 180 50], 'FontSize', 26,...
                      'HorizontalAlignment','center',...
                      'Callback', @noCallback);

gph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Intensity Plot'],...
             'Position', [0 200 pos(3)/2 pos(4)-200]);
gax = axes('Parent', gph, 'Box', 'on', 'Units', 'pixels');


montph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Montage'],...
             'Position', [0 0 pos(3) 100]);
         
movieh = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Movie'],...
             'Position', [pos(3)/2 100 pos(3)/2 pos(4)-100]);
movieax = axes('Parent', movieh, 'Box', 'on', 'Units', 'pixels');

if isempty(oldIsPuff)
  isPuff = nan(numel(tracks),1);
else
  isPuff = oldIsPuff;
end
refreshTrack();

%======================================
% Helper functions:
%======================================
function quitCallback(varargin)
    cla(movieax, 'reset');
    assignin('base','isPuffNew',isPuff);
    if evalin('base','length(find(~isnan(isPuffNew))) > length(find(~isnan(isPuff)))')
        disp('Saving over old isPuff');
        assignin('base','isPuff',isPuff);
    else
        disp('isPuff in base workspace appears newer, saving as isPuffNew only.');
    end
    delete(gcf);
    return;
end

function refreshTrack()
    tcur = datasample(find(isnan(isPuff)), 1);
    cla(gax, 'reset');
    cla(movieax, 'reset');
    plotTrack(handles.data, tracks(tcur), 'Handle', gax);
    set(gph,'Title',['Intensity Plot Track: ' num2str(tcur)]);

    delete(get(montph,'Children'));
    [itrack, xa, ya] = getTrackStack(tcur, 10, 'track');
    set(montph,'Title',['Montage Track: ' num2str(tcur)]);
    set(movieh,'Title',['Movie Track ' num2str(tcur)]);
    plotTrackMontageLocal(tracks(tcur), itrack, xa, ya, montph, 1200, data.markers, movieax)
end

function yesCallback(varargin)
    isPuff = returnResults(tcur, 1, isPuff);
    refreshTrack();
end

function maybeCallback(varargin)
    isPuff = returnResults(tcur, 3, isPuff);
    refreshTrack();
end

function noCallback(varargin)
    isPuff = returnResults(tcur, 2, isPuff);
    refreshTrack();
end

function isPuffResult = returnResults(tn, res, isPuffBefore)
    isPuffBefore(tn) = res;
    isPuffResult = isPuffBefore;
    numPuffs = length(find(isPuffResult==1));
    disp(['Track ' num2str(tn) ' = ' num2str(res) ...
        ' (Scored ' num2str(length(find(~isnan(isPuffResult)))) '/' num2str(length(isPuffResult)) ...
        ', ' num2str(numPuffs) ' Puffs)']);
end

function plotTrackMontageLocal(track, trackStack, xa, ya, ph, width, labels, movieax)

[nc, nf] = size(trackStack);

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

dRange = zeros(nc,2);
c = 0;
for c = 1:nc
    cstack = cat(3, trackStack{c,:});
    dRange(c,:) = [min(cstack(:)) max(cstack(:))];
end

% display with fixed #frames/width
framesPerRow=40;
[wxi, dxi, dci, nr, height] = getProportions(width, framesPerRow, nf, nc);

% stack index si: x + (rowi-1)*nc*nx + (c-1)*nx
im_for_movie(nf+1) = struct('cdata', [], 'colormap', []);
for si = 1:nf
    x = rem(si-1,framesPerRow)+1;
    rowi = ceil(si/framesPerRow); % each row contains all channels

    for c = 1:nc
        ha(c,si) = axes('Units', 'pixels', 'Parent', ph,...
            'Position', [(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
            'XLim', [0 wxi], 'YLim', [0 wxi]);
        frame = imagesc(xa{si}, ya{si}, trackStack{c, si}); axis image off; caxis(dRange(c,:));
        frame = double(frame.CData-dRange(c,1))/double(dRange(c,2)-dRange(c,1));
        im_for_movie(si).cdata = imresize(uint8(round(frame*255)),20);
        im_for_movie(si).colormap = gray(256);
        hold on;
    end
end
colormap(gray(256));
im_for_movie(nf + 1) = struct('cdata', im_for_movie(1).cdata*0, 'colormap', gray(256));
set(movieax, 'Visible', 'Off')
movie(movieax, im_for_movie, 100, 5);
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

end
