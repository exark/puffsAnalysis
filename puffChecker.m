function validList = puffChecker(data, unknownList, validList, varargin)

if ~isempty(varargin)
    cmp = varargin{1};
    cmp = strcmp('compare',cmp);
    if cmp
        disp('You said compare, but that aint implemented yet');
    end
end

if isempty(validList)
    validList = NaN(size(unknownList));
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
% Load cell mask
%-------------------------------------------------------------------------------
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
else
    msgID = 'puffChecker:NoTrackingFile';
    msg = 'Unable to find processedTracks.mat for specified data';
    baseException = MException(msgID,msg);
    throw(baseException);
end
fprintf('done.\n');

hfig = figure('Units', 'pixels', 'Position', [250 250 1200 600],...
    'PaperPositionMode', 'auto', 'Toolbar', 'none', 'Resize', 'on',...
    'CloseRequestFcn', @quitCallback,...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', 'PuffScorer');
pos = get(hfig, 'Position');

gcph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Scoring:',...
               'Position', [0 0 pos(3)/2 100]);
yesButton = uicontrol(gcph, 'Style', 'pushbutton', 'String', 'PUFF!',...
                      'Position', [50 25 100 50], 'FontSize', 26,...
                      'HorizontalAlignment','center',...
                      'Callback', @yesCallback);
maybeButton = uicontrol(gcph, 'Style', 'pushbutton', 'String', 'MAYBE!',...
                      'Position', [200 25 120 50], 'FontSize', 26,...
                      'HorizontalAlignment','center',...
                      'Callback', @maybeCallback);
noButton = uicontrol(gcph, 'Style', 'pushbutton', 'String', 'NONPUFF!',...
                      'Position', [370 25 180 50], 'FontSize', 26,...
                      'HorizontalAlignment','center',...
                      'Callback', @noCallback);

gph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Intensity Plot'],...
             'Position', [0 100 pos(3)/2 pos(4)-100]);
gax = axes('Parent', gph, 'Box', 'on', 'Units', 'pixels');


montph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Montage'],...
             'Position', [pos(3)/2 0 pos(3)/2 pos(4)]);

tidx = NaN;
refreshTrack();

%======================================
% Helper functions:
%======================================
function quitCallback(varargin)
    assignin('base','validListNew',validList);
    if evalin('base','length(find(~isnan(validListNew))) > length(find(~isnan(validList)))')
        disp('Saving over old validList');
        assignin('base','validList',validList);
    else
        disp('validList in base workspace appears newer, saving as validListNew only.');
    end
    delete(gcf);
    return;
end

function refreshTrack()
    try
        tidx = datasample(find(isnan(validList)), 1);
    catch excep
        disp('No more tracks! Time to relax, babe!');
        return
    end
    tcur = unknownList(tidx);
    cla(gax, 'reset');
    plotTrack(handles.data, tracks(tcur), 'Handle', gax);
    set(gph,'Title',['Intensity Plot Track: ' num2str(tcur)]);

    delete(get(montph,'Children'));
    [itrack, xa, ya] = getTrackStack(tcur, 6, 'track');
    plotTrackMontageLocal(tracks(tcur), itrack, xa, ya, montph, 600, data.markers)
    set(gph,'Title',['Montage Track: ' num2str(tcur)]);
end

function yesCallback(varargin)
    validList = returnResults(tidx, 1, validList);
    refreshTrack();
end

function maybeCallback(varargin)
    validList = returnResults(tidx, 3, validList);
    refreshTrack();
end

function noCallback(varargin)
    validList = returnResults(tidx, 2, validList);
    refreshTrack();
end

function validListResult = returnResults(idx, res, validListBefore)
    validListBefore(idx) = res;
    validListResult = validListBefore;
    disp(['Track ' num2str(tcur) ' = ' num2str(res) ' (Valid/Maybe/Invalid/Total: ' num2str(length(find(validListResult==1))) '/' num2str(length(find(validListResult==3))) '/' num2str(length(find(validListResult==2))) '/' num2str(length(validListResult))]);
end

function plotTrackMontageLocal(track, trackStack, xa, ya, ph, width, labels)

[nc, nf] = size(trackStack);

fontName = 'Helvetica';
fontSize = 14;
framesPerRow = 20;

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
[wxi, dxi, dci, nr, height] = getProportions(width, framesPerRow, nf, nc);

% stack index si: x + (rowi-1)*nc*nx + (c-1)*nx
for si = 1:nf
    x = rem(si-1,framesPerRow)+1;
    rowi = ceil(si/framesPerRow); % each row contains all channels

    for c = 1:nc
        ha(c,si) = axes('Units', 'pixels', 'Parent', ph,...
            'Position', [(x-1)*(wxi+dxi) height-wxi-((rowi-1)*(nc*wxi+(nc-1)*dxi+dci)+(c-1)*(wxi+dxi)) wxi wxi],...
            'XLim', [0 wxi], 'YLim', [0 wxi]);
        imagesc(xa{si}, ya{si}, trackStack{c, si}); axis image off; caxis(dRange(c,:));
        hold on;
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

end
