function plotTrackAndMontage(data, tracks)

mCh = find(strcmp(data.source, data.channels));
nx = data.imagesize(2);
ny = data.imagesize(1);
nCh = 1;

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

hfig = figure('Units', 'pixels', 'Position', [250 250 1200 600],...
    'PaperPositionMode', 'auto', 'Toolbar', 'none', 'Resize', 'on',...
    'DefaultUicontrolUnits', 'pixels', 'Units', 'pixels', 'Name', 'PuffScorer');
pos = get(hfig, 'Position');

gcph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', 'Track #:',...
               'Position', [0 0 pos(3)/2 100]);
trackPicker = uicontrol(gcph, 'Style','edit', 'String', '1',...
    'Units','pixels','Position', [50 30 100 40], 'HorizontalAlignment','center', 'Callback', @refreshTrack);

gph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Intensity Plot'],...
             'Position', [0 100 pos(3)/2 pos(4)-100]);
gax = axes('Parent', gph, 'Box', 'on', 'Units', 'pixels');


montph = uipanel('Parent', hfig, 'Units', 'pixels', 'Title', ['Montage'],...
             'Position', [pos(3)/2 0 pos(3)/2 pos(4)]);

function refreshTrack(source, ~)
    tcur = str2num(source.String);
    cla(gax, 'reset');
    plotTrack(data, tracks(tcur), 'Handle', gax);
    set(gph,'Title',['Intensity Plot Track: ' num2str(tcur)]);

    delete(get(montph,'Children'));
    [itrack, xa, ya] = getTrackStack(tracks(tcur), 6, 'track');
    plotTrackMontageLocal(tracks(tcur), itrack, xa, ya, montph, 600, data.markers)
    set(gph,'Title',['Montage Track: ' num2str(tcur)]);
end

function [tstack, xa, ya] = getTrackStack(t, w, reference)

    w = ceil(w*1.1);

    % coordinate matrices
    x0 = t.x;
    y0 = t.y;

    % start and end buffer sizes
    if ~isempty(t.startBuffer)
        sb = numel(t.startBuffer.t);
        x0 = [t.startBuffer.x x0];
        y0 = [t.startBuffer.y y0];
    else
        sb = 0;
    end
    if ~isempty(t.endBuffer)
        eb = numel(t.endBuffer.t);
        x0 = [x0 t.endBuffer.x];
        y0 = [y0 t.endBuffer.y];
    else
        eb = 0;
    end

    % frame index
    tfi = t.start-sb:t.end+eb;
    tnf = length(tfi);


    if t.nSeg==1 && strcmpi(reference, 'track') % align frames to track
        xi = round(x0(mCh,:));
        yi = round(y0(mCh,:));
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

end