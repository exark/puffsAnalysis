function [tstack, xa, ya] = getTrackStackStandalone(t, w, reference, sigma ,tracks, data)

    %bullshit placeholder stuff
    handles.mCh = 1;
    nx = 512;
    ny = 512;
    nCh = 1;
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

    w = ceil(w*sigma);

    % coordinate matrices
    x0 = tracks(t).x;
    y0 = tracks(t).y;

    % start and end buffer sizes
   sb=0;
   eb=0;

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