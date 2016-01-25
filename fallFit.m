function [fitted_fall fgof] = fallFit(track) %(TP) n should be track number in ProcessedTracks
  iv = track.A;
  x1 = [iv((find(iv==max(iv))):end)];
  x1 = find(~isnan(x1));
  y1 = [1:numel(x1)]*0.1;
  y1 = y1(find(~isnan(x1)));

  ymax = mean(y1(numel(x1):end));
  ybase = min(y1(numel(x1):end));
  

  ft = fittype( 'power1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';
  opts.StartPoint = [ybase, ymax];

  if numel(x1)<2
      fitted_fall = [];
      fgof = struct('rsquare', 'NaN');
  else
      [fitted_fall fgof] = fit( x1', y1', ft, opts);
  end