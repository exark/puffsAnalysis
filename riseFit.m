function [fitted_rise rgof] = riseFit(track) %(TP) n should be track number in ProcessedTracks
  iv = track.A;
  x1 = [iv(1:(find(iv==max(iv))))];
  x1 = find(~isnan(x1));
  y1 = [1:numel(x1)]*0.1;
  y1 = y1(find(~isnan(x1)));

  ymax = mean(y1(1:numel(x1)));
  ybase = min(y1(1:numel(x1)));
  
  ymax = find(~isnan(ymax));
  ybase = ybase(find(~isnan(ymax)));

  ft = fittype( 'exp1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';
  opts.StartPoint = [ybase, ymax]
  
  if numel(x1)<2
      fitted_rise = [];
      rgof = struct('rsquare', NaN);
  else
      [fitted_rise rgof] = fit( x1', y1', ft, opts );
  end