function [fitted_fall fgof] = fallFit(track) %(TP) n should be track number in ProcessedTracks
  iv = track.A;
  y1 = [iv((find(iv==max(iv))):end)];
  x1 = [1:numel(y1)]*0.1;
  y1 = y1(find(~isnan(y1)));
  x1 = x1(find(~isnan(y1)));

  ft = fittype( 'power1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'Off';

  if numel(x1)<2
      fitted_fall = [];
      fgof = struct('rsquare', NaN);
  else
      [fitted_fall fgof] = fit( x1', y1', ft, opts);
  end
