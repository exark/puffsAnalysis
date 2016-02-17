function [fitted_rise rgof numRise] = riseFit(track) %(TP) n should be track number in ProcessedTracks
  iv = track.A;
  y1 = [iv(1:(find(iv==max(iv))))];
  x1 = [1:numel(y1)]*0.1;
  y1 = y1(find(~isnan(y1)));
  x1 = x1(find(~isnan(y1)));
  numRise = numel(y1);

  ft = fittype( 'exp1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'final';

  if numel(x1)<2
      fitted_rise = [];
      rgof = struct('rsquare', NaN);
  else
      [fitted_rise rgof] = fit(x1', y1', ft, opts );
  end
