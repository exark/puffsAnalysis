function [fitted_rise rgof numRise] = riseFit(track) %(TP) the entire struct of the desired track
  iv = [track.A] + [track.c];
  [~,i] = max(iv)
  y1 = [iv(1:i)];
  x1 = [1:numel(y1)]*0.1;
  y1 = y1(find(~isnan(y1)));
  x1 = x1(find(~isnan(y1)));
  numRise = numel(y1);

  ft = fittype( 'exp1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares','Display','off');

  if numel(x1)<2
      fitted_rise = [];
      rgof = struct('rsquare', -1);
  else
      [fitted_rise rgof] = fit(x1', y1', ft, opts );
  end
