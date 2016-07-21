function [fitted_fall fgof numFall] = fallFit(track)
  iv = [track.A] + [track.c];
  [~,i] = max(iv);
  y1 = [iv(i:end)];
  x1 = [1:numel(y1)]*0.1;
  y1 = y1(find(~isnan(y1)));
  x1 = x1(find(~isnan(y1)));
  numFall = numel(y1);
  
  ftp = fittype( 'power1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares','Display','off');

  if numel(x1)<2
      fitted_fall = [];
      fgof = struct('rsquare', -1);
  else
      [fitted_fall fgof] = fit(x1', y1', ftp, opts); % powerfit
  end
