function [fitted_fall fgof numFall] = fallFit(track) 
  fitted_fall = cell( 2, 1 );
  fgof = struct( 'sse', cell(2,1), ...
    'rsquare', [], 'dfe', [], 'adjrsquare', [], 'rmse', [] );

  iv = track.A;
  y1 = [iv((find(iv==max(iv))):end)];
  x1 = [1:numel(y1)]*0.1;
  y1 = y1(find(~isnan(y1)));
  x1 = x1(find(~isnan(y1)));
  numFall = numel(y1);

  ftp = fittype( 'power1' );
  fte = fittype( 'exp1' );
  opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
  opts.Display = 'final';
  
  if numel(x1)<2
      fitted_fall = [];
      fgof(1).rsquare = -1;
      fgof(2).rsquare = -1;
  else
      [fitted_fall{1} fgof(1)] = fit( x1', y1', ftp, opts); % powerfit
      [fitted_fall{2} fgof(2)] = fit(x1', y1', fte, opts); %expfit
  end
