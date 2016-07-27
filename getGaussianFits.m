function res = getGaussianFits(tcur, tracks, data)

[tstack, ~, ~] = getTrackStackStandalone(tcur, 6, 'track', 1.1, tracks, data);
cur_track = tracks(tcur);

for i = 1:length(tstack)
    
    [~, ~, ~, res(i)] = fitGaussian2D(double(tstack{i}), [0 0 cur_track.A(i) 1.1 cur_track.c(i)], 'As');
    
end

end

