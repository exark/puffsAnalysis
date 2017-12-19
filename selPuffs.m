function [trainStruct, heldOutStruct] = selPuffs(scoredTracks, trainStruct, heldOutStruct)

b2cell12puffs = find([scoredTracks.isPuff]==1);
b2cell12selpuffs = datasample(b2cell12puffs, round(numel(b2cell12puffs)/2), 'Replace',false);
b2cell12heldpuffs = setdiff(b2cell12puffs,b2cell12selpuffs);

b2cell12nonpuffs = find([scoredTracks.isPuff]==2);
b2cell12selnonpuffs = datasample(b2cell12nonpuffs, round(numel(b2cell12nonpuffs)/2), 'Replace',false);
b2cell12heldnonpuffs = setdiff(b2cell12nonpuffs,b2cell12selnonpuffs);

trainStruct = [trainStruct scoredTracks(b2cell12selpuffs) scoredTracks(b2cell12selnonpuffs)];
heldOutStruct = [heldOutStruct scoredTracks(b2cell12heldpuffs) scoredTracks(b2cell12heldnonpuffs)];