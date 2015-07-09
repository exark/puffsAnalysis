cmeAnalysis Workflow
====================

After executing loadConditionData(), all the magic starts with cmeAnalysis(). cmeAnalysis has essentially four steps:
* Detection: find point sources in the image
* Tracking: links point sources between frames
* Track Processing: Analyzes tracks for gaps and category.
* Lifetime analysis: Calculate statistics for CCP lifetimes

Below, each step is split into its component pieces.


## Track Processing
**runTrackProcessing.m** starts by reading in frameInfo from Detection. It proceeds by setting up a window for viewing pit lifetimes. This is completed by looking at the sigma generated during detection and expanding it to cover 4 standard deviations outside of that initial range. An *annularMask* is then generated which is used when the track is classified later. Next, the output of runTracking.m is read in.

Preprocessing (marked in the script) takes care of several manipulations required before actually classifying any of the tracks...

Processing cleans up each track, looking for gaps and invalidly detected tracks.

Post-processing uses the metadata generated from the processing step to categorize tracks into the following categories:

* Ia)  Single tracks with valid gaps
* Ib)  Single tracks with invalid gaps
* Ic)  Single tracks cut at beginning or end
* Id)  Single tracks, persistent
* IIa) Compound tracks with valid gaps
* IIb) Compound tracks with invalid gaps
* IIc) Compound tracks cut at beginning or end
* IId) Compound tracks, persistent
