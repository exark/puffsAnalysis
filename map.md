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

Preprocessing (marked in the script) takes care of several manipulations required before actually classifying any of the tracks
