function [avbetaxBBC_V,avbetayBBC_V,lenBBC_V,sbegBBC_V,sendBBC_V]=meanbeta_BBC_V(twissfilename,beam)

% function to compute the average beta functions of vertical beam-beam wire compensator (BBC).
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the vertical BBC
% and beginning and end (in terms of longitudinal position s) of each of them.

% we put it at 147m from IP1 (on each side)

ringlength=26658.8832;

% 147m from IP
sBBC=[147];
sbegBBC_V=sort([sBBC ringlength-sBBC]); % 2 BBC
sendBBC_V=sbegBBC_V+1; % approximate length of 1m
lenBBC_V=sum(sendBBC_V-sbegBBC_V);
[avbetaxBBC_V,avbetayBBC_V]=meanbetatwiss2(twissfilename,sbegBBC_V,sendBBC_V)

