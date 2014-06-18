function [avbetaxBBC_H,avbetayBBC_H,lenBBC_H,sbegBBC_H,sendBBC_H]=meanbeta_BBC_H(twissfilename,beam)

% function to compute the average beta functions of beam-beam horizontal wire compensator (BBC).
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the
% horizontal BBC and beginning and end (in terms of longitudinal position s) of each of them.

% we put it at 147m from the IP5 (on each side)

ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% 147m from IP
sBBC=[147];
sbegBBC_H=sort([sIP5-sBBC sIP5+sBBC]); % 2 BBC
sendBBC_H=sbegBBC_H+1; % approximate length of 1m
lenBBC_H=sum(sendBBC_H-sbegBBC_H);
[avbetaxBBC_H,avbetayBBC_H]=meanbetatwiss2(twissfilename,sbegBBC_H,sendBBC_H)

