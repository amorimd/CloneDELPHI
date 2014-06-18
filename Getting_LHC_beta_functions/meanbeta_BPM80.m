function [avbetaxBPM80,avbetayBPM80,lenBPM80,sbegBPM80,sendBPM80]=meanbeta_BPM80(twissfilename,beam)

% function to compute the average beta functions of combined BPMs CW in IR1 & 5 triplets region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

% BPMS: positions from R. Jones, HL-LHC PLC Jan 2013 (also in G. Arduini slides, 
% 9th WP2 task leaders' meeting, 06/06/2013)

ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% third BPM, stripline, after Q3, of diameter between electrodes 80mm
% positions from IPs:
sBPM80=[70.456];
sbegBPM80=sort([sBPM80 ringlength-sBPM80 sIP5-sBPM80 sIP5+sBPM80]) % 4 BPMs
sendBPM80=sbegBPM80+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM80=sum(sendBPM80-sbegBPM80);
[avbetaxBPM80,avbetayBPM80]=meanbetatwiss2(twissfilename,sbegBPM80,sendBPM80);

