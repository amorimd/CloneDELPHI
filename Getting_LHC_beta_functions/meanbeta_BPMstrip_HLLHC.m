function [avbetaxBPMs,avbetayBPMs,lenBPMs,sbegBPMs,sendBPMs]=meanbeta_BPMstrip_HLLHC(twissfilename,beam)

% function to compute the average beta functions of stripline BPMs in 
% IR1 & 5 triplets region, for HL-LHC.
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

% other BPMs (stripline)
% positions from IPs:
sBPMs=[32 42.667 53.263 62.737 70.456 81.193];
sbegBPMs=sort([sBPMs ringlength-sBPMs sIP5-sBPMs sIP5+sBPMs]) % 24 BPMs
sendBPMs=sbegBPMs+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPMs=sum(sendBPMs-sbegBPMs);
[avbetaxBPMs,avbetayBPMs]=meanbetatwiss2(twissfilename,sbegBPMs,sendBPMs);
