function [avbetaxBPMCW,avbetayBPMCW,lenBPMCW,sbegBPMCW,sendBPMCW]=meanbeta_BPMCW_HLLHC(twissfilename,beam)

% function to compute the average beta functions of combined BPMs CW in 
% IR1 & 5 triplets region, for HL-LHC.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

ringlength=26658.8832;

% BPMS: positions from R. Jones, HL-LHC PLC Jan 2013 (also in G. Arduini slides, 
% 9th WP2 task leaders' meeting, 06/06/2013)

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (beam==1)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% first BPM (CW kind - combined stripline and buttons)
% positions from IPs:
sBPMCW=[21.474];
sbegBPMCW=sort([sBPMCW ringlength-sBPMCW sIP5-sBPMCW sIP5+sBPMCW]) % 4 BPMs
sendBPMCW=sbegBPMCW+0.3; % approximate length of 0.3m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPMCW=sum(sendBPMCW-sbegBPMCW);
[avbetaxBPMCW,avbetayBPMCW]=meanbetatwiss2(twissfilename,sbegBPMCW,sendBPMCW);
