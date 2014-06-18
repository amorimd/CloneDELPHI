function [avbetaxBPM61CW,avbetayBPM61CW,lenBPM61CW,sbegBPM61CW,sendBPM61CW]=meanbeta_BPM61CW(twissfilename,beam)

% function to compute the average beta functions of combined BPMs CW in IR1 & 5 triplets region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

% BPMS: positions from R. Jones, HL-LHC PLC Jan 2013 (also in G. Arduini slides, 
% 9th WP2 task leaders' meeting, 06/06/2013)

ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (strfind(twissfilename,'b1')>0)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% BPM before Q1 (BPM of CW kind - combined button + stripline), of diameter between electrodes 61mm
% positions from IPs:
sBPM61CW=[21.474];
sbegBPM61CW=sort([sBPM61CW ringlength-sBPM61CW sIP5-sBPM61CW sIP5+sBPM61CW]) % 4 BPMs
sendBPM61CW=sbegBPM61CW+0.3; % approximate length of 0.3m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM61CW=sum(sendBPM61CW-sbegBPM61CW);
[avbetaxBPM61CW,avbetayBPM61CW]=meanbetatwiss2(twissfilename,sbegBPM61CW,sendBPM61CW);
