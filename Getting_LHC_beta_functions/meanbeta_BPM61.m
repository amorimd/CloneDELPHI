function [avbetaxBPM61,avbetayBPM61,lenBPM61,sbegBPM61,sendBPM61]=meanbeta_BPM61(twissfilename,beam)

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

% second BPM, stripline, between Q1 and Q2, of diameter between electrodes 60mm
% positions from IPs:
sBPM61=[32];
sbegBPM61=sort([sBPM61 ringlength-sBPM61 sIP5-sBPM61 sIP5+sBPM61]) % 4 BPMs
sendBPM61=sbegBPM61+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM61=sum(sendBPM61-sbegBPM61);
[avbetaxBPM61,avbetayBPM61]=meanbetatwiss2(twissfilename,sbegBPM61,sendBPM61);

