function [avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V_HLLHC(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 53V in HL-LHC (take out those in high-beta region close to IP1 & 5).
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 53V, with vertical half gap b=40.35/2
%%%%

% Q1 in IR8
[sbegQ1_IR8,sendQ1_IR8,namesQ1_IR8,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][8]',8.6815-0.77-0.03);
%size(sbegQ1_IR8),namesQ1_IR8,sum(len),sum(sendQ1_IR8-sbegQ1_IR8)
% 2 BS. total BS length=15.793

% only those for this type of BS
sbeg53V=sbegQ1_IR8;send53V=sendQ1_IR8;
len53V=sum(send53V-sbeg53V)


% RESULT: average beta functions for the type 53V vertical BS
[avbetax53V,avbetay53V]=meanbetatwiss2(twissfilename,sbeg53V,send53V);

