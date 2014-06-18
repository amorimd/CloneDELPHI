function [avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H_HLLHC(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 53H in HL-LHC (take out those in high-beta region close to IP1 & 5).
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 53H, with horizontal half gap b=40.35/2
%%%%

% Q1 in IR2
[sbegQ1_IR2,sendQ1_IR2,namesQ1_IR2,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][2]',8.6815-0.77-0.03);
%size(sbegQ1_IR2),namesQ1_IR2,sum(len),sum(sendQ1_IR2-sbegQ1_IR2)
% 2 BS. total BS length=15.793

% only those for this type of BS
sbeg53H=sbegQ1_IR2;send53H=sendQ1_IR2;
len53H=sum(send53H-sbeg53H)

% RESULT: average beta functions for the type 53H vertical BS
[avbetax53H,avbetay53H]=meanbetatwiss2(twissfilename,sbeg53H,send53H);

