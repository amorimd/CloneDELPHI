function [avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 53H.
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 53H, with horizontal half gap b=40.35/2
%%%%

% Q1 in IR1 and IR2
[sbegQ1_IR1_2,sendQ1_IR1_2,namesQ1_IR1_2,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][12]',8.6815-0.77-0.03);
%size(sbegQ1_IR1_2),namesQ1_IR1_2,sum(len),sum(sendQ1_IR1_2-sbegQ1_IR1_2)
% 4 BS. total BS length=31.526

% only those for this type of BS
sbeg53H=sbegQ1_IR1_2;send53H=sendQ1_IR1_2;
len53H=sum(send53H-sbeg53H) % total length should be as for N. Kos: 31.586

% RESULT: average beta functions for the type 53H vertical BS
[avbetax53H,avbetay53H]=meanbetatwiss2(twissfilename,sbeg53H,send53H);
% avbetax53H=77.0849, avbetay53H=77.0849 (injection b1)
% avbetax53H=349.9507, avbetay53H=349.9506 (collision with sq 2m b1)

