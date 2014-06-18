function [avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 53V.
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 53V, with vertical half gap b=40.35/2
%%%%

% Q1 in IR5 and IR8
[sbegQ1_IR5_8,sendQ1_IR5_8,namesQ1_IR5_8,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][58]',8.6815-0.77-0.03);
%size(sbegQ1_IR5_8),namesQ1_IR5_8,sum(len),sum(sendQ1_IR5_8-sbegQ1_IR5_8)
% 4 BS. total BS length=31.526

% only those for this type of BS
sbeg53V=sbegQ1_IR5_8;send53V=sendQ1_IR5_8;
len53V=sum(send53V-sbeg53V) % total length should be as for N. Kos: 31.586


% RESULT: average beta functions for the type 53V vertical BS
[avbetax53V,avbetay53V]=meanbetatwiss2(twissfilename,sbeg53V,send53V);
% avbetax53V=77.0849 , avbetay53V=77.0849 (injection b1)
% avbetax53V=349.8843 , avbetay53V=349.8842 (collision with sq 2m b1)

