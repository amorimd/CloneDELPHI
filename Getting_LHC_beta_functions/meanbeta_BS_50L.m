function [avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 50L.
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.

%%%%
% beam screen type 50L, with vertical half gap b=37.55/2
%%%%

% Q5 and Q6 in IR1 and IR5
[sbegQ5_6_IR1_5,sendQ5_6_IR1_5,names,len]=extract_beamscreens(twissfilename,'MQML\.[56][RL][15]',7.0847);
%size(sbegQ5_6_IR1_5),names,sum(len),sum(sendQ5_6_IR1_5-sbegQ5_6_IR1_5)
% 8 BS. total BS length=56.677599999996573

% Q6 in IR2 and IR8 (two quads in each)
[sbegQ6_IR2_8,sendQ6_IR2_8,names,len]=extract_beamscreens(twissfilename,{'MQM\.6[RL][28]','MQML\.6[RL][28]'},10.8647);
%size(sbegQ6_IR2_8),names,sum(len),sum(sendQ6_IR2_8-sbegQ6_IR2_8)
% 4 BS. total BS length=43.458799999996700

% Q6 in IR3 and IR7 (six quads in each)
[sbegQ6_IR3_7,sendQ6_IR3_7,names,len]=extract_beamscreens(twissfilename,{'MQTLH\.A6[RL][37]','MQTLH\.B6[RL][37]','MQTLH\.C6[RL][37]','MQTLH\.D6[RL][37]','MQTLH\.E6[RL][37]','MQTLH\.F6[RL][37]'},10.8647);
%size(sbegQ6_IR3_7),names,sum(len),sum(sendQ6_IR3_7-sbegQ6_IR3_7)
% 4 BS. total BS length=43.458799999996700

% Q5 in IR2 (right) and IR8 (left) (two quads in each)
[sbegQ5_IR2R_8L,sendQ5_IR2R_8L,names,len]=extract_beamscreens(twissfilename,{'MQM\.A5[RL][28]','MQM\.B5[RL][28]'},11.8197);
%size(sbegQ5_IR2R_8L),names,sum(len),sum(sendQ5_IR2R_8L-sbegQ5_IR2R_8L)
% 2 BS. total BS length=23.639400000000023

% total length: 167.2346 as for N. Kos
%sum(sendQ5_6_IR1_5-sbegQ5_6_IR1_5)+sum(sendQ6_IR2_8-sbegQ6_IR2_8)+sum(sendQ6_IR3_7-sbegQ6_IR3_7)+sum(sendQ5_IR2R_8L-sbegQ5_IR2R_8L)


% compiling all the type 50L beam screens
[sbeg50L,ind]=sort([sbegQ5_6_IR1_5 sbegQ6_IR2_8 sbegQ6_IR3_7 sbegQ5_IR2R_8L]);
send50L=[sendQ5_6_IR1_5 sendQ6_IR2_8 sendQ6_IR3_7 sendQ5_IR2R_8L];
send50L=send50L(ind);
length(sbeg50L),length(send50L) % =8+8+2
len50L=sum(send50L-sbeg50L) % total length should be as for N. Kos: 167.2346
ringlength=26658.8832;
% overlapping test
if ( min(sbeg50L(2:end)-send50L(1:end-1))<0 ) || ( (send50L(end)-sbeg50L(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% 

% RESULT: average beta functions for the type 50L vertical BS
[avbetax50L,avbetay50L]=meanbetatwiss2(twissfilename,sbeg50L,send50L);
% avbetax50L= 145.6699 , avbetay50L= 147.9059 (injection b1)
% avbetax50L= 165.471 , avbetay50L= 136.4996 (collision with sq 2m b1)
