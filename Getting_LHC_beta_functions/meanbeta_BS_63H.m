function [avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 63H.
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.



%%%%
% beam screen type 63H, with horizontal half gap b=50.45/2
%%%%

% extension of the 53 horizontal beam screen in Q1 (IR1 and 2)
% it is between Q1 and Q2
% re-extract first some 53H beam screens
[sbegQ1_IR1_2,sendQ1_IR1_2,namesQ1_IR1_2,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][12]',8.6815-0.77-0.03);
%size(sbegQ1_IR1_2),namesQ1_IR1_2,sum(len),sum(sendQ1_IR1_2-sbegQ1_IR1_2)
% 4 BS. total BS length=31.586
BSlengt=0.77;IClengt=0.03;
%namesQ1_IR1_2
sbegQ1_ext_IR1_2=[sbegQ1_IR1_2([2 4])-BSlengt-IClengt sendQ1_IR1_2([1 3])+IClengt];
sendQ1_ext_IR1_2=sbegQ1_ext_IR1_2+BSlengt;
% 4 BS. total length=3.08

% Q3 in IR1 and IR2
[sbegQ3_IR1_2,sendQ3_IR1_2,namesQ3_IR1_2,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL][12]',9.744);
%size(sbegQ3_IR1_2),namesQ3_IR1_2,sum(len),sum(sendQ3_IR1_2-sbegQ3_IR1_2)
% 4 BS. total BS length=38.975999999998152

% Q2 in IR1 and IR2 (two quads in each)
[sbegQ2_IR1_2,sendQ2_IR1_2,names,len]=extract_beamscreens(twissfilename,{'MQXB\.A2[RL][12]','MQXB\.B2[RL][12]'},13.23);
%size(sbegQ2_IR1_2),names,sum(len),sum(sendQ2_IR1_2-sbegQ2_IR1_2)
% 4 BS. total BS length=52.920000000002332

% compiling all the type 63 horizontal beam screens
[sbeg63H,ind]=sort([sbegQ1_ext_IR1_2 sbegQ3_IR1_2 sbegQ2_IR1_2]);
send63H=[sendQ1_ext_IR1_2 sendQ3_IR1_2 sendQ2_IR1_2];
send63H=send63H(ind);
length(sbeg63H),length(send63H) % 12
len63H=sum(send63H-sbeg63H) % total length: 94.976 (note: it seems that it does not correspond exactly to the calculation of N. Kos, probably it's OK)
ringlength=26658.8832;
% overlapping test
if ( min(sbeg63H(2:end)-send63H(1:end-1))<0 ) || ( (send63H(end)-sbeg63H(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end


% RESULT: average beta functions for the type 63 horizontal BS
[avbetax63H,avbetay63H]=meanbetatwiss2(twissfilename,sbeg63H,send63H);
% avbetax63H=148.5693, avbetay63H=148.5692 (injection b1)
% avbetax63H=745.2177, avbetay63H=745.2174 (collision sq 2m b1)
