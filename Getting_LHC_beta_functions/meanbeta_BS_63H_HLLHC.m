function [avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H_HLLHC(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 63H in HL-LHC (take out those in high-beta region close to IP1 & 5).
% in input: name of the twiss file to use for this.
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.



%%%%
% beam screen type 63H, with horizontal half gap b=50.45/2
%%%%

% extension of the 53 horizontal beam screen in Q1 (IR2)
% it is between Q1 and Q2
% re-extract first some 53H beam screens
[sbegQ1_IR2,sendQ1_IR2,namesQ1_IR2,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][2]',8.6815-0.77-0.03);
%size(sbegQ1_IR2),namesQ1_IR2,sum(len),sum(sendQ1_IR2-sbegQ1_IR2)
% 2 BS. total BS length=15.793
BSlengt=0.77;IClengt=0.03;
%namesQ1_IR2
sbegQ1_ext_IR2=[sbegQ1_IR2(2)-BSlengt-IClengt sendQ1_IR2(1)+IClengt];
sendQ1_ext_IR2=sbegQ1_ext_IR2+BSlengt;

% Q3 in IR2
[sbegQ3_IR2,sendQ3_IR2,namesQ3_IR2,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL][2]',9.744);
%size(sbegQ3_IR2),namesQ3_IR2,sum(len),sum(sendQ3_IR2-sbegQ3_IR2)
% 2 BS. total BS length=19.488

% Q2 in IR2 (two quads in each)
[sbegQ2_IR2,sendQ2_IR2,names,len]=extract_beamscreens(twissfilename,{'MQXB\.A2[RL][2]','MQXB\.B2[RL][2]'},13.23);
%size(sbegQ2_IR2),names,sum(len),sum(sendQ2_IR2-sbegQ2_IR2)
% 2 BS. total BS length=26.46

% compiling all the type 63 horizontal beam screens
[sbeg63H,ind]=sort([sbegQ1_ext_IR2 sbegQ3_IR2 sbegQ2_IR2]);
send63H=[sendQ1_ext_IR2 sendQ3_IR2 sendQ2_IR2];
send63H=send63H(ind);
length(sbeg63H),length(send63H) % 6
len63H=sum(send63H-sbeg63H) % total length: 47.488 (note: it seems that it does not correspond exactly to the calculation of N. Kos, probably it's OK)
ringlength=26658.8832;
% overlapping test
if ( min(sbeg63H(2:end)-send63H(1:end-1))<0 ) || ( (send63H(end)-sbeg63H(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end


% RESULT: average beta functions for the type 63 horizontal BS
[avbetax63H,avbetay63H]=meanbetatwiss2(twissfilename,sbeg63H,send63H);
