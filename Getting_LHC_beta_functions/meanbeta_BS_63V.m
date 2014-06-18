function [avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 63V.
% in input: name of the twiss file to use for this, beam = 1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 63V, with vertical half gap b=50.45/2
%%%%

% extension of the 53 vertical beam screen in Q1 (IR5 and 8)
% it is between Q1 and Q2
% re-extract first some 53V beam screens
[sbegQ1_IR5_8,sendQ1_IR5_8,namesQ1_IR5_8,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL][58]',8.6815-0.77-0.03);
%size(sbegQ1_IR5_8),namesQ1_IR5_8,sum(len),sum(sendQ1_IR5_8-sbegQ1_IR5_8)
% 4 BS. total BS length=31.586
BSlengt=0.77;IClengt=0.03;
%namesQ1_IR5_8
sbegQ1_ext_IR5_8=[sbegQ1_IR5_8([1 3])-BSlengt-IClengt sendQ1_IR5_8([2 4])+IClengt];
sendQ1_ext_IR5_8=sbegQ1_ext_IR5_8+BSlengt;
% 4 BS. total length=3.08

% undulator right (beam1) or left (beam2) in IR4
% it is next to the MBRS and indicated by MU in drawings and in madx (4 elements)
if (beam==1) 
[sbegUnd_IR4_R,sendUnd_IR4_R,names,len]=extract_beamscreens(twissfilename,{'MU\.A5R4','MU\.B5R4','MU\.C5R4','MU\.D5R4'},1.3405);
elseif (beam==2)
[sbegUnd_IR4_R,sendUnd_IR4_R,names,len]=extract_beamscreens(twissfilename,{'MU\.A5L4','MU\.B5L4','MU\.C5L4','MU\.D5L4'},1.3405);
end  
%size(sbegUnd_IR4_R),names,sum(len),sum(sendUnd_IR4_R-sbegUnd_IR4_R)
% 1 BS. total BS length=1.3405

% DFBA left near Q7 at all IRs
% it is right after (i.e. nearer the IP) the last quad in Q7, except for
% IR6: there it is right after the last beam screen at the missing dipole
% position.
% first: find back some beam screens of type 50A (we need their positions)
[sbegQ7_IR1_2_5_8,sendQ7_IR1_2_5_8,namesQ7_IR1_2_5_8,len]=extract_beamscreens(twissfilename,{'MQM\.A7[RL][1258]','MQM\.B7[RL][1258]'},9.4597);
[sbegQ7_8_10_IR3_7,sendQ7_8_10_IR3_7,namesQ7_8_10_IR3_7,len]=extract_beamscreens(twissfilename,{'MQ\.(([78][RL])|(10[RL]))','MQTLI\.(([78][RL])|(10[RL]))'},7.0847);
[sbegQ7_IR4,sendQ7_IR4,namesQ7_IR4,len]=extract_beamscreens(twissfilename,'MQM\.7[RL]',5.8097);
if (beam==1) 
  sbegQ7L_IR6=16392.31082;
elseif (beam==2)
  sbegQ7L_IR6=16392.46318;
end
sendQ7L_IR6=sbegQ7L_IR6+9.955;

BSlengt=1.7065;IClengt=0.165;% BS length and interconnect length between BS (the
% latter is an estimate; we took as for dipoles, see N. Kos)
%namesQ7_IR1_2_5_8(2:2:8),namesQ7_8_10_IR3_7([3 9]),namesQ7_IR4(1)
sbegDFBA_Q7L=[sendQ7_IR1_2_5_8(2:2:8) sendQ7_8_10_IR3_7([3 9]) sendQ7_IR4(1) ...
    sendQ7L_IR6] +IClengt;
sendDFBA_Q7L=sbegDFBA_Q7L+BSlengt;
%sum(sendDFBA_Q7L-sbegDFBA_Q7L)
% 8 BS. total BS length: 13.652

% DFBA right near Q7 at all IRs
% it is right before (i.e. nearer the IP) the last quad in Q7, except for
% IR6: there it is before the MBA of Q8, namely MB.A8R6.B1 (we take some margin:
% the end of the BS is at the beginning of DRIFT_670, before the two multipoles 
% MCO.8R6.B1 and MCD.8R6.B1. This is two prevent overlapping with the BS of the MBA).
BSlengt=2.208;IClengt=0.165;% BS length and interconnect length between BS (the
% latter is an estimation; we took as for dipoles, see N. Kos)
%namesQ7_IR1_2_5_8(1:2:7),namesQ7_8_10_IR3_7([4 10]),namesQ7_IR4(2)
if (beam==1)
  sendDFBA_Q7R=[sbegQ7_IR1_2_5_8(1:2:7) sbegQ7_8_10_IR3_7([4 10]) sbegQ7_IR4(2) ...
    16931.14082+IClengt] -IClengt;
elseif (beam==2)
  sendDFBA_Q7R=[sbegQ7_IR1_2_5_8(1:2:7) sbegQ7_8_10_IR3_7([4 10]) sbegQ7_IR4(2) ...
    16931.29318+IClengt] -IClengt;
end
sbegDFBA_Q7R=sendDFBA_Q7R-BSlengt;
%sum(sendDFBA_Q7R-sbegDFBA_Q7R)
% 8 BS. total BS length: 17.664

% Q4, Q5 in IR6 and Q5, Q6 in IR4
[sbegQ4_5_IR6_Q5_6_IR4,sendQ4_5_IR6_Q5_6_IR4,names,len]=extract_beamscreens(twissfilename,'MQY\.[456][RL][46]',6.1015);
%size(sbegQ4_5_IR6_Q5_6_IR4),names,sum(len),sum(sendQ4_5_IR6_Q5_6_IR4-sbegQ4_5_IR6_Q5_6_IR4)
% 8 BS. total BS length=48.811999999990803

% Q4 in IR1 and IR5
[sbegQ4_IR1_5,sendQ4_IR1_5,names,len]=extract_beamscreens(twissfilename,'MQY\.4[RL][15]',8.7765);
%size(sbegQ4_IR1_5),names,sum(len),sum(sendQ4_IR1_5-sbegQ4_IR1_5)
% 4 BS. total BS length=35.105999999999824

% Q3 in IR5 and IR8
[sbegQ3_IR5_8,sendQ3_IR5_8,namesQ3_IR5_8,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL][58]',9.744);
%size(sbegQ3_IR5_8),namesQ3_IR5_8,sum(len),sum(sendQ3_IR5_8-sbegQ3_IR5_8)
% 4 BS. total BS length=38.975999999995111

% MQY at Q4 and Q5 in IR2 and IR8 (for Q5, only left in IR2 and right in
% IR8) (two quads at each Q)
[sbegQ4_5_IR2_8,sendQ4_5_IR2_8,names,len]=extract_beamscreens(twissfilename,{'MQY\.A[45][RL][28]','MQY\.B[45][RL][28]'},12.1115);
%size(sbegQ4_5_IR2_8),names,sum(len),sum(sendQ4_5_IR2_8-sbegQ4_5_IR2_8)
% 6 BS. total BS length=72.668999999997141

% Q2 in IR5 and IR8 (two quads in each)
[sbegQ2_IR5_8,sendQ2_IR5_8,names,len]=extract_beamscreens(twissfilename,{'MQXB\.A2[RL][58]','MQXB\.B2[RL][58]'},13.23);
%size(sbegQ2_IR5_8),names,sum(len),sum(sendQ2_IR5_8-sbegQ2_IR5_8)
% 4 BS. total BS length=52.920

% compiling all the type 63 vertical beam screens
[sbeg63V,ind]=sort([sbegQ1_ext_IR5_8 sbegUnd_IR4_R sbegDFBA_Q7L sbegDFBA_Q7R sbegQ4_5_IR6_Q5_6_IR4 ...
    sbegQ4_IR1_5 sbegQ3_IR5_8 sbegQ4_5_IR2_8 sbegQ2_IR5_8]);
send63V=[sendQ1_ext_IR5_8 sendUnd_IR4_R sendDFBA_Q7L sendDFBA_Q7R sendQ4_5_IR6_Q5_6_IR4 ...
    sendQ4_IR1_5 sendQ3_IR5_8 sendQ4_5_IR2_8 sendQ2_IR5_8];
send63V=send63V(ind);
length(sbeg63V),length(send63V) % 47=4+1+2*8+8+4+4+6+4
len63V=sum(send63V-sbeg63V) % total length: 284.2195 (note: it seems that it does not correspond exactly to the calculation of N. Kos, probably it's OK)
ringlength=26658.8832;
% overlapping test
if ( min(sbeg63V(2:end)-send63V(1:end-1))<0 ) || ( (send63V(end)-sbeg63V(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end



% RESULT: average beta functions for the type 63 vertical BS
[avbetax63V,avbetay63V]=meanbetatwiss2(twissfilename,sbeg63V,send63V);
% avbetax63V=159.5394, avbetay63V=160.1524 (injection b1)
% avbetax63V=412.2242, avbetay63V=418.7521 (collision sq 2m b1)

