function [avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 50A.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.

%%%%
% beam screen type 50 A, with vertical half gap b=36.75/2
% see N. Kos doc. "BS lengths for E_Metral v3.xls
%%%%

% dipoles beam screen
% around dipoles 
[sbegdip,senddip,names,len]=extract_beamscreens(twissfilename,'MB\.',15.495);
%size(sbegdip),names,sum(len),sum(senddip-sbegdip)
% 1232 beam screens, total BS length=1.908984000000095e+04


% quadrupoles beam screen
% around quadrupoles in FODO cells
[sbegQfodo,sendQfodo,names,len]=extract_beamscreens(twissfilename,'MQ\.(([23].)|(1[^01]))',5.8097);
%size(sbegQfodo),names,sum(len),sum(sendQfodo-sbegQfodo)
% 360 BS (= 362-2 in IR4 Q7, see N. Kos)
% total BS length=2.091491999999892e+03

% Q8 and Q10 in all IRs except IR3 and IR7
[sbegQ8_10,sendQ8_10,names,len]=extract_beamscreens(twissfilename,'MQ.*\.((8[RL][^37])|(10[RL][^37]))',7.0847);
%size(sbegQ8_10),names,sum(len),sum(sendQ8_10-sbegQ8_10)
% 24 BS. total BS length=1.700327999999911e+02

% Q9 in all IRs except IR3 and IR7 (two quads in each)
[sbegQ9,sendQ9,names,len]=extract_beamscreens(twissfilename,{'MQM\.9[RL][^37]','MQMC\.9[RL][^37]'},8.4847);
%size(sbegQ9),names,sum(len),sum(sendQ9-sbegQ9)
% 12 BS. total BS length=1.018164000000079e+02

% Q11 in IR1, 3, 5 and 7, right side (two quads in each)
[sbegQ11_IR1_3_5_7_R,sendQ11_IR1_3_5_7_R,names,len]=extract_beamscreens(twissfilename,{'MQ\.11R[1357]','MQTLI\.11R[1357]'},7.0847);
%size(sbegQ11_IR1_3_5_7_R),names,sum(len),sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)
% 4 BS. total BS length=28.338799999998287

% Q11 in IR1, 3, 5 and 7, left side (two quads in each)
[sbegQ11_IR1_3_5_7_L,sendQ11_IR1_3_5_7_L,names,len]=extract_beamscreens(twissfilename,{'MQ\.11L[1357]','MQTLI\.11L[1357]'},7.0847);
%size(sbegQ11_IR1_3_5_7_L),names,sum(len),sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)
% 4 BS. total BS length=28.338799999997718

% Q11 in IR2, 4, 6 and 8, right side (two quads in each)
[sbegQ11_IR2_4_6_8_R,sendQ11_IR2_4_6_8_R,names,len]=extract_beamscreens(twissfilename,{'MQ\.11R[2468]','MQTLI\.11R[2468]'},7.0847);
%size(sbegQ11_IR2_4_6_8_R),names,sum(len),sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)
% 4 BS. total BS length=28.338799999998628

% Q11 in IR2, 4, 6 and 8, left side (two quads in each)
[sbegQ11_IR2_4_6_8_L,sendQ11_IR2_4_6_8_L,names,len]=extract_beamscreens(twissfilename,{'MQ\.11L[2468]','MQTLI\.11L[2468]'},7.0847);
%size(sbegQ11_IR2_4_6_8_L),names,sum(len),sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)
% 4 BS. total BS length=28.338799999998628

% Q7, Q8 and Q10 in IR3 and IR7 (two quads in each)
[sbegQ7_8_10_IR3_7,sendQ7_8_10_IR3_7,namesQ7_8_10_IR3_7,len]=extract_beamscreens(twissfilename,{'MQ\.(([78][RL])|(10[RL]))','MQTLI\.(([78][RL])|(10[RL]))'},7.0847);
%size(sbegQ7_8_10_IR3_7),namesQ7_8_10_IR3_7,sum(len),sum(sendQ7_8_10_IR3_7-sbegQ7_8_10_IR3_7)
% 12 BS. total BS length=85.016399999993155

% Q9 in IR3 and IR7 (three quads in each)
[sbegQ9_IR3_7,sendQ9_IR3_7,names,len]=extract_beamscreens(twissfilename,{'MQ\.9[RL]','MQTLI\.A9[RL]','MQTLI\.B9[RL]'},8.4847);
%size(sbegQ9_IR3_7),names,sum(len),sum(sendQ9_IR3_7-sbegQ9_IR3_7)
% 4 BS. total BS length=33.938800000003539

% Q7 in IR4
[sbegQ7_IR4,sendQ7_IR4,namesQ7_IR4,len]=extract_beamscreens(twissfilename,'MQM\.7[RL]',5.8097);
%size(sbegQ7_IR4),namesQ7_IR4,sum(len),sum(sendQ7_IR4-sbegQ7_IR4)
% 2 BS. total BS length=11.619400000003225

% Q7 in IR1, 2 5 and 8 (two quads in each)
[sbegQ7_IR1_2_5_8,sendQ7_IR1_2_5_8,namesQ7_IR1_2_5_8,len]=extract_beamscreens(twissfilename,{'MQM\.A7[RL][1258]','MQM\.B7[RL][1258]'},9.4597);
%size(sbegQ7_IR1_2_5_8),namesQ7_IR1_2_5_8,sum(len),sum(sendQ7_IR1_2_5_8-sbegQ7_IR1_2_5_8)
% 8 BS. total BS length=75.677599999997824
% total number of BS in quads: 360+24+12+4*4+12+4+2+8=438 (as N. Kos)
% total lengthof BS in quads : 2.682948599999881e+03=
%sum(sendQfodo-sbegQfodo)+sum(sendQ8_10-sbegQ8_10)+sum(sendQ9-sbegQ9)+ ...
%    sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)+sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)+ ...
%    sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)+sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)+ ...
%    sum(sendQ7_8_10_IR3_7-sbegQ7_8_10_IR3_7)+sum(sendQ9_IR3_7-sbegQ9_IR3_7)+ ...
%    sum(sendQ7_IR4-sbegQ7_IR4)+sum(sendQ7_IR1_2_5_8-sbegQ7_IR1_2_5_8)


% drift type beam screens (still 50A)
% at Q11 (IP side), between MBB/MBA and Q11, in each IR (missing dipoles).
% we simply decrease sbegQ11_IR1_3_5_7_R and sbegQ11_IR2_4_6_8_R, and increase
% sendQ11_IR1_3_5_7_L and sendQ11_IR2_4_6_8_L, by an amount equal to the length
% of the beam screens
sbegQ11_IR1_3_5_7_R=sbegQ11_IR1_3_5_7_R-13.5517;
sbegQ11_IR2_4_6_8_R=sbegQ11_IR2_4_6_8_R-12.6097;
sendQ11_IR1_3_5_7_L=sendQ11_IR1_3_5_7_L+13.5517;
sendQ11_IR2_4_6_8_L=sendQ11_IR2_4_6_8_L+12.6097;
%sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)+sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)+...
%    sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)+sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)
% new total length: 3.226463999999984e+02

% LSS at the place of Q7 (missing) in IR6, left side (the side leading to IP5)
% we take it from the end of the MBA magnet in Q8 (left side), toward IP6
% MBA element name: MB.A8L6.B1 (1 Beam Screen only)
% We take as begining the end of the second drift DRIFT_342 after that MBA 
% (this drift is just after the sextupole MCS.A8L6.B1), to prevent some
% overlapping between the dipole beam screen and the LSS one.
if (beam==1) 
  sbegQ7L_IR6=16392.31082;
elseif (beam==2)
  sbegQ7L_IR6=16392.46318;
end
sendQ7L_IR6=sbegQ7L_IR6+9.955;
% 1 BS, length 9.955


% compiling all the type 50A beam screens
[sbeg50A,ind]=sort([sbegdip sbegQfodo sbegQ8_10 sbegQ9 sbegQ11_IR1_3_5_7_R ...
    sbegQ11_IR1_3_5_7_L sbegQ11_IR2_4_6_8_R sbegQ11_IR2_4_6_8_L...
    sbegQ7_8_10_IR3_7 sbegQ9_IR3_7 sbegQ7_IR4 sbegQ7_IR1_2_5_8 sbegQ7L_IR6]);
send50A=[senddip sendQfodo sendQ8_10 sendQ9 sendQ11_IR1_3_5_7_R ...
    sendQ11_IR1_3_5_7_L sendQ11_IR2_4_6_8_R sendQ11_IR2_4_6_8_L ...
    sendQ7_8_10_IR3_7 sendQ9_IR3_7 sendQ7_IR4 sendQ7_IR1_2_5_8 sendQ7L_IR6];
send50A=send50A(ind);
length(sbeg50A),length(send50A) % =1232+438+1 (16 drift BS are together with quad. ones)
len50A=sum(send50A-sbeg50A) % total length should be as for N. Kos: 2.19920348e+04
ringlength=26658.8832;
% overlapping test
if ( min(sbeg50A(2:end)-send50A(1:end-1))<0 ) || ( (send50A(end)-sbeg50A(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
%


% computes then the beta functions averages at all these s
% test first on dipoles only, using the dipole length instead of the BS length,
% comparing with meanbetatwiss.m
%[sbegtest,sendtest,names,len]=extract_beamscreens(twissfilename,'MB\.',14.3);
%[avbetaxtest,avbetaytest]=meanbetatwiss2(twissfilename,sbegtest,sendtest)
%[s,betax,betay,avbetax,avbetay,names,len]=meanbetatwiss(twissfilename,'MB.');
%abs(avbetax-avbetaxtest)/avbetax,abs(avbetay-avbetaytest)/avbetax
% almost no difference (~2e-11)



% now real calculation: 
% RESULT: average beta functions for the type 50A vertical BS
[avbetax50A,avbetay50A]=meanbetatwiss2(twissfilename,sbeg50A,send50A);
% avbetax50A= 86.2662 , avbetay50A= 91.085 (injection b1)
% avbetax50A= 86.2495 , avbetay50A= 91.1333 (collision b1)
% avbetax50A= 86.4865 , avbetay50A= 91.7239 (collision with sq 2m b1)


