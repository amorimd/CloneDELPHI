%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% various components beta functions (beam screen, warm pipe, 
% MBW, MQW) at 3.5TeV collision with squeeze of 2m
% Energy is not the point, the important thing is the squeeze to 2m 
% (it is the only significant difference between coll3.5tev.2m and injection)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% beam screen beta functions at collision
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% NOTE: nomenclature in madx (and in drawings) (see also Elias annotations on 
% N. Kos documents):
% MQML.10R6.... : MQML is the magnet or element name, 10 the number of the
% Q near which it is (there can be a letter just before that, if there are
% several same elements in the same Q), R/L for right/left of the IP, and 6
% is the number of interaction region considered (IR1 to IR8)

% We try first to extract and compute the beginnings and ends of
% all the beam screen section, using their approximate position and their lengths
% provided by N. Kos (see BS lengths for_E_Metral v3.xls).
% Beware, the real total length of beam screen is given by sum(send-sbeg),
% not by sum(len) which is the total length of the elements (dipoles,
% quads, etc.)
% NOTE: when there are several quads per Q with the same beam screen, two
% or more columns appear in "names", one needs to check in "names" that in each line we
% are at the same Q, same IR, and same side of the IR (right or left)

%%%%
% beam screen type 50 A, with vertical half gap b=36.75/2
% see N. Kos doc. "BS lengths for E_Metral v3.xls
%%%%

% dipoles beam screen
% around dipoles 
[sbegdip,senddip,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MB\.',15.495);
size(sbegdip),names,sum(len),sum(senddip-sbegdip)
% 1232 beam screens, total BS length=1.908984000000095e+04


% quadrupoles beam screen
% around quadrupoles in FODO cells
[sbegQfodo,sendQfodo,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQ\.(([23].)|(1[^01]))',5.8097);
size(sbegQfodo),names,sum(len),sum(sendQfodo-sbegQfodo)
% 360 BS (= 362-2 in IR4 Q7, see N. Kos)
% total BS length=2.091491999999892e+03

% Q8 and Q10 in all IRs except IR3 and IR7
[sbegQ8_10,sendQ8_10,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQ.*\.((8[RL][^37])|(10[RL][^37]))',7.0847);
size(sbegQ8_10),names,sum(len),sum(sendQ8_10-sbegQ8_10)
% 24 BS. total BS length=1.700327999999911e+02

% Q9 in all IRs except IR3 and IR7 (two quads in each)
[sbegQ9,sendQ9,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQM\.9[RL][^37]','MQMC\.9[RL][^37]'},8.4847);
size(sbegQ9),names,sum(len),sum(sendQ9-sbegQ9)
% 12 BS. total BS length=1.018164000000079e+02

% Q11 in IR1, 3, 5 and 7, right side (two quads in each)
[sbegQ11_IR1_3_5_7_R,sendQ11_IR1_3_5_7_R,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.11R[1357]','MQTLI\.11R[1357]'},7.0847);
size(sbegQ11_IR1_3_5_7_R),names,sum(len),sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)
% 4 BS. total BS length=28.338799999998287

% Q11 in IR1, 3, 5 and 7, left side (two quads in each)
[sbegQ11_IR1_3_5_7_L,sendQ11_IR1_3_5_7_L,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.11L[1357]','MQTLI\.11L[1357]'},7.0847);
size(sbegQ11_IR1_3_5_7_L),names,sum(len),sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)
% 4 BS. total BS length=28.338799999997718

% Q11 in IR2, 4, 6 and 8, right side (two quads in each)
[sbegQ11_IR2_4_6_8_R,sendQ11_IR2_4_6_8_R,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.11R[2468]','MQTLI\.11R[2468]'},7.0847);
size(sbegQ11_IR2_4_6_8_R),names,sum(len),sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)
% 4 BS. total BS length=28.338799999998628

% Q11 in IR2, 4, 6 and 8, left side (two quads in each)
[sbegQ11_IR2_4_6_8_L,sendQ11_IR2_4_6_8_L,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.11L[2468]','MQTLI\.11L[2468]'},7.0847);
size(sbegQ11_IR2_4_6_8_L),names,sum(len),sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)
% 4 BS. total BS length=28.338799999998628

% Q7, Q8 and Q10 in IR3 and IR7 (two quads in each)
[sbegQ7_8_10_IR3_7,sendQ7_8_10_IR3_7,namesQ7_8_10_IR3_7,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.(([78][RL])|(10[RL]))','MQTLI\.(([78][RL])|(10[RL]))'},7.0847);
size(sbegQ7_8_10_IR3_7),namesQ7_8_10_IR3_7,sum(len),sum(sendQ7_8_10_IR3_7-sbegQ7_8_10_IR3_7)
% 12 BS. total BS length=85.016399999993155

% Q9 in IR3 and IR7 (three quads in each)
[sbegQ9_IR3_7,sendQ9_IR3_7,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQ\.9[RL]','MQTLI\.A9[RL]','MQTLI\.B9[RL]'},8.4847);
size(sbegQ9_IR3_7),names,sum(len),sum(sendQ9_IR3_7-sbegQ9_IR3_7)
% 4 BS. total BS length=33.938800000003539

% Q7 in IR4
[sbegQ7_IR4,sendQ7_IR4,namesQ7_IR4,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQM\.7[RL]',5.8097);
size(sbegQ7_IR4),namesQ7_IR4,sum(len),sum(sendQ7_IR4-sbegQ7_IR4)
% 2 BS. total BS length=11.619400000003225

% Q7 in IR1, 2 5 and 8 (two quads in each)
[sbegQ7_IR1_2_5_8,sendQ7_IR1_2_5_8,namesQ7_IR1_2_5_8,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQM\.A7[RL][1258]','MQM\.B7[RL][1258]'},9.4597);
size(sbegQ7_IR1_2_5_8),namesQ7_IR1_2_5_8,sum(len),sum(sendQ7_IR1_2_5_8-sbegQ7_IR1_2_5_8)
% 8 BS. total BS length=75.677599999997824
% total number of BS in quads: 360+24+12+4*4+12+4+2+8=438 (as N. Kos)
% total lengthof BS in quads : 2.682948599999881e+03=
sum(sendQfodo-sbegQfodo)+sum(sendQ8_10-sbegQ8_10)+sum(sendQ9-sbegQ9)+ ...
    sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)+sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)+ ...
    sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)+sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)+ ...
    sum(sendQ7_8_10_IR3_7-sbegQ7_8_10_IR3_7)+sum(sendQ9_IR3_7-sbegQ9_IR3_7)+ ...
    sum(sendQ7_IR4-sbegQ7_IR4)+sum(sendQ7_IR1_2_5_8-sbegQ7_IR1_2_5_8)


% drift type beam screens (still 50A)
% at Q11 (IP side), between MBB/MBA and Q11, in each IR (missing dipoles).
% we simply decrease sbegQ11_IR1_3_5_7_R and sbegQ11_IR2_4_6_8_R, and increase
% sendQ11_IR1_3_5_7_L and sendQ11_IR2_4_6_8_L, by an amount equal to the length
% of the beam screens
sbegQ11_IR1_3_5_7_R=sbegQ11_IR1_3_5_7_R-13.5517
sbegQ11_IR2_4_6_8_R=sbegQ11_IR2_4_6_8_R-12.6097
sendQ11_IR1_3_5_7_L=sendQ11_IR1_3_5_7_L+13.5517
sendQ11_IR2_4_6_8_L=sendQ11_IR2_4_6_8_L+12.6097
sum(sendQ11_IR2_4_6_8_R-sbegQ11_IR2_4_6_8_R)+sum(sendQ11_IR1_3_5_7_R-sbegQ11_IR1_3_5_7_R)+...
    sum(sendQ11_IR2_4_6_8_L-sbegQ11_IR2_4_6_8_L)+sum(sendQ11_IR1_3_5_7_L-sbegQ11_IR1_3_5_7_L)
% new total length: 3.226463999999984e+02

% LSS at the place of Q7 (missing) in IR6, left side (the side leading to IP5)
% we take it from the end of the MBA magnet in Q8 (left side), toward IP6
% MBA element name: MB.A8L6.B1 (1 Beam Screen only)
% We take as begining the end of the second drift DRIFT_342 after that MBA 
% (this drift is just after the sextupole MCS.A8L6.B1), to prevent some
% overlapping between the dipole beam screen and the LSS one.
sbegQ7L_IR6=16392.31082,sendQ7L_IR6=sbegQ7L_IR6+9.955
% 1 BS, length 9.955


% compiling all the type 50A beam screens
[sbeg50A,ind]=sort([sbegdip sbegQfodo sbegQ8_10 sbegQ9 sbegQ11_IR1_3_5_7_R ...
    sbegQ11_IR1_3_5_7_L sbegQ11_IR2_4_6_8_R sbegQ11_IR2_4_6_8_L...
    sbegQ7_8_10_IR3_7 sbegQ9_IR3_7 sbegQ7_IR4 sbegQ7_IR1_2_5_8 sbegQ7L_IR6]);
send50A=[senddip sendQfodo sendQ8_10 sendQ9 sendQ11_IR1_3_5_7_R ...
    sendQ11_IR1_3_5_7_L sendQ11_IR2_4_6_8_R sendQ11_IR2_4_6_8_L ...
    sendQ7_8_10_IR3_7 sendQ9_IR3_7 sendQ7_IR4 sendQ7_IR1_2_5_8 sendQ7L_IR6];
send50A=send50A(ind);
length(sbeg50A),length(send50A) % 1671=1232+438+1 (16 drift BS are together with quad. ones)
sum(send50A-sbeg50A) % total length as for N. Kos: 2.19920348e+04
ringlength=26658.8832;
% overlapping test
if ( min(sbeg50A(2:end)-send50A(1:end-1))<0 ) || ( (send50A(end)-sbeg50A(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% success !

% computes then the beta functions averages at all these s
% test first on dipoles only, using the dipole length instead of the BS length,
% comparing with meanbetatwiss.m
%[sbegtest,sendtest,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MB\.',14.3);
%[avbetaxtest,avbetaytest]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbegtest,sendtest)
%[s,betax,betay,avbetax,avbetay,names,len]=meanbetatwiss('twiss.asbuilt.b1.coll3.5tev.2m.thick','MB.');
%abs(avbetax-avbetaxtest)/avbetax,abs(avbetay-avbetaytest)/avbetax
% almost no difference (~2e-11)



% now real calculation: 
% RESULT: average beta functions for the type 50A vertical BS
[avbetax50A,avbetay50A]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg50A,send50A)
% avbetax50A= 86.4865 , avbetay50A= 91.7239



%%%%
% beam screen type 50L, with vertical half gap b=37.55/2
%%%%

% Q5 and Q6 in IR1 and IR5
[sbegQ5_6_IR1_5,sendQ5_6_IR1_5,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQML\.[56][RL][15]',7.0847);
size(sbegQ5_6_IR1_5),names,sum(len),sum(sendQ5_6_IR1_5-sbegQ5_6_IR1_5)
% 8 BS. total BS length=56.677599999996573

% Q6 in IR2 and IR8 (two quads in each)
[sbegQ6_IR2_8,sendQ6_IR2_8,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQM\.6[RL][28]','MQML\.6[RL][28]'},10.8647);
size(sbegQ6_IR2_8),names,sum(len),sum(sendQ6_IR2_8-sbegQ6_IR2_8)
% 4 BS. total BS length=43.458799999996700

% Q6 in IR3 and IR7 (six quads in each)
[sbegQ6_IR3_7,sendQ6_IR3_7,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQTLH\.A6[RL][37]','MQTLH\.B6[RL][37]','MQTLH\.C6[RL][37]','MQTLH\.D6[RL][37]','MQTLH\.E6[RL][37]','MQTLH\.F6[RL][37]'},10.8647);
size(sbegQ6_IR3_7),names,sum(len),sum(sendQ6_IR3_7-sbegQ6_IR3_7)
% 4 BS. total BS length=43.458799999996700

% Q5 in IR2 (right) and IR8 (left) (two quads in each)
[sbegQ5_IR2R_8L,sendQ5_IR2R_8L,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQM\.A5[RL][28]','MQM\.B5[RL][28]'},11.8197);
size(sbegQ5_IR2R_8L),names,sum(len),sum(sendQ5_IR2R_8L-sbegQ5_IR2R_8L)
% 2 BS. total BS length=23.639400000000023

% total length: 167.2346 as for N. Kos
sum(sendQ5_6_IR1_5-sbegQ5_6_IR1_5)+sum(sendQ6_IR2_8-sbegQ6_IR2_8)+sum(sendQ6_IR3_7-sbegQ6_IR3_7)+sum(sendQ5_IR2R_8L-sbegQ5_IR2R_8L)


% compiling all the type 50L beam screens
[sbeg50L,ind]=sort([sbegQ5_6_IR1_5 sbegQ6_IR2_8 sbegQ6_IR3_7 sbegQ5_IR2R_8L]);
send50L=[sendQ5_6_IR1_5 sendQ6_IR2_8 sendQ6_IR3_7 sendQ5_IR2R_8L];
send50L=send50L(ind);
length(sbeg50L),length(send50L) % =8+8+2
sum(send50L-sbeg50L) % total length as for N. Kos: 167.2346
ringlength=26658.8832;
% overlapping test
if ( min(sbeg50L(2:end)-send50L(1:end-1))<0 ) || ( (send50L(end)-sbeg50L(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% OK

% compiling all the beam screens up to now, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50L sbeg50A]);
sendtot=[sendtot send50L send50A];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 50L vertical BS
[avbetax50L,avbetay50L]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg50L,send50L)
% avbetax50L= 165.471 , avbetay50L= 136.4996


%%%%
% beam screen type 53, with vertical half gap b=40.35/2
%%%%

% Q1 in IR5 and IR8
[sbegQ1_IR5_8,sendQ1_IR5_8,namesQ1_IR5_8,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQXA\.1[RL][58]',8.6815-0.77-0.015);
size(sbegQ1_IR5_8),namesQ1_IR5_8,sum(len),sum(sendQ1_IR5_8-sbegQ1_IR5_8)
% 4 BS. total BS length=31.586

% only those for this type of BS
sbeg53V=sbegQ1_IR5_8;send53V=sendQ1_IR5_8;

% compiling all the beam screens up to now, and operlapping test
[sbegtot,ind]=sort([sbegtot sbeg53V]);
sendtot=[sendtot send53V];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 53V vertical BS
[avbetax53V,avbetay53V]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg53V,send53V)
% avbetax53V=349.8843 , avbetay53V=349.8842

%%%%
% beam screen type 53, with horizontal half gap b=40.35/2
%%%%

% Q1 in IR1 and IR2
[sbegQ1_IR1_2,sendQ1_IR1_2,namesQ1_IR1_2,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQXA\.1[RL][12]',8.6815-0.77-0.015);
size(sbegQ1_IR1_2),namesQ1_IR1_2,sum(len),sum(sendQ1_IR1_2-sbegQ1_IR1_2)
% 4 BS. total BS length=31.586

% only those for this type of BS
sbeg53H=sbegQ1_IR1_2;send53H=sendQ1_IR1_2;

% compiling all the beam screens up to now, and operlapping test
[sbegtot,ind]=sort([sbegtot sbeg53H]);
sendtot=[sendtot send53H];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 53H vertical BS
[avbetax53H,avbetay53H]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg53H,send53H)
% avbetax53H=349.9507, avbetay53H=349.9506


%%%%
% beam screen type 63, with vertical half gap b=50.45/2
%%%%

% extension of the 53 vertical beam screen in Q1 (IR5 and 8)
% it is between Q1 and Q2
BSlengt=0.77
namesQ1_IR5_8
sbegQ1_ext_IR5_8=[sbegQ1_IR5_8([1 3])-BSlengt sendQ1_IR5_8([2 4])]
sendQ1_ext_IR5_8=sbegQ1_ext_IR5_8+BSlengt
% 4 BS. total length=3.08

% undulator right in IR4
% it is next to the MBRS and indicated by MU in drawings and in madx (4 elements)
[sbegUnd_IR4_R,sendUnd_IR4_R,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MU\.A5R4','MU\.B5R4','MU\.C5R4','MU\.D5R4'},1.3405);
size(sbegUnd_IR4_R),names,sum(len),sum(sendUnd_IR4_R-sbegUnd_IR4_R)
% 1 BS. total BS length=1.3405

% DFBA left near Q7 at all IRs
% it is right after (i.e. nearer the IP) the last quad in Q7, except for
% IR6: there it is right after the last beam screen at the missing dipole
% position.
BSlengt=1.7065;IClengt=0.165;% BS length and interconnect length between BS (the
% latter is an estimation; we took as for dipoles, see N. Kos)
namesQ7_IR1_2_5_8(2:2:8),namesQ7_8_10_IR3_7([3 9]),namesQ7_IR4(1)
sbegDFBA_Q7L=[sendQ7_IR1_2_5_8(2:2:8) sendQ7_8_10_IR3_7([3 9]) sendQ7_IR4(1) ...
    sendQ7L_IR6] +IClengt
sendDFBA_Q7L=sbegDFBA_Q7L+BSlengt
sum(sendDFBA_Q7L-sbegDFBA_Q7L)
% 8 BS. total BS length: 13.652

% DFBA right near Q7 at all IRs
% it is right before (i.e. nearer the IP) the last quad in Q7, except for
% IR6: there it is before the MBA of Q8, namely MB.A8R6.B1 (we take some margin:
% the end of the BS is at the beginning of DRIFT_670, before the two multipoles 
% MCO.8R6.B1 and MCD.8R6.B1. This is two prevent overlapping with the BS of the MBA).
BSlengt=2.208;IClengt=0.165;% BS length and interconnect length between BS (the
% latter is an estimation; we took as for dipoles, see N. Kos)
namesQ7_IR1_2_5_8(1:2:7),namesQ7_8_10_IR3_7([4 10]),namesQ7_IR4(2)
sendDFBA_Q7R=[sbegQ7_IR1_2_5_8(1:2:7) sbegQ7_8_10_IR3_7([4 10]) sbegQ7_IR4(2) ...
    16931.14082+IClengt] -IClengt
sbegDFBA_Q7R=sendDFBA_Q7R-BSlengt
sum(sendDFBA_Q7R-sbegDFBA_Q7R)
% 8 BS. total BS length: 17.664

% Q4, Q5 in IR6 and Q5, Q6 in IR4
[sbegQ4_5_IR6_Q5_6_IR4,sendQ4_5_IR6_Q5_6_IR4,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQY\.[456][RL][46]',6.1015);
size(sbegQ4_5_IR6_Q5_6_IR4),names,sum(len),sum(sendQ4_5_IR6_Q5_6_IR4-sbegQ4_5_IR6_Q5_6_IR4)
% 8 BS. total BS length=48.811999999990803

% Q4 in IR1 and IR5
[sbegQ4_IR1_5,sendQ4_IR1_5,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQY\.4[RL][15]',8.7765);
size(sbegQ4_IR1_5),names,sum(len),sum(sendQ4_IR1_5-sbegQ4_IR1_5)
% 4 BS. total BS length=35.105999999999824

% Q3 in IR5 and IR8
[sbegQ3_IR5_8,sendQ3_IR5_8,namesQ3_IR5_8,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQXA\.3[RL][58]',9.744);
size(sbegQ3_IR5_8),namesQ3_IR5_8,sum(len),sum(sendQ3_IR5_8-sbegQ3_IR5_8)
% 4 BS. total BS length=38.975999999995111

% MQY at Q4 and Q5 in IR2 and IR8 (for Q5, only left in IR2 and right in
% IR8) (two quads at each Q)
[sbegQ4_5_IR2_8,sendQ4_5_IR2_8,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQY\.A[45][RL][28]','MQY\.B[45][RL][28]'},12.1115);
size(sbegQ4_5_IR2_8),names,sum(len),sum(sendQ4_5_IR2_8-sbegQ4_5_IR2_8)
% 6 BS. total BS length=72.668999999997141

% Q2 in IR5 and IR8 (two quads in each)
[sbegQ2_IR5_8,sendQ2_IR5_8,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQXB\.A2[RL][58]','MQXB\.B2[RL][58]'},13.23);
size(sbegQ2_IR5_8),names,sum(len),sum(sendQ2_IR5_8-sbegQ2_IR5_8)
% 4 BS. total BS length=52.920

% compiling all the type 63 vertical beam screens
[sbeg63V,ind]=sort([sbegQ1_ext_IR5_8 sbegUnd_IR4_R sbegDFBA_Q7L sbegDFBA_Q7R sbegQ4_5_IR6_Q5_6_IR4 ...
    sbegQ4_IR1_5 sbegQ3_IR5_8 sbegQ4_5_IR2_8 sbegQ2_IR5_8]);
send63V=[sendQ1_ext_IR5_8 sendUnd_IR4_R sendDFBA_Q7L sendDFBA_Q7R sendQ4_5_IR6_Q5_6_IR4 ...
    sendQ4_IR1_5 sendQ3_IR5_8 sendQ4_5_IR2_8 sendQ2_IR5_8];
send63V=send63V(ind);
length(sbeg63V),length(send63V) % 47=4+1+2*8+8+4+4+6+4
sum(send63V-sbeg63V) % total length: 284.2195
ringlength=26658.8832;
% overlapping test
if ( min(sbeg63V(2:end)-send63V(1:end-1))<0 ) || ( (send63V(end)-sbeg63V(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% OK


% compiling all the beam screens up to now, and operlapping test
[sbegtot,ind]=sort([sbegtot sbeg63V]);
sendtot=[sendtot send63V];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 63 vertical BS
[avbetax63V,avbetay63V]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg63V,send63V)
% avbetax63V=412.2242, avbetay63V=418.7521



%%%%
% beam screen type 63, with horizontal half gap b=50.45/2
%%%%

% extension of the 53 horizontal beam screen in Q1 (IR1 and 2)
% it is between Q1 and Q2
BSlengt=0.77
namesQ1_IR1_2
sbegQ1_ext_IR1_2=[sbegQ1_IR1_2([2 4])-BSlengt sendQ1_IR1_2([1 3])]
sendQ1_ext_IR1_2=sbegQ1_ext_IR1_2+BSlengt
% 4 BS. total length=3.08

% Q3 in IR1 and IR2
[sbegQ3_IR1_2,sendQ3_IR1_2,namesQ3_IR1_2,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQXA\.3[RL][12]',9.744);
size(sbegQ3_IR1_2),namesQ3_IR1_2,sum(len),sum(sendQ3_IR1_2-sbegQ3_IR1_2)
% 4 BS. total BS length=38.975999999998152

% Q2 in IR1 and IR2 (two quads in each)
[sbegQ2_IR1_2,sendQ2_IR1_2,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick',{'MQXB\.A2[RL][12]','MQXB\.B2[RL][12]'},13.23);
size(sbegQ2_IR1_2),names,sum(len),sum(sendQ2_IR1_2-sbegQ2_IR1_2)
% 4 BS. total BS length=52.920000000002332

% compiling all the type 63 horizontal beam screens
[sbeg63H,ind]=sort([sbegQ1_ext_IR1_2 sbegQ3_IR1_2 sbegQ2_IR1_2]);
send63H=[sendQ1_ext_IR1_2 sendQ3_IR1_2 sendQ2_IR1_2];
send63H=send63H(ind);
length(sbeg63H),length(send63H) % 12
sum(send63H-sbeg63H) % total length: 94.976
ringlength=26658.8832;
% overlapping test
if ( min(sbeg63H(2:end)-send63H(1:end-1))<0 ) || ( (send63H(end)-sbeg63H(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% OK


% compiling all the beam screens up to now, and operlapping test
[sbegtot,ind]=sort([sbegtot sbeg63H]);
sendtot=[sendtot send63H];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 63 horizontal BS
[avbetax63H,avbetay63H]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg63H,send63H)
% avbetax63H=745.2177, avbetay63H=745.2174


%%%%
% beam screen type 69, with vertical half gap b=56.15/2
%%%%

% D2 in IR1, IR2, IR5 and IR8 (in Q4) and D4 in IR4 (in Q5)
[sbegD2_4_IR1_2_4_5_8,sendD2_4_IR1_2_4_5_8,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MBR((C\.4[RL][1258])|(B\.5[RL][4]))',10.7325);
size(sbegD2_4_IR1_2_4_5_8),names,sum(len),sum(sendD2_4_IR1_2_4_5_8-sbegD2_4_IR1_2_4_5_8)
% 10 BS. total BS length=1.073249999999879e+02

% D3 in IR4 (in Q5)
[sbegD3_IR4,sendD3_IR4,namesD3_IR4,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MBRS\.5[RL]4',10.9895);
size(sbegD3_IR4),namesD3_IR4,sum(len),sum(sendD3_IR4-sbegD3_IR4)
% we offset the right one a bit toward the right to avoid overlapping 
% with the beam screen of the undulator
sbegD3_IR4(2)=sbegD3_IR4(2)+0.3,sendD3_IR4(2)=sendD3_IR4(2)+0.3
% 2 BS. total BS length=21.978999999999360

% undulator left in IR4
% it is next to the MBRS of D3 (going toward IP4)
BSlengt=1.3405,IClengt=0.165 % interconnect length is a rough estimation
namesD3_IR4
sbegUnd_IR4_L=sendD3_IR4(1)+IClengt
sendUnd_IR4_L=sbegUnd_IR4_L+BSlengt
size(sbegUnd_IR4_L),sum(sendUnd_IR4_L-sbegUnd_IR4_L)
% 1 BS. total BS length=1.3405

% compiling all the type 69 vertical beam screens
[sbeg69,ind]=sort([sbegD2_4_IR1_2_4_5_8 sbegD3_IR4 sbegUnd_IR4_L]);
send69=[sendD2_4_IR1_2_4_5_8 sendD3_IR4 sendUnd_IR4_L];
send69=send69(ind);
length(sbeg69),length(send69) % 13
sum(send69-sbeg69) % total length: 130.6445 
ringlength=26658.8832;
% overlapping test
if ( min(sbeg69(2:end)-send69(1:end-1))<0 ) || ( (send69(end)-sbeg69(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% OK


% compiling all the beam screens up to now, and operlapping test
[sbegtot,ind]=sort([sbegtot sbeg69]);
sendtot=[sendtot send69];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length up to now
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% super


% RESULT: average beta functions for the type 69 vertical BS
[avbetax69,avbetay69]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg69,send69)
% avbetax69=295.2808, avbetay69=275.0131


%%%%
% beam screen type 74, with vertical half gap b=60.95/2
%%%%

% DFBX near Q3 in IR2 and IR8 (bet'n D1 and Q3)
% We take it right outside the quad of Q3 (MQXA)
BSlengt=2.6235;IClengt=0.0;% BS length and interconnect length between BS (the
% latter is a rough estimation)
namesQ3_IR1_2([2 3]),namesQ3_IR5_8([3 4])
sbegDFBX_Q3_IR2_8=[ [sbegQ3_IR1_2(2) sbegQ3_IR5_8(3)]-IClengt-BSlengt ...
    [sendQ3_IR1_2(3) sendQ3_IR5_8(4)]+IClengt]
sendDFBX_Q3_IR2_8=sbegDFBX_Q3_IR2_8+BSlengt
sum(sendDFBX_Q3_IR2_8-sbegDFBX_Q3_IR2_8)
% 4 BS. total length: 10.494

% DFBX near Q3 in IR1 and IR5 (bet'n D1 and Q3)
% We take it right outside the quad of Q3 (MQXA)
BSlengt=2.6965;IClengt=0.165;% BS length and interconnect length between BS (the
% latter is a rough estimation)
namesQ3_IR1_2([1 4]),namesQ3_IR5_8([1 2])
sbegDFBX_Q3_IR1_5=[ [sbegQ3_IR1_2(4) sbegQ3_IR5_8(1)]-IClengt-BSlengt ...
    [sendQ3_IR1_2(1) sendQ3_IR5_8(2)]+IClengt]
sendDFBX_Q3_IR1_5=sbegDFBX_Q3_IR1_5+BSlengt
sum(sendDFBX_Q3_IR1_5-sbegDFBX_Q3_IR1_5)
% 4 BS. total length: 10.786

% D1 in IR2 and IR8 (in Q4)
[sbegD1_IR2_8,sendD1_IR2_8,names,len]=extract_beamscreens('twiss.asbuilt.b1.coll3.5tev.2m.thick','MBX\.4[RL][28]',10.8325);
size(sbegD1_IR2_8),names,sum(len),sum(sendD1_IR2_8-sbegD1_IR2_8)
% 4 BS. total BS length=43.33

% compiling all the type 74 vertical beam screens
[sbeg74,ind]=sort([sbegDFBX_Q3_IR2_8 sbegDFBX_Q3_IR1_5 sbegD1_IR2_8]);
send74=[sendDFBX_Q3_IR2_8 sendDFBX_Q3_IR1_5 sendD1_IR2_8];
send74=send74(ind);
length(sbeg74),length(send74) % 12
sum(send74-sbeg74) % total length: 64.61
ringlength=26658.8832;
% overlapping test
if ( min(sbeg74(2:end)-send74(1:end-1))<0 ) || ( (send74(end)-sbeg74(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% OK


% RESULT: average beta functions for the type 74 vertical BS
[avbetax74,avbetay74]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbeg74,send74)
% avbetax74=736.9396, avbetay74=736.9392



%%%%%%%%%%%%%%%%%%%
% compiling all the beam screens, and overlapping test
%%%%%%%%%%%%%%%%%%%
[sbegtot,ind]=sort([sbegtot sbeg74]);
sendtot=[sendtot send74];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) % 1781=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % total length : 22796.8914
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end
% yeah !
sbegtot(1),sendtot(end),ringlength



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% summary of the data on beam screens
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% type 50A
halfgap50A=36.75/2;lengt50A=2.19920348e+04;avbetax50A=86.4865;avbetay50A=91.7239;
% type 50L
halfgap50L=37.55/2;lengt50L=167.2346;avbetax50L=165.471;avbetay50L=136.4996;
% type 53V
halfgap53V=40.35/2;lengt53V=31.586;avbetax53V=349.8843;avbetay53V=349.8842;
% type 53H
halfgap53H=40.35/2;lengt53H=31.586;avbetax53H=349.9507;avbetay53H=349.9506;
% type 63V
halfgap63V=50.45/2;lengt63V=284.2195;avbetax63V=412.2242;avbetay63V=418.7521;
% type 63H
halfgap63H=50.45/2;lengt63H=94.976;avbetax63H=745.2177;avbetay63H=745.2174;
% type 69
halfgap69=56.15/2;lengt69=130.6445;avbetax69=295.2808;avbetay69=275.0131;
% type 74
halfgap74=60.95/2;lengt74=64.61;avbetax74=736.9396;avbetay74=736.9392;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MBW and MQW (copper, warm) at coll3.5tev with 2m squeeze
% From LHC design report, chap.5:
% MBW: circular, radius 22mm
% MQW: elliptic, half gap 14.5mm, (a-b)/(a+b)=0.275 (half are vertical and
% the other half horizontal)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% MBW
[s,betax,betay,avbetaxMBW,avbetayMBW,namesMBW,lenMBW,sbegMBW,sendMBW]=meanbetatwiss('twiss.asbuilt.b1.coll3.5tev.2m.thick','MBW');
length(namesMBW),sum(lenMBW),sum(sendMBW-sbegMBW)
% 21 elements. 
% avbetaxMBW=140.045 and avbetayMBW=136.6282 (in m)
% total length: 70.6 (m). Elias had 70m (see Mathematica notebook,
% 17-11-2009)
% small check that same average beta:
%[avbetax,avbetay]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbegMBW,sendMBW)
% same result

% MQW
[s,betax,betay,avbetaxMQW,avbetayMQW,namesMQW,lenMQW,sbegMQW,sendMQW]=meanbetatwiss('twiss.asbuilt.b1.coll3.5tev.2m.thick','MQW');
length(namesMQW),sum(lenMQW),sum(sendMQW-sbegMQW)
% 48 elements. 
% avbetaxMQW=192.2401 and avbetayMQW=210.154 (in m)
% total length: 149.184 (m). Elias had 155m (see Mathematica notebook,
% 17-11-2009)


% compiling with all the beam screens, and overlapping test
[sbegtot,ind]=sort([sbegtot sbegMBW sbegMQW]);
sendtot=[sendtot sendMBW sendMQW];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) % 1850=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8+21+48 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % total length : 23016.6754
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping');
end
% oh yeah !
sbegtot(1),sendtot(end),ringlength



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rest of the machine (where warm pipe) at coll3.5tev with 2m squeeze
% Note: we should take out the collimators from this, but we won't
% (pessismistic approach)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
sbegrest=[0 sendtot];
sendrest=[sbegtot ringlength];
sum(sendrest-sbegrest) % total length: 3642.2078
% tests
sum(sendrest-sbegrest)+sum(sendtot-sbegtot)-ringlength % should be zero
if ( min(sbegrest(2:end)-sendrest(1:end-1))<0 ) || ( (sendrest(end)-sbegrest(1)-ringlength)>0 )
    disp('Aie... Overlapping');
end
% RESULT: average beta functions for the warm pipe (rest of the machine)
[avbetaxrest,avbetayrest]=meanbetatwiss2('twiss.asbuilt.b1.coll3.5tev.2m.thick',sbegrest,sendrest)
% avbetaxrest=217.0115, avbetayrest=216.0646


%%%%
% small check on the average beta functions
%%%%

[s,betax,betay,avbetax,avbetay]=meanbetatwiss('twiss.asbuilt.b1.coll3.5tev.2m.thick');

avbetax2=(avbetax50A*sum(send50A-sbeg50A)+avbetax50L*sum(send50L-sbeg50L)+ ...
    avbetax53V*sum(send53V-sbeg53V)+avbetax53H*sum(send53H-sbeg53H)+ ...
    avbetax63V*sum(send63V-sbeg63V)+avbetax63H*sum(send63H-sbeg63H)+ ...
    avbetax69*sum(send69-sbeg69)+avbetax74*sum(send74-sbeg74)+ ...
    avbetaxMBW*sum(sendMBW-sbegMBW)+avbetaxMQW*sum(sendMQW-sbegMQW)+ ...
    avbetaxrest*sum(sendrest-sbegrest))/ringlength
avbetay2=(avbetay50A*sum(send50A-sbeg50A)+avbetay50L*sum(send50L-sbeg50L)+ ...
    avbetay53V*sum(send53V-sbeg53V)+avbetay53H*sum(send53H-sbeg53H)+ ...
    avbetay63V*sum(send63V-sbeg63V)+avbetay63H*sum(send63H-sbeg63H)+ ...
    avbetay69*sum(send69-sbeg69)+avbetay74*sum(send74-sbeg74)+ ...
    avbetayMBW*sum(sendMBW-sbegMBW)+avbetayMQW*sum(sendMQW-sbegMQW)+ ...
    avbetayrest*sum(sendrest-sbegrest))/ringlength

(avbetax2-avbetax)/avbetax
(avbetay2-avbetay)/avbetay
% perfect (-1e-11 relative difference)


