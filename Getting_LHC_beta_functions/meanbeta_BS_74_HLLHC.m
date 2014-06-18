function [avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74_HLLHC(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 74 in HL-LHC (take out those in high-beta region close to IP1 & 5).
% in input: name of the twiss file to use for this. beam=1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.


%%%%
% beam screen type 74, with vertical half gap b=60.95/2
%%%%

% DFBX near Q3 in IR2 and IR8 (bet'n D1 and Q3)
% We take it right outside the quad of Q3 (MQXA)
% re-extract first some 63V and 63H beam screens
[sbegQ3_IR8,sendQ3_IR8,namesQ3_IR8,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL][8]',9.744);
[sbegQ3_IR2,sendQ3_IR2,namesQ3_IR2,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL][2]',9.744);

BSlengt=2.6235;IClengt=0.0;% BS length and interconnect length between BS (the
% latter is a rough estimation)
%namesQ3_IR1_2([2 3]),namesQ3_IR5_8([3 4])
sbegDFBX_Q3_IR2_8=[ [sbegQ3_IR2(1) sbegQ3_IR8(1)]-IClengt-BSlengt ...
    [sendQ3_IR2(2) sendQ3_IR8(2)]+IClengt];
sendDFBX_Q3_IR2_8=sbegDFBX_Q3_IR2_8+BSlengt;
%sum(sendDFBX_Q3_IR2_8-sbegDFBX_Q3_IR2_8)
% 4 BS. total length: 10.494

% D1 in IR2 and IR8 (in Q4)
[sbegD1_IR2_8,sendD1_IR2_8,names,len]=extract_beamscreens(twissfilename,'MBX\.4[RL][28]',10.8325);
%size(sbegD1_IR2_8),names,sum(len),sum(sendD1_IR2_8-sbegD1_IR2_8)
% 4 BS. total BS length=43.33

% compiling all the type 74 vertical beam screens
[sbeg74,ind]=sort([sbegDFBX_Q3_IR2_8 sbegD1_IR2_8]);
send74=[sendDFBX_Q3_IR2_8 sendD1_IR2_8];
send74=send74(ind);
length(sbeg74),length(send74) % 8
len74=sum(send74-sbeg74) % total length: 53.824
ringlength=26658.8832;
% overlapping test
if ( min(sbeg74(2:end)-send74(1:end-1))<0 ) || ( (send74(end)-sbeg74(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end


% RESULT: average beta functions for the type 74 vertical BS
[avbetax74,avbetay74]=meanbetatwiss2(twissfilename,sbeg74,send74);
% avbetax74=124.3366, avbetay74=124.3365 (injection b1)
% avbetax74=736.9396, avbetay74=736.9392 (collision sq 2m b1)

