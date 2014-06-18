function [avbetaxtap63_74,avbetaytap63_74,lentap63_74,sbegtap63_74,sendtap63_74]=meanbeta_taper63_74(twissfilename,beam)

% function to compute the average beta functions of tapers between 63 
% and 74 type beam screens in IR1 & 5 triplets region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

% BS 53H: Q1 in IR1
[sbeg53H,send53H,names,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL]1',8.6815-0.77-0.03);

% BS 53V: Q1 in IR5
[sbeg53V,send53V,names,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL]5',8.6815-0.77-0.03);

% BS 63H
% extension of the 53 horizontal beam screen in Q1 (IR1)
% it is between Q1 and Q2
BSlengt=0.77;IClengt=0.03;
sbegQ1_ext_IR1=[sbeg53H(2)-BSlengt-IClengt send53H(1)+IClengt];
sendQ1_ext_IR1=sbegQ1_ext_IR1+BSlengt;
% Q3 in IR1
[sbegQ3_IR1,sendQ3_IR1,names,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL]1',9.744);
% Q2 in IR1 (two quads)
[sbegQ2_IR1,sendQ2_IR1,names,len]=extract_beamscreens(twissfilename,{'MQXB\.A2[RL]1','MQXB\.B2[RL]1'},13.23);
% compiling the type 63 horizontal beam screens
[sbeg63H,ind]=sort([sbegQ1_ext_IR1 sbegQ3_IR1 sbegQ2_IR1]);
send63H=[sendQ1_ext_IR1 sendQ3_IR1 sendQ2_IR1];
send63H=send63H(ind);

% BS 63V
% extension of the 53 vertical beam screen in Q1 (IR5)
% it is between Q1 and Q2
BSlengt=0.77;IClengt=0.03;
sbegQ1_ext_IR5=[sbeg53V(1)-BSlengt-IClengt send53V(2)+IClengt];
sendQ1_ext_IR5=sbegQ1_ext_IR5+BSlengt;
% Q3 in IR5
[sbegQ3_IR5,sendQ3_IR5,names,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL]5',9.744);
% Q2 in IR5 (two quads)
[sbegQ2_IR5,sendQ2_IR5,names,len]=extract_beamscreens(twissfilename,{'MQXB\.A2[RL]5','MQXB\.B2[RL]5'},13.23);
% compiling the type 63 vertical beam screens
[sbeg63V,ind]=sort([sbegQ1_ext_IR5 sbegQ3_IR5 sbegQ2_IR5]);
send63V=[sendQ1_ext_IR5 sendQ3_IR5 sendQ2_IR5];
send63V=send63V(ind);

% BS 74
% DFBX near Q3 in IR1 and IR5 (bet'n D1 and Q3)
% We take it right outside the quad of Q3 (MQXA)
BSlengt=2.6965;IClengt=0.077;% BS length and interconnect length between BS (the
% latter is the taper length)
sbegDFBX_Q3_IR1_5=[ [sbegQ3_IR1(2) sbegQ3_IR5(1)]-IClengt-BSlengt ...
    [sendQ3_IR1(1) sendQ3_IR5(2)]+IClengt];
sendDFBX_Q3_IR1_5=sbegDFBX_Q3_IR1_5+BSlengt;
% compiling the type 74 vertical beam screens
[sbeg74,ind]=sort([sbegDFBX_Q3_IR1_5]);
send74=[sendDFBX_Q3_IR1_5];
send74=send74(ind);


% tapers from BS 63 to BS 74 (note: tapers are round -> do not depend on BS
% orientation)
sbegtap63_74=[send63H(3) send74(4) send63V(6) send74(2)]
sendtap63_74=[sbeg74(1) sbeg63H(4) sbeg74(3) sbeg63V(1)]
sendtap63_74-sbegtap63_74
lentap63_74=sum(sendtap63_74-sbegtap63_74) % 4 tapers of 0.077m each
[avbetaxtap63_74,avbetaytap63_74]=meanbetatwiss2(twissfilename,sbegtap63_74,sendtap63_74);

