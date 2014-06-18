function [avbetaxtap74_WP,avbetaytap74_WP,lentap74_WP,sbegtap74_WP,sendtap74_WP]=meanbeta_taper74_WP(twissfilename,beam)

% function to compute the average beta functions of tapers between 74 
% and warm pipe type beam screens in IR1 & 5 triplets region.
% in input: name of the twiss file to use for this, beam=1 or 2
% in output: average beta functions in x and y, total length of the tapers
% and beginning and end (in terms of longitudinal position s) of each of them.

% Q3 in IR1
[sbegQ3_IR1,sendQ3_IR1,names,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL]1',9.744);
% Q3 in IR5
[sbegQ3_IR5,sendQ3_IR5,names,len]=extract_beamscreens(twissfilename,'MQXA\.3[RL]5',9.744);

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

% warm pipe in D1
% we take it right outside the BS 74 (until the end of D1) (length computed
% from plan)
WPlengt=27.25;IClengt=0.035;
sbegWP=[send74(1)+IClengt sbeg74(2)-WPlengt-IClengt send74(3)+IClengt sbeg74(4)-WPlengt-IClengt];
sendWP=sbegWP+WPlengt;

% tapers from BS 74 to warm pipe
sbegtap74_WP=[send74(1) sendWP(4) send74(3) sendWP(2)]
sendtap74_WP=[sbegWP(1) sbeg74(4) sbegWP(3) sbeg74(2)]
sendtap74_WP-sbegtap74_WP
lentap74_WP=sum(sendtap74_WP-sbegtap74_WP) % 4 tapers of 0.035m each
[avbetaxtap74_WP,avbetaytap74_WP]=meanbetatwiss2(twissfilename,sbegtap74_WP,sendtap74_WP);

