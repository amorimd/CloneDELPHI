function [avbetaxtap53_63,avbetaytap53_63,lentap53_63,sbegtap53_63,sendtap53_63]=meanbeta_taper53_63(twissfilename,beam)

% function to compute the average beta functions of all tapers between 53 
% and 63 type beam screens in IR1 & 5 triplets region.
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


% tapers from BS 53 to BS 63 (note: tapers are round -> do not depend on BS
% orientation)
sbegtap53_63=[send53H(1) send63H(6) send53V(2) send63V(3)]
sendtap53_63=[sbeg63H(1) sbeg53H(2) sbeg63V(4) sbeg53V(1)]
sendtap53_63-sbegtap53_63
lentap53_63=sum(sendtap53_63-sbegtap53_63) % 4 tapers of 0.03m each
[avbetaxtap53_63,avbetaytap53_63]=meanbetatwiss2(twissfilename,sbegtap53_63,sendtap53_63);

