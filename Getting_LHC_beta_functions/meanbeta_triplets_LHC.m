function [elem_names,betax,betay,len]=meanbeta_triplets_LHC(twissfilename);

% compute all beta functions of the elements (beam screens, tapers, warm pipe, BPM) in the LHC triplets (IR1 & 5).
% twissfilename is the twiss file name
% in output: names of each kind of element, beta functions and length (tables)

elem_names={'BS_53H','BS_53V','BS_63H','BS_63V','BS_74','Warm_Pipe','taperBS53_BS63','taperBS63_BS74','taperBS74_warm','BPM61CW','BPM61strip','BPM80strip'};
ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (strfind(twissfilename,'b1')>0)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% BS 53H: Q1 in IR1
[sbeg53H,send53H,names,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL]1',8.6815-0.77-0.03);
[avbetax53H,avbetay53H]=meanbetatwiss2(twissfilename,sbeg53H,send53H);
len53H=sum(send53H-sbeg53H)

% BS 53V: Q1 in IR5
[sbeg53V,send53V,names,len]=extract_beamscreens(twissfilename,'MQXA\.1[RL]5',8.6815-0.77-0.03);
[avbetax53V,avbetay53V]=meanbetatwiss2(twissfilename,sbeg53V,send53V);
len53V=sum(send53V-sbeg53V)

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
length(sbeg63H),length(send63H) % 6
len63H=sum(send63H-sbeg63H)
[avbetax63H,avbetay63H]=meanbetatwiss2(twissfilename,sbeg63H,send63H);

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
length(sbeg63V),length(send63V) % 6
len63V=sum(send63V-sbeg63V)
[avbetax63V,avbetay63V]=meanbetatwiss2(twissfilename,sbeg63V,send63V);

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
length(sbeg74),length(send74) % 4
len74=sum(send74-sbeg74)
[avbetax74,avbetay74]=meanbetatwiss2(twissfilename,sbeg74,send74);

% warm pipe in D1
% we take it right outside the BS 74 (until the end of D1) (length computed
% from plan)
WPlengt=27.25;IClengt=0.035;
sbegWP=[send74(1)+IClengt sbeg74(2)-WPlengt-IClengt send74(3)+IClengt sbeg74(4)-WPlengt-IClengt];
sendWP=sbegWP+WPlengt;
lenWP=sum(sendWP-sbegWP)
[avbetaxWP,avbetayWP]=meanbetatwiss2(twissfilename,sbegWP,sendWP);

% tapers from BS 53 to BS 63 (note: tapers are round -> do not depend on BS
% orientation)
sbegtap53_63=[send53H(1) send63H(6) send53V(2) send63V(3)]
sendtap53_63=[sbeg63H(1) sbeg53H(2) sbeg63V(4) sbeg53V(1)]
sendtap53_63-sbegtap53_63
lentap53_63=sum(sendtap53_63-sbegtap53_63) % 4 tapers of 0.03m each
[avbetaxtap53_63,avbetaytap53_63]=meanbetatwiss2(twissfilename,sbegtap53_63,sendtap53_63);

% tapers from BS 63 to BS 74 (note: tapers are round -> do not depend on BS
% orientation)
sbegtap63_74=[send63H(3) send74(4) send63V(6) send74(2)]
sendtap63_74=[sbeg74(1) sbeg63H(4) sbeg74(3) sbeg63V(1)]
sendtap63_74-sbegtap63_74
lentap63_74=sum(sendtap63_74-sbegtap63_74) % 4 tapers of 0.077m each
[avbetaxtap63_74,avbetaytap63_74]=meanbetatwiss2(twissfilename,sbegtap63_74,sendtap63_74);

% tapers from BS 74 to warm pipe
sbegtap74_WP=[send74(1) sendWP(4) send74(3) sendWP(2)]
sendtap74_WP=[sbegWP(1) sbeg74(4) sbegWP(3) sbeg74(2)]
sendtap74_WP-sbegtap74_WP
lentap74_WP=sum(sendtap74_WP-sbegtap74_WP) % 4 tapers of 0.035m each
[avbetaxtap74_WP,avbetaytap74_WP]=meanbetatwiss2(twissfilename,sbegtap74_WP,sendtap74_WP);

% BPMS: positions from R. Jones, HL-LHC PLC Jan 2013 (also in G. Arduini slides, 
% 9th WP2 task leaders' meeting, 06/06/2013)

% BPM before Q1 (BPM of CW kind - combined button + stripline), of diameter between electrodes 61mm
% positions from IPs:
sBPM61CW=[21.474];
sbegBPM61CW=sort([sBPM61CW ringlength-sBPM61CW sIP5-sBPM61CW sIP5+sBPM61CW]) % 4 BPMs
sendBPM61CW=sbegBPM61CW+0.3; % approximate length of 0.3m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM61CW=sum(sendBPM61CW-sbegBPM61CW);
[avbetaxBPM61CW,avbetayBPM61CW]=meanbetatwiss2(twissfilename,sbegBPM61CW,sendBPM61CW);

% second BPM, stripline, between Q1 and Q2, of diameter between electrodes 60mm
% positions from IPs:
sBPM61=[32];
sbegBPM61=sort([sBPM61 ringlength-sBPM61 sIP5-sBPM61 sIP5+sBPM61]) % 4 BPMs
sendBPM61=sbegBPM61+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM61=sum(sendBPM61-sbegBPM61);
[avbetaxBPM61,avbetayBPM61]=meanbetatwiss2(twissfilename,sbegBPM61,sendBPM61);

% third BPM, stripline, after Q3, of diameter between electrodes 80mm
% positions from IPs:
sBPM80=[70.456];
sbegBPM80=sort([sBPM80 ringlength-sBPM80 sIP5-sBPM80 sIP5+sBPM80]) % 4 BPMs
sendBPM80=sbegBPM80+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPM80=sum(sendBPM80-sbegBPM80);
[avbetaxBPM80,avbetayBPM80]=meanbetatwiss2(twissfilename,sbegBPM80,sendBPM80);


% compiling all the elements except BPMs, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg53H sbeg53V sbeg63H sbeg63V sbeg74 sbegWP sbegtap53_63 sbegtap63_74 sbegtap74_WP]);
sendtot=[sendtot send53H send53V send63H send63V send74 sendWP sendtap53_63 sendtap63_74 sendtap74_WP];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping b1');
end
sbegtot(1),sendtot(end),ringlength

betax=[avbetax53H,avbetax53V,avbetax63H,avbetax63V,avbetax74,avbetaxWP,avbetaxtap53_63,avbetaxtap63_74,avbetaxtap74_WP,avbetaxBPM61CW,avbetaxBPM61,avbetaxBPM80];
betay=[avbetay53H,avbetay53V,avbetay63H,avbetay63V,avbetay74,avbetayWP,avbetaytap53_63,avbetaytap63_74,avbetaytap74_WP,avbetayBPM61CW,avbetayBPM61,avbetayBPM80];

len=[len53H,len53V,len63H,len63V,len74,lenWP,lentap53_63,lentap63_74,lentap74_WP,lenBPM61CW,lenBPM61,lenBPM80];


% write everything in a file
fid=fopen([twissfilename,'_triplets_beta_elements.dat'],'wt')
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(elem_names)
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',elem_names{i},betax(i),betay(i),len(i));
end
fclose(fid);
