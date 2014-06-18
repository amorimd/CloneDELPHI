function [elem_names,betax,betay,len]=meanbeta_triplets_HLLHC(twissfilename);

% compute all beta functions of the elements (beam screens, tapers, warm pipe, BPM) in the HL-LHC triplets (IR1 & 5).
% twissfilename is the twiss file name
% in output: names of each kind of element, beta functions and length (tables)

elem_names={'BS_101','taperBS101_BS121','BS_121','BPMCW','BPMs'};
ringlength=26658.8832;

% IP5 position (different for b1 w.r.t b2 by ~30cm)
if (strfind(twissfilename,'b1')>0)
    sIP5=13329.28923;
else
    sIP5=13329.59397;
end

% BS 101: Q1 in IR1 & 5
BSlengt=9.5;sfromIP1=22.4; % length and position from IP1 of one BS of diameter 101mm (from E. Todesco, seen in 
% G. Arduini slides, 9th WP2 task leaders' meeting, 06/06/2013)
sbeg101=[sfromIP1 sIP5-sfromIP1-BSlengt sIP5+sfromIP1 ringlength-sfromIP1-BSlengt];
send101=sbeg101+BSlengt;
len101=sum(send101-sbeg101)
[avbetax101,avbetay101]=meanbetatwiss2(twissfilename,sbeg101,send101);

% tapers from BS 101 to BS 121
IClengt=0.1; % taper length for ~6 deg angle (theta=0.1)
sbegtap101_121=[send101(1) sbeg101(2)-IClengt send101(3) sbeg101(4)-IClengt];
sendtap101_121=sbegtap101_121+IClengt;
lentap101_121=sum(sendtap101_121-sbegtap101_121)
[avbetaxtap101_121,avbetaytap101_121]=meanbetatwiss2(twissfilename,sbegtap101_121,sendtap101_121);

% BS 121: Q2 & Q3 in IR1 & 5
BSlengt=61.7-9.5-IClengt;
sbeg121=[sendtap101_121(1) sbegtap101_121(2)-BSlengt sendtap101_121(3) sbegtap101_121(4)-BSlengt];
send121=sbeg121+BSlengt;
len121=sum(send121-sbeg121)
[avbetax121,avbetay121]=meanbetatwiss2(twissfilename,sbeg121,send121);

% BPMS: take positions from R. Jones, HL-LHC PLC Jan 2013 (also in G. Arduini slides, 
% 9th WP2 task leaders' meeting, 06/06/2013)

% first BPM (CW kind - combined stripline and buttons)
% positions from IPs:
sBPMCW=[21.474];
sbegBPMCW=sort([sBPMCW ringlength-sBPMCW sIP5-sBPMCW sIP5+sBPMCW]) % 4 BPMs
sendBPMCW=sbegBPMCW+0.3; % approximate length of 0.3m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPMCW=sum(sendBPMCW-sbegBPMCW);
[avbetaxBPMCW,avbetayBPMCW]=meanbetatwiss2(twissfilename,sbegBPMCW,sendBPMCW);

% other BPMs (stripline)
% positions from IPs:
sBPMs=[32 42.667 53.263 62.737 70.456 81.193];
sbegBPMs=sort([sBPMs ringlength-sBPMs sIP5-sBPMs sIP5+sBPMs]) % 24 BPMs
sendBPMs=sbegBPMs+0.25; % approximate length of 0.25m -> see R. Jones slides Jan. 2013 HL-LHC meeting
lenBPMs=sum(sendBPMs-sbegBPMs);
[avbetaxBPMs,avbetayBPMs]=meanbetatwiss2(twissfilename,sbegBPMs,sendBPMs);


% compiling all the elements except BPMs, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg101 sbegtap101_121 sbeg121]);
sendtot=[sendtot send101 sendtap101_121 send121];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping b1');
end
sbegtot(1),sendtot(end),ringlength

betax=[avbetax101,avbetaxtap101_121,avbetax121,avbetaxBPMCW,avbetaxBPMs];
betay=[avbetay101,avbetaytap101_121,avbetay121,avbetayBPMCW,avbetayBPMs];

len=[len101,lentap101_121,len121,lenBPMCW,lenBPMs];


% write everything in a file
fid=fopen([twissfilename,'_triplets_beta_elements.dat'],'wt')
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(elem_names)
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',elem_names{i},betax(i),betay(i),len(i));
end
fclose(fid);
