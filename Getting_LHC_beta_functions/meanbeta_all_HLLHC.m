function [names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_HLLHC(twissb1,twissb2);

% compute all beta functions (HL-LHC) of the beam screens (8 different types), MBWs,
% MQWs, warm pipe, plus some tapers and BPMS in high-beta region, for both beams.
% collimators total length is taken out from the warm pipe length.
% in output: names of each kind of element, beta functions and length (tables)

names={'BS_50A','BS_50L','BS_53V','BS_53H','BS_63V','BS_63H','BS_69','BS_74','MBW','MQW_V','MQW_H','Warm_Pipe','BS_101','taperBS101_BS121','BS_121','BPMCW','BPMs','RF','CMS','ATLAS','ALICE','LHCb','BBC_H','BBC_V','Crab'};
ringlength=26658.8832;

% beam 1
[avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twissb1,1);
[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L_HLLHC(twissb1,1);
[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V_HLLHC(twissb1,1);
[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H_HLLHC(twissb1,1);
[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V_HLLHC(twissb1,1);
[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H_HLLHC(twissb1,1);
[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69_HLLHC(twissb1,1);
[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74_HLLHC(twissb1,1);
[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twissb1,1);
[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twissb1,1);
[avbetax101,avbetay101,len101,sbeg101,send101]=meanbeta_BS101(twissb1,1);
[avbetaxtap101_121,avbetaytap101_121,lentap101_121,sbegtap101_121,sendtap101_121]=meanbeta_taper101_121(twissb1,1);
[avbetax121,avbetay121,len121,sbeg121,send121]=meanbeta_BS121(twissb1,1);
[avbetaxBPMCW,avbetayBPMCW,lenBPMCW,sbegBPMCW,sendBPMCW]=meanbeta_BPMCW_HLLHC(twissb1,1);
[avbetaxBPMs,avbetayBPMs,lenBPMs,sbegBPMs,sendBPMs]=meanbeta_BPMstrip_HLLHC(twissb1,1);
[avbetaxRF,avbetayRF,lenRF,sbegRF,sendRF]=meanbeta_RF(twissb1,1);
[avbetaxCMS,avbetayCMS,lenCMS,sbegCMS,sendCMS]=meanbeta_CMS(twissb1,1);
[avbetaxATLAS,avbetayATLAS,lenATLAS,sbegATLAS,sendATLAS]=meanbeta_ATLAS(twissb1,1);
[avbetaxALICE,avbetayALICE,lenALICE,sbegALICE,sendALICE]=meanbeta_ALICE(twissb1,1);
[avbetaxLHCb,avbetayLHCb,lenLHCb,sbegLHCb,sendLHCb]=meanbeta_LHCb(twissb1,1);
[avbetaxBBC_H,avbetayBBC_H,lenBBC_H,sbegBBC_H,sendBBC_H]=meanbeta_BBC_H(twissb1,1);
[avbetaxBBC_V,avbetayBBC_V,lenBBC_V,sbegBBC_V,sendBBC_V]=meanbeta_BBC_V(twissb1,1);
[avbetaxCrab,avbetayCrab,lenCrab,sbegCrab,sendCrab]=meanbeta_Crab(twissb1,1);


% compiling all the elements except BPMs, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW sbeg101 sbegtap101_121 sbeg121]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW send101 sendtap101_121 send121];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot)
sum(sendtot-sbegtot) % total length
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping b1');
end
sbegtot(1),sendtot(end),ringlength

% rest of the machine (where warm pipe) at coll3.5tev with 2m squeeze
% Note: we take out the collimators from this (but not when computing the average beta...)

sbegrest=[0 sendtot];
sendrest=[sbegtot ringlength];
lenrest=sum(sendrest-sbegrest) % total length with collimators
% tests
if (sum(sendrest-sbegrest)+sum(sendtot-sbegtot)-ringlength~=0)
  disp(['Length problem ',numstr(sum(sendrest-sbegrest)+sum(sendtot-sbegtot)-ringlength)]);
end
if ( min(sbegrest(2:end)-sendrest(1:end-1))<0 ) || ( (sendrest(end)-sbegrest(1)-ringlength)>0 )
    disp('Aie... Overlapping b1 (with rest)');
end
% RESULT: average beta functions for the warm pipe (rest of the machine)
[avbetaxrest,avbetayrest]=meanbetatwiss2(twissb1,sbegrest,sendrest);
% avbetaxrest=145.3546, avbetayrest=146.2521 (injection b1)
% avbetaxrest=217.0115, avbetayrest=216.0646 (collision sq 2m b1)
lenrest=lenrest-55.5; % 55.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss(twissb1);

avbetax2=(avbetax50A*sum(send50A-sbeg50A)+avbetax50L*sum(send50L-sbeg50L)+ ...
    avbetax53V*sum(send53V-sbeg53V)+avbetax53H*sum(send53H-sbeg53H)+ ...
    avbetax63V*sum(send63V-sbeg63V)+avbetax63H*sum(send63H-sbeg63H)+ ...
    avbetax69*sum(send69-sbeg69)+avbetax74*sum(send74-sbeg74)+ ...
    avbetaxMBW*sum(sendMBW-sbegMBW)+avbetaxMQW*sum(sendMQW-sbegMQW)+ ...
    avbetax101*(send101-sbeg101)+avbetax121*(send121-sbeg121)+ ...
    avbetaxrest*sum(sendrest-sbegrest))/ringlength;
avbetay2=(avbetay50A*sum(send50A-sbeg50A)+avbetay50L*sum(send50L-sbeg50L)+ ...
    avbetay53V*sum(send53V-sbeg53V)+avbetay53H*sum(send53H-sbeg53H)+ ...
    avbetay63V*sum(send63V-sbeg63V)+avbetay63H*sum(send63H-sbeg63H)+ ...
    avbetay69*sum(send69-sbeg69)+avbetay74*sum(send74-sbeg74)+ ...
    avbetayMBW*sum(sendMBW-sbegMBW)+avbetayMQW*sum(sendMQW-sbegMQW)+ ...
    avbetay101*(send101-sbeg101)+avbetay121*(send121-sbeg121)+ ...
    avbetayrest*sum(sendrest-sbegrest))/ringlength;

if ((avbetax2-avbetax)/avbetax>1e-9)
  disp(['difference between averages (x, b1): ',num2str((avbetax2-avbetax)/avbetax)]);
end
if ((avbetay2-avbetay)/avbetay>1e-9)
  disp(['difference between averages (y, b1): ',num2str((avbetay2-avbetay)/avbetay)]);
end

betaxb1=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetax101,avbetaxtap101_121,avbetax121,avbetaxBPMCW,avbetaxBPMs,avbetaxRF,avbetaxCMS,avbetaxATLAS,avbetaxALICE,avbetaxLHCb,avbetaxBBC_H,avbetaxBBC_V,avbetaxCrab];
betayb1=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetay101,avbetaytap101_121,avbetay121,avbetayBPMCW,avbetayBPMs,avbetayRF,avbetayCMS,avbetayATLAS,avbetayALICE,avbetayLHCb,avbetayBBC_H,avbetayBBC_V,avbetayCrab];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,len101,lentap101_121,len121,lenBPMCW,lenBPMs,lenRF,lenCMS,lenATLAS,lenALICE,lenLHCb,lenBBC_H,lenBBC_V,lenCrab];


% write everything in a file
fid=fopen([twissb1,'_beta_elements.dat'],'wt')
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(names)
  %fprintf(fid,'%s \t %lf \t %lf \t %lf\n',names{i},betaxb1(i),betayb1(i),len(i));
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',names{i},betaxb1(i),betayb1(i),len(i));
end
fclose(fid);



% beam 2
[avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twissb2,2);
[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L_HLLHC(twissb2,2);
[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V_HLLHC(twissb2,2);
[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H_HLLHC(twissb2,2);
[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V_HLLHC(twissb2,2);
[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H_HLLHC(twissb2,2);
[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69_HLLHC(twissb2,2);
[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74_HLLHC(twissb2,2);
[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twissb2,2);
[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twissb2,2);
[avbetax101,avbetay101,len101,sbeg101,send101]=meanbeta_BS101(twissb2,2);
[avbetaxtap101_121,avbetaytap101_121,lentap101_121,sbegtap101_121,sendtap101_121]=meanbeta_taper101_121(twissb2,2);
[avbetax121,avbetay121,len121,sbeg121,send121]=meanbeta_BS121(twissb2,2);
[avbetaxBPMCW,avbetayBPMCW,lenBPMCW,sbegBPMCW,sendBPMCW]=meanbeta_BPMCW_HLLHC(twissb2,2);
[avbetaxBPMs,avbetayBPMs,lenBPMs,sbegBPMs,sendBPMs]=meanbeta_BPMstrip_HLLHC(twissb2,2);
[avbetaxRF,avbetayRF,lenRF,sbegRF,sendRF]=meanbeta_RF(twissb2,2);
[avbetaxCMS,avbetayCMS,lenCMS,sbegCMS,sendCMS]=meanbeta_CMS(twissb2,2);
[avbetaxATLAS,avbetayATLAS,lenATLAS,sbegATLAS,sendATLAS]=meanbeta_ATLAS(twissb2,2);
[avbetaxALICE,avbetayALICE,lenALICE,sbegALICE,sendALICE]=meanbeta_ALICE(twissb2,2);
[avbetaxLHCb,avbetayLHCb,lenLHCb,sbegLHCb,sendLHCb]=meanbeta_LHCb(twissb2,2);
[avbetaxBBC_H,avbetayBBC_H,lenBBC_H,sbegBBC_H,sendBBC_H]=meanbeta_BBC_H(twissb2,2);
[avbetaxBBC_V,avbetayBBC_V,lenBBC_V,sbegBBC_V,sendBBC_V]=meanbeta_BBC_V(twissb2,2);
[avbetaxCrab,avbetayCrab,lenCrab,sbegCrab,sendCrab]=meanbeta_Crab(twissb2,2);

% compiling all the elements, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW  sbeg101 sbegtap101_121 sbeg121]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW send101 sendtap101_121 send121];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) %1850=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8+21+48 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % % total length : 23016.6754
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping b2');
end
sbegtot(1),sendtot(end),ringlength

% rest of the machine (where warm pipe) at coll3.5tev with 2m squeeze
% Note: we take out the collimators from this (but not when computing the average beta...)

sbegrest=[0 sendtot];
sendrest=[sbegtot ringlength];
lenrest=sum(sendrest-sbegrest) % total length with collimators: 3642.2078
% tests
if (sum(sendrest-sbegrest)+sum(sendtot-sbegtot)-ringlength~=0)
  disp(['Length problem ',numstr(sum(sendrest-sbegrest)+sum(sendtot-sbegtot)-ringlength)]);
end
if ( min(sbegrest(2:end)-sendrest(1:end-1))<0 ) || ( (sendrest(end)-sbegrest(1)-ringlength)>0 )
    disp('Aie... Overlapping b2 (with rest)');
end
% RESULT: average beta functions for the warm pipe (rest of the machine)
[avbetaxrest,avbetayrest]=meanbetatwiss2(twissb2,sbegrest,sendrest);
lenrest=lenrest-55.5; % 55.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss(twissb2);

avbetax2=(avbetax50A*sum(send50A-sbeg50A)+avbetax50L*sum(send50L-sbeg50L)+ ...
    avbetax53V*sum(send53V-sbeg53V)+avbetax53H*sum(send53H-sbeg53H)+ ...
    avbetax63V*sum(send63V-sbeg63V)+avbetax63H*sum(send63H-sbeg63H)+ ...
    avbetax69*sum(send69-sbeg69)+avbetax74*sum(send74-sbeg74)+ ...
    avbetaxMBW*sum(sendMBW-sbegMBW)+avbetaxMQW*sum(sendMQW-sbegMQW)+ ...
    avbetax101*(send101-sbeg101)+avbetax121*(send121-sbeg121)+ ...
    avbetaxrest*sum(sendrest-sbegrest))/ringlength;
avbetay2=(avbetay50A*sum(send50A-sbeg50A)+avbetay50L*sum(send50L-sbeg50L)+ ...
    avbetay53V*sum(send53V-sbeg53V)+avbetay53H*sum(send53H-sbeg53H)+ ...
    avbetay63V*sum(send63V-sbeg63V)+avbetay63H*sum(send63H-sbeg63H)+ ...
    avbetay69*sum(send69-sbeg69)+avbetay74*sum(send74-sbeg74)+ ...
    avbetayMBW*sum(sendMBW-sbegMBW)+avbetayMQW*sum(sendMQW-sbegMQW)+ ...
    avbetay101*(send101-sbeg101)+avbetay121*(send121-sbeg121)+ ...
    avbetayrest*sum(sendrest-sbegrest))/ringlength;

if ((avbetax2-avbetax)/avbetax>1e-9)
  disp(['difference between averages (x, b2): ',num2str((avbetax2-avbetax)/avbetax)]);
end
if ((avbetay2-avbetay)/avbetay>1e-9)
  disp(['difference between averages (y, b2): ',num2str((avbetay2-avbetay)/avbetay)]);
end

betaxb2=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetax101,avbetaxtap101_121,avbetax121,avbetaxBPMCW,avbetaxBPMs,avbetaxRF,avbetaxCMS,avbetaxATLAS,avbetaxALICE,avbetaxLHCb,avbetaxBBC_H,avbetaxBBC_V,avbetaxCrab];
betayb2=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetay101,avbetaytap101_121,avbetay121,avbetayBPMCW,avbetayBPMs,avbetayRF,avbetayCMS,avbetayATLAS,avbetayALICE,avbetayLHCb,avbetayBBC_H,avbetayBBC_V,avbetayCrab];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,len101,lentap101_121,len121,lenBPMCW,lenBPMs,lenRF,lenCMS,lenATLAS,lenALICE,lenLHCb,lenBBC_H,lenBBC_V,lenCrab];


% write everything in a file
fid=fopen([twissb2,'_beta_elements.dat'],'wt');
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(names)
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',names{i},betaxb2(i),betayb2(i),len(i));
end
fclose(fid);
