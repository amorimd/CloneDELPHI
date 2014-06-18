function [names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all(twissb1,twissb2);

% compute all beta functions of the beam screens (8 different types), MBWs,
% MQWs, warm pipe, plus some tapers and BPMS in high-beta region, for both beams.
% collimators total length is taken out from the warm pipe length.
% in output: names of each kind of element, beta functions and length (tables)

names={'BS_50A','BS_50L','BS_53V','BS_53H','BS_63V','BS_63H','BS_69','BS_74','MBW','MQW_V','MQW_H','Warm_Pipe','taperBS53_BS63','taperBS63_BS74','taperBS74_warm','BPM61CW','BPM61strip','BPM80strip','RF','CMS','ALICE','LHCb'};
ringlength=26658.8832;

% beam 1
[avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twissb1,1);
[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L(twissb1,1);
[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V(twissb1,1);
[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H(twissb1,1);
[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V(twissb1,1);
[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H(twissb1,1);
[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69(twissb1,1);
[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74(twissb1,1);
[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twissb1,1);
[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twissb1,1);
[avbetaxtap53_63,avbetaytap53_63,lentap53_63,sbegtap53_63,sendtap53_63]=meanbeta_taper53_63(twissb1,1);
[avbetaxtap63_74,avbetaytap63_74,lentap63_74,sbegtap63_74,sendtap63_74]=meanbeta_taper63_74(twissb1,1);
[avbetaxtap74_WP,avbetaytap74_WP,lentap74_WP,sbegtap74_WP,sendtap74_WP]=meanbeta_taper74_WP(twissb1,1);
[avbetaxBPM61CW,avbetayBPM61CW,lenBPM61CW,sbegBPM61CW,sendBPM61CW]=meanbeta_BPM61CW(twissb1,1);
[avbetaxBPM61,avbetayBPM61,lenBPM61,sbegBPM61,sendBPM61]=meanbeta_BPM61(twissb1,1);
[avbetaxBPM80,avbetayBPM80,lenBPM80,sbegBPM80,sendBPM80]=meanbeta_BPM80(twissb1,1);
[avbetaxRF,avbetayRF,lenRF,sbegRF,sendRF]=meanbeta_RF(twissb1,1);
[avbetaxCMS,avbetayCMS,lenCMS,sbegCMS,sendCMS]=meanbeta_CMS(twissb1,1);
[avbetaxALICE,avbetayALICE,lenALICE,sbegALICE,sendALICE]=meanbeta_ALICE(twissb1,1);
[avbetaxLHCb,avbetayLHCb,lenLHCb,sbegLHCb,sendLHCb]=meanbeta_LHCb(twissb1,1);

% compiling all the elements except BPMs, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW sbegtap53_63 sbegtap63_74 sbegtap74_WP]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW sendtap53_63 sendtap63_74 sendtap74_WP];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) %1862=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8+21+48+4*3 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % total length : 23017
if ( min(sbegtot(2:end)-sendtot(1:end-1))<0 ) || ( (sendtot(end)-sbegtot(1)-ringlength)>0 )
    disp('Aie... Overlapping b1');
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
    disp('Aie... Overlapping b1 (with rest)');
end
% RESULT: average beta functions for the warm pipe (rest of the machine)
[avbetaxrest,avbetayrest]=meanbetatwiss2(twissb1,sbegrest,sendrest);
% avbetaxrest=145.3546, avbetayrest=146.2521 (injection b1)
% avbetaxrest=217.0115, avbetayrest=216.0646 (collision sq 2m b1)
lenrest=lenrest-49.5; % 49.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss(twissb1);

avbetax2=(avbetax50A*sum(send50A-sbeg50A)+avbetax50L*sum(send50L-sbeg50L)+ ...
    avbetax53V*sum(send53V-sbeg53V)+avbetax53H*sum(send53H-sbeg53H)+ ...
    avbetax63V*sum(send63V-sbeg63V)+avbetax63H*sum(send63H-sbeg63H)+ ...
    avbetax69*sum(send69-sbeg69)+avbetax74*sum(send74-sbeg74)+ ...
    avbetaxMBW*sum(sendMBW-sbegMBW)+avbetaxMQW*sum(sendMQW-sbegMQW)+ ...
    avbetaxrest*sum(sendrest-sbegrest))/ringlength;
avbetay2=(avbetay50A*sum(send50A-sbeg50A)+avbetay50L*sum(send50L-sbeg50L)+ ...
    avbetay53V*sum(send53V-sbeg53V)+avbetay53H*sum(send53H-sbeg53H)+ ...
    avbetay63V*sum(send63V-sbeg63V)+avbetay63H*sum(send63H-sbeg63H)+ ...
    avbetay69*sum(send69-sbeg69)+avbetay74*sum(send74-sbeg74)+ ...
    avbetayMBW*sum(sendMBW-sbegMBW)+avbetayMQW*sum(sendMQW-sbegMQW)+ ...
    avbetayrest*sum(sendrest-sbegrest))/ringlength;

if ((avbetax2-avbetax)/avbetax>1e-9)
  disp(['difference between averages (x, b1): ',num2str((avbetax2-avbetax)/avbetax)]);
end
if ((avbetay2-avbetay)/avbetay>1e-9)
  disp(['difference between averages (y, b1): ',num2str((avbetay2-avbetay)/avbetay)]);
end

betaxb1=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetaxtap53_63,avbetaxtap63_74,avbetaxtap74_WP,avbetaxBPM61CW,avbetaxBPM61,avbetaxBPM80,avbetaxRF,avbetaxCMS,avbetaxALICE,avbetaxLHCb];
betayb1=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetaytap53_63,avbetaytap63_74,avbetaytap74_WP,avbetayBPM61CW,avbetayBPM61,avbetayBPM80,avbetayRF,avbetayCMS,avbetayALICE,avbetayLHCb];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,lentap53_63,lentap63_74,lentap74_WP,lenBPM61CW,lenBPM61,lenBPM80,lenRF,lenCMS,lenALICE,lenLHCb];


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
[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L(twissb2,2);
[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V(twissb2,2);
[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H(twissb2,2);
[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V(twissb2,2);
[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H(twissb2,2);
[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69(twissb2,2);
[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74(twissb2,2);
[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twissb2,2);
[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twissb2,2);
[avbetaxtap53_63,avbetaytap53_63,lentap53_63,sbegtap53_63,sendtap53_63]=meanbeta_taper53_63(twissb2,2);
[avbetaxtap63_74,avbetaytap63_74,lentap63_74,sbegtap63_74,sendtap63_74]=meanbeta_taper63_74(twissb2,2);
[avbetaxtap74_WP,avbetaytap74_WP,lentap74_WP,sbegtap74_WP,sendtap74_WP]=meanbeta_taper74_WP(twissb2,2);
[avbetaxBPM61CW,avbetayBPM61CW,lenBPM61CW,sbegBPM61CW,sendBPM61CW]=meanbeta_BPM61CW(twissb2,2);
[avbetaxBPM61,avbetayBPM61,lenBPM61,sbegBPM61,sendBPM61]=meanbeta_BPM61(twissb2,2);
[avbetaxBPM80,avbetayBPM80,lenBPM80,sbegBPM80,sendBPM80]=meanbeta_BPM80(twissb2,2);
[avbetaxRF,avbetayRF,lenRF,sbegRF,sendRF]=meanbeta_RF(twissb2,2);
[avbetaxCMS,avbetayCMS,lenCMS,sbegCMS,sendCMS]=meanbeta_CMS(twissb2,2);
[avbetaxALICE,avbetayALICE,lenALICE,sbegALICE,sendALICE]=meanbeta_ALICE(twissb2,2);
[avbetaxLHCb,avbetayLHCb,lenLHCb,sbegLHCb,sendLHCb]=meanbeta_LHCb(twissb2,2);

% compiling all the elements, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW sbegtap53_63 sbegtap63_74 sbegtap74_WP]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW sendtap53_63 sendtap63_74 sendtap74_WP];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) %1862=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8+21+48+4*3 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % % total length : 23017
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
lenrest=lenrest-49.5; % 49.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss(twissb2);

avbetax2=(avbetax50A*sum(send50A-sbeg50A)+avbetax50L*sum(send50L-sbeg50L)+ ...
    avbetax53V*sum(send53V-sbeg53V)+avbetax53H*sum(send53H-sbeg53H)+ ...
    avbetax63V*sum(send63V-sbeg63V)+avbetax63H*sum(send63H-sbeg63H)+ ...
    avbetax69*sum(send69-sbeg69)+avbetax74*sum(send74-sbeg74)+ ...
    avbetaxMBW*sum(sendMBW-sbegMBW)+avbetaxMQW*sum(sendMQW-sbegMQW)+ ...
    avbetaxrest*sum(sendrest-sbegrest))/ringlength;
avbetay2=(avbetay50A*sum(send50A-sbeg50A)+avbetay50L*sum(send50L-sbeg50L)+ ...
    avbetay53V*sum(send53V-sbeg53V)+avbetay53H*sum(send53H-sbeg53H)+ ...
    avbetay63V*sum(send63V-sbeg63V)+avbetay63H*sum(send63H-sbeg63H)+ ...
    avbetay69*sum(send69-sbeg69)+avbetay74*sum(send74-sbeg74)+ ...
    avbetayMBW*sum(sendMBW-sbegMBW)+avbetayMQW*sum(sendMQW-sbegMQW)+ ...
    avbetayrest*sum(sendrest-sbegrest))/ringlength;

if ((avbetax2-avbetax)/avbetax>1e-9)
  disp(['difference between averages (x, b2): ',num2str((avbetax2-avbetax)/avbetax)]);
end
if ((avbetay2-avbetay)/avbetay>1e-9)
  disp(['difference between averages (y, b2): ',num2str((avbetay2-avbetay)/avbetay)]);
end

betaxb2=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetaxtap53_63,avbetaxtap63_74,avbetaxtap74_WP,avbetaxBPM61CW,avbetaxBPM61,avbetaxBPM80,avbetaxRF,avbetaxCMS,avbetaxALICE,avbetaxLHCb];
betayb2=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetaytap53_63,avbetaytap63_74,avbetaytap74_WP,avbetayBPM61CW,avbetayBPM61,avbetayBPM80,avbetayRF,avbetayCMS,avbetayALICE,avbetayLHCb];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,lentap53_63,lentap63_74,lentap74_WP,lenBPM61CW,lenBPM61,lenBPM80,lenRF,lenCMS,lenALICE,lenLHCb];


% write everything in a file
fid=fopen([twissb2,'_beta_elements.dat'],'wt');
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(names)
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',names{i},betaxb2(i),betayb2(i),len(i));
end
fclose(fid);
