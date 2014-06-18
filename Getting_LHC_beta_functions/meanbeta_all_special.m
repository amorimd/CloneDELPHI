function [names,betaxb1,betayb1,betaxb2,betayb2,len]=meanbeta_all_special(twissb1,twissb2,twiss_ref_b1,twiss_ref_b2);

% compute all beta functions of the beam screens (8 different types), MBWs, MQWs
% and warm pipe, and for both beams.
% collimators total length is taken out from the warm pipe length.
% in output: names of each kind of element, beta functions and length (tables)

% special case for HL-LHC very preliminary model: use another twiss (from
% 2012 machine typically) to get beginning and end of each device
% (twiss_ref_b[12]) and twissb[12] to get beta functions there.
% This is to be used with great CAUTION.

names={'BS_50A','BS_50L','BS_53V','BS_53H','BS_63V','BS_63H','BS_69','BS_74','MBW','MQW_V','MQW_H','Warm_Pipe','Total'};
ringlength=26658.8832;

% beam 1
[avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twiss_ref_b1,1);
[avbetax50A,avbetay50A]=meanbetatwiss2_special(twissb1,sbeg50A,send50A);

[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L(twiss_ref_b1,1);
[avbetax50L,avbetay50L]=meanbetatwiss2_special(twissb1,sbeg50L,send50L);

[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V(twiss_ref_b1,1);
[avbetax53V,avbetay53V]=meanbetatwiss2_special(twissb1,sbeg53V,send53V);

[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H(twiss_ref_b1,1);
[avbetax53H,avbetay53H]=meanbetatwiss2_special(twissb1,sbeg53H,send53H);

[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V(twiss_ref_b1,1);
[avbetax63V,avbetay63V]=meanbetatwiss2_special(twissb1,sbeg63V,send63V);

[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H(twiss_ref_b1,1);
[avbetax63H,avbetay63H]=meanbetatwiss2_special(twissb1,sbeg63H,send63H);

[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69(twiss_ref_b1,1);
[avbetax69,avbetay69]=meanbetatwiss2_special(twissb1,sbeg69,send69);

[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74(twiss_ref_b1,1);
[avbetax74,avbetay74]=meanbetatwiss2_special(twissb1,sbeg74,send74);

[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twiss_ref_b1,1);
[avbetaxMBW,avbetayMBW]=meanbetatwiss2_special(twissb1,sbegMBW,sendMBW);

[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twiss_ref_b1,1);
[avbetaxMQW,avbetayMQW]=meanbetatwiss2_special(twissb1,sbegMQW,sendMQW);

% compiling all the elements, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW];
sendtot=sendtot(ind);
length(sbegtot),length(sendtot) %1850=1232+362+52+16+8+8+8+1+8+8+2+4+4+1+8+8+8+4+4+4+6+4+4+1+10+2+2+2+4+4-16+8+21+48 elements
%(16 drift BS are together with quad. ones and 8 Q1 extensions are artificially separated)
sum(sendtot-sbegtot) % % total length : 23016.6754
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
[avbetaxrest,avbetayrest]=meanbetatwiss2_special(twissb1,sbegrest,sendrest);
% avbetaxrest=145.3546, avbetayrest=146.2521 (injection b1)
% avbetaxrest=217.0115, avbetayrest=216.0646 (collision sq 2m b1)
lenrest=lenrest-49.5; % 49.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss_special(twissb1);

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

betaxb1=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetax];
betayb1=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetay];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,ringlength];


% write everything in a file
fid=fopen([twissb1,'_beta_elements.dat'],'wt');
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(names)
  %fprintf(fid,'%s \t %lf \t %lf \t %lf\n',names{i},betaxb1(i),betayb1(i),len(i));
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',names{i},betaxb1(i),betayb1(i),len(i));
end
fclose(fid);



% beam 2
[avbetax50A,avbetay50A,len50A,sbeg50A,send50A]=meanbeta_BS_50A(twiss_ref_b2,2);
[avbetax50A,avbetay50A]=meanbetatwiss2_special(twissb2,sbeg50A,send50A);

[avbetax50L,avbetay50L,len50L,sbeg50L,send50L]=meanbeta_BS_50L(twiss_ref_b2,2);
[avbetax50L,avbetay50L]=meanbetatwiss2_special(twissb2,sbeg50L,send50L);

[avbetax53V,avbetay53V,len53V,sbeg53V,send53V]=meanbeta_BS_53V(twiss_ref_b2,2);
[avbetax53V,avbetay53V]=meanbetatwiss2_special(twissb2,sbeg53V,send53V);

[avbetax53H,avbetay53H,len53H,sbeg53H,send53H]=meanbeta_BS_53H(twiss_ref_b2,2);
[avbetax53H,avbetay53H]=meanbetatwiss2_special(twissb2,sbeg53H,send53H);

[avbetax63V,avbetay63V,len63V,sbeg63V,send63V]=meanbeta_BS_63V(twiss_ref_b2,2);
[avbetax63V,avbetay63V]=meanbetatwiss2_special(twissb2,sbeg63V,send63V);

[avbetax63H,avbetay63H,len63H,sbeg63H,send63H]=meanbeta_BS_63H(twiss_ref_b2,2);
[avbetax63H,avbetay63H]=meanbetatwiss2_special(twissb2,sbeg63H,send63H);

[avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69(twiss_ref_b2,2);
[avbetax69,avbetay69]=meanbetatwiss2_special(twissb2,sbeg69,send69);

[avbetax74,avbetay74,len74,sbeg74,send74]=meanbeta_BS_74(twiss_ref_b2,2);
[avbetax74,avbetay74]=meanbetatwiss2_special(twissb2,sbeg74,send74);

[avbetaxMBW,avbetayMBW,lenMBW,sbegMBW,sendMBW]=meanbeta_MBW(twiss_ref_b2,2);
[avbetaxMBW,avbetayMBW]=meanbetatwiss2_special(twissb2,sbegMBW,sendMBW);

[avbetaxMQW,avbetayMQW,lenMQW,sbegMQW,sendMQW]=meanbeta_MQW(twiss_ref_b2,2);
[avbetaxMQW,avbetayMQW]=meanbetatwiss2_special(twissb2,sbegMQW,sendMQW);


% compiling all the elements, and operlapping test
sbegtot=[];sendtot=[];
[sbegtot,ind]=sort([sbegtot sbeg50A sbeg50L sbeg53H sbeg53V sbeg63H sbeg63V sbeg69 sbeg74 sbegMBW sbegMQW]);
sendtot=[sendtot send50A send50L send53H send53V send63H send63V send69 send74 sendMBW sendMQW];
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
[avbetaxrest,avbetayrest]=meanbetatwiss2_special(twissb2,sbegrest,sendrest);
lenrest=lenrest-49.5; % 49.5m is the total length of collimators

% small check on the average beta functions
[s,betax,betay,avbetax,avbetay]=meanbetatwiss_special(twissb2);

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

betaxb2=[avbetax50A,avbetax50L,avbetax53V,avbetax53H,avbetax63V,avbetax63H,avbetax69,avbetax74,avbetaxMBW,avbetaxMQW,avbetaxMQW,avbetaxrest,avbetax];
betayb2=[avbetay50A,avbetay50L,avbetay53V,avbetay53H,avbetay63V,avbetay63H,avbetay69,avbetay74,avbetayMBW,avbetayMQW,avbetayMQW,avbetayrest,avbetay];

len=[len50A,len50L,len53V,len53H,len63V,len63H,len69,len74,lenMBW,lenMQW/2,lenMQW/2,lenrest,ringlength];


% write everything in a file
fid=fopen([twissb2,'_beta_elements.dat'],'wt');
fprintf(fid,'name \t betax \t betay \t length\n');
for i=1:length(names)
  fprintf(fid,'%s\t%10.5e\t%10.5e\t%13.8e\n',names{i},betaxb2(i),betayb2(i),len(i));
end
fclose(fid);
