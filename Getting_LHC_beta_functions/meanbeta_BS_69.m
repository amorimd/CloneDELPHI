function [avbetax69,avbetay69,len69,sbeg69,send69]=meanbeta_BS_69(twissfilename,beam);


% function to compute the average beta functions of all beam screens of
% type 69.
% in input: name of the twiss file to use for this. beam=1 or 2
% in output: average beta functions in x and y, total length of the beam
% screens, and beginning and end (in terms of longitudinal position s) of each of them.



%%%%
% beam screen type 69, with vertical half gap b=56.15/2
%%%%

% D2 in IR1, IR2, IR5 and IR8 (in Q4) and D4 in IR4 (in Q5)
[sbegD2_4_IR1_2_4_5_8,sendD2_4_IR1_2_4_5_8,names,len]=extract_beamscreens(twissfilename,'MBR((C\.4[RL][1258])|(B\.5[RL][4]))',10.7325);
%size(sbegD2_4_IR1_2_4_5_8),names,sum(len),sum(sendD2_4_IR1_2_4_5_8-sbegD2_4_IR1_2_4_5_8)
% 10 BS. total BS length=1.073249999999879e+02

% D3 in IR4 (in Q5)
[sbegD3_IR4,sendD3_IR4,namesD3_IR4,len]=extract_beamscreens(twissfilename,'MBRS\.5[RL]4',10.9895);
%size(sbegD3_IR4),namesD3_IR4,sum(len),sum(sendD3_IR4-sbegD3_IR4)

if (beam==1)
  % we offset the right one a bit toward the right to avoid overlapping 
  % with the beam screen of the undulator
  sbegD3_IR4(2)=sbegD3_IR4(2)+0.3;sendD3_IR4(2)=sendD3_IR4(2)+0.3;
  % 2 BS. total BS length=21.978999999999360
  % undulator left in IR4
  % it is next to the MBRS of D3 (going toward IP4)
  BSlengt=1.3405;IClengt=0.165; % interconnect length is a rough estimate
  %namesD3_IR4
  sbegUnd_IR4=sendD3_IR4(1)+IClengt;
  sendUnd_IR4=sbegUnd_IR4+BSlengt;
elseif (beam==2)
  % we offset the left one a bit toward the left to avoid overlapping 
  % with the beam screen of the undulator
  sbegD3_IR4(1)=sbegD3_IR4(1)-0.3;sendD3_IR4(1)=sendD3_IR4(1)-0.3;
  % 2 BS. total BS length=21.978999999999360
  % undulator left in IR4
  % it is next to the MBRS of D3 (going toward IP4)
  BSlengt=1.3405;IClengt=0.165; % interconnect length is a rough estimate
  %namesD3_IR4
  sendUnd_IR4=sbegD3_IR4(2)-IClengt;
  sbegUnd_IR4=sendUnd_IR4-BSlengt;  
end
%size(sbegUnd_IR4),sum(sendUnd_IR4-sbegUnd_IR4)
% 1 BS. total BS length=1.3405

% compiling all the type 69 vertical beam screens
[sbeg69,ind]=sort([sbegD2_4_IR1_2_4_5_8 sbegD3_IR4 sbegUnd_IR4]);
send69=[sendD2_4_IR1_2_4_5_8 sendD3_IR4 sendUnd_IR4];
send69=send69(ind);
length(sbeg69),length(send69) % 13
len69=sum(send69-sbeg69) % total length: 130.6445 
ringlength=26658.8832;
% overlapping test
if ( min(sbeg69(2:end)-send69(1:end-1))<0 ) || ( (send69(end)-sbeg69(1)-ringlength)>0 )
    disp('Aie... Overlapping beam screens');
end


% RESULT: average beta functions for the type 69 vertical BS
[avbetax69,avbetay69]=meanbetatwiss2(twissfilename,sbeg69,send69);
% avbetax69=151.9529, avbetay69=131.691 (injection b1)
% avbetax69=295.2808, avbetay69=275.0131 (collision sq 2m b1)

