function [a,qsec]=octupole(twissname);

% computes the maximum detuning coefficients from the octupoles and maximum
% Q'' due to the octupoles, at 7 TeV/c (these are inversely proportional to the energy).

% a is a matrix 3*2: [axF  axD;
%                     ayF  ayD;
%                     axyF axyD]
% Those are the detuning coefficient with maximum current (550A) in the focusing
% octupoles and zero in the defocusing ones (with letter F), or zero in the focusing ones and maximum (550A) 
% in the defocusing octupoles (with letter D).
% components ax are multiplied by Jx and detuning Qx
% components ay are multiplied by Jy and detuning Qy
% components axy are multiplied by Jy (resp. Jx) and detuning Qx (resp. Qy)

% qsec is a matrix 2*2: [Q''xF Q''xD;
%                        Q''yF Q''yD]
% Those are the Q'' at 7TeV, with maximum current (550A) in the focusing
% octupoles and zero in the defocusing ones (with letter F), or zero in the focusing ones and maximum (550A)
% in the defocusing octupoles (with letter D).

% name of the twiss file is in twissname
%

ncol=14; % not used but helps for the definitions of name and C below
nlineheader=45; % number of lines in header section
O3=63100; % maximum absolute octupolar strength in T/m^3 (from MAD-X)
mom=7000; % momentum in GeV/c
K3=6*O3/(mom/0.299792458); % K3+ (maximum normalized octupolar strength)

fid = fopen(twissname);
header = textscan(fid,'%s%s%s%f','HeaderLines',4);
fclose(fid);
fid = fopen(twissname);
name = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',nlineheader);
fclose(fid);
fid = fopen(twissname);
C = textscan(fid,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
fclose(fid);

headercolname=header{2};
ringlength=header{4}(find(strcmp(headercolname,'LENGTH'))); % length of the ring
Qx=header{4}(find(strcmp(headercolname,'Q1'))); % tune x
Qy=header{4}(find(strcmp(headercolname,'Q2'))); % tune y

indl=[];
for i=1:ncol+1
  tmp=name{i};
  colname{i}=char(tmp(1));
  if strcmp(colname{i},'L')
      indl=i-1;
  elseif strcmp(colname{i},'S')
      inds=i-1;
  elseif strcmp(colname{i},'BETX')
      indbetax=i-1;
  elseif strcmp(colname{i},'BETY')
      indbetay=i-1;
  elseif strcmp(colname{i},'DX')
      inddx=i-1;
  elseif strcmp(colname{i},'DY')
      inddy=i-1;
  end
end
nameselem=strrep(C{1},'"','');
lengt=C{indl};
s=C{inds};
betax=C{indbetax};
betay=C{indbetay};
dx=C{inddx};
dy=C{inddy};


% plots
%figure;plot(s,betax,'b','LineWidth',2);hold on;plot(s,betay,'r','LineWidth',2);

% select elements
sel='MO'; % element(s) to select
indsel=strmatch(sel,nameselem);
indF=indsel(find(abs(betax(indsel)-180)<10)); % focusing octupoles
indD=indsel(find(abs(betax(indsel)-30)<10)); % defocusing octupoles

% not: additional minus sign for the defocusing octupoles because O3D=-O3F for
% the same current in foc. and defoc. octupoles

axF=sum(lengt(indF).*betax(indF).^2)*K3/(16*pi);
axD=-sum(lengt(indD).*betax(indD).^2)*K3/(16*pi);
ayF=sum(lengt(indF).*betay(indF).^2)*K3/(16*pi);
ayD=-sum(lengt(indD).*betay(indD).^2)*K3/(16*pi);
axyF=-sum(lengt(indF).*betax(indF).*betay(indF))*K3/(8*pi);
axyD=sum(lengt(indD).*betax(indD).*betay(indD))*K3/(8*pi);

a=[axF  axD; ayF  ayD; axyF axyD];

qsecxF=sum(lengt(indF).*betax(indF).*dx(indF).^2)*K3/(4*pi);
qsecxD=-sum(lengt(indD).*betax(indD).*dx(indD).^2)*K3/(4*pi);
qsecyF=-sum(lengt(indF).*betay(indF).*dx(indF).^2)*K3/(4*pi);
qsecyD=sum(lengt(indD).*betay(indD).*dx(indD).^2)*K3/(4*pi);

qsec=[qsecxF qsecxD; qsecyF qsecyD];

    