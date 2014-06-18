function [a,qsec]=octupole_special(twissname,twissname_ref,varargin);

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
% Those are the Q'' with maximum current (550A) in the focusing
% octupoles and zero in the defocusing ones (with letter F), or zero in the focusing ones and maximum (550A)
% in the defocusing octupoles (with letter D).

% name of the twiss file is in twissname
%
% special case of HL-LHC twiss file from S. Fartoukh:
% twissname_ref contains an LHC twiss file (this is to disentangle foc. /
% defoc. octupoles)

% varagin{1}: number of columns in optics file (only 2 cases possible : 10
% or 14)

if (length(varargin)>0)
    ncol=varargin{1} % number of columns in tiwss file. if changed, change lines 46 and 49 accordingly
else
    ncol=10; % number of columns in tiwss file. if changed, change lines 51 and 54 accordingly
end
nlineheader=45; % number of lines in header section
O3=63100; % maximum absolute octupolar strength in T/m^3 (from MAD-X)
mom=7000; % momentum in GeV/c
K3=6*O3/(mom/0.299792458); % K3+ (maximum normalized octupolar strength)

fid = fopen(twissname);
header = textscan(fid,'%s%s%s%f','HeaderLines',4);
fclose(fid);
fid = fopen(twissname);
if ncol==10
    name = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',nlineheader);
    fclose(fid);
    fid = fopen(twissname);
    C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
elseif ncol==14
    name = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',nlineheader);
    fclose(fid);
    fid = fopen(twissname);
    C = textscan(fid,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
end    
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
  end
end
nameselem=strrep(C{1},'"','');
if (length(indl)~=0)
    lengt=C{indl};
else
    lengt=diff([0;C{inds}]);
end
s=C{inds};
betax=C{indbetax};
betay=C{indbetay};
dx=C{inddx};

% read now reference file
fid = fopen(twissname_ref);
header_ref = textscan(fid,'%s%s%s%f','HeaderLines',4);
fclose(fid);
fid = fopen(twissname_ref);
name_ref = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',nlineheader);
fclose(fid);
fid = fopen(twissname_ref);
C_ref = textscan(fid,'%s%s%f%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
fclose(fid);
nameselem_ref=strrep(C_ref{1},'"','');
for i=1:15
    tmp_ref=name_ref{i};
    colname_ref{i}=char(tmp_ref(1));
    if strcmp(colname_ref{i},'BETX')
        indbetax_ref=i-1;
    end
end
betax_ref=C_ref{indbetax_ref};

% plots
%figure;plot(s,betax,'b','LineWidth',2);hold on;plot(s,betay,'r','LineWidth',2);

% select elements
sel='MO'; % element(s) to select
indsel_ref=strmatch(sel,nameselem_ref);
indF_ref=indsel_ref(find(abs(betax_ref(indsel_ref)-180)<10)); % focusing octupoles
indD_ref=indsel_ref(find(abs(betax_ref(indsel_ref)-30)<10)); % defocusing octupoles

for i=1:length(indF_ref)
    indF(i)=strmatch(nameselem_ref(indF_ref(i)),nameselem);
end
for i=1:length(indD_ref)
    indD(i)=strmatch(nameselem_ref(indD_ref(i)),nameselem);
end
%length(indF),length(indD)

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

    