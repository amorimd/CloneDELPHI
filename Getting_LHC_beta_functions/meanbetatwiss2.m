function [avbetax,avbetay]=meanbetatwiss2(twissname,sbeg,send);

% computes mean beta functions of some elements, from a twiss file (madx)
% and their coordinates s (beginning in sbeg, end in send)
%
% name of the twiss file is in twissname
%
% it uses a spline interpolation (every 0.01m) and a trapezoidal integration

ncol=12; % number of columns in twiss file. if changed, change also lines 17, 20 and 43 accordingly
nlineheader=45; % number of lines in header section

fid = fopen(twissname);
header = textscan(fid,'%s%s%s%f','HeaderLines',4);
fclose(fid);
fid = fopen(twissname);
name = textscan(fid,'%s%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',nlineheader);
fclose(fid);
fid = fopen(twissname);
C = textscan(fid,'%s%s%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
fclose(fid);

headercolname=header{2};
ringlength=header{4}(find(strcmp(headercolname,'LENGTH'))); % length of the ring

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
  end
end

if (inds==2)
    fid = fopen(twissname);
    C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
    fclose(fid);
end
s=C{inds};
betax=C{indbetax};
betay=C{indbetay};

if (length(indl)~=0)
    lengt=C{indl};
else
    lengt=diff([0;C{inds}]);
end


% delete duplicates
ind=find(diff(s)==0);s2=s;betax2=betax;betay2=betay;
s2(ind)=[];betax2(ind)=[];betay2(ind)=[];

% spline interpolation
% sampling
sint=sort([[0:0.01:ringlength] s2.' sbeg send]);sint(find(diff(sint==0)))=[];
% selecting only the points between sbeg and send
indsel=[];betaxint=[];betayint=[];
for i=1:length(sbeg)
    % select only the points between sbeg and send
    indtmp=find((sint-sbeg(i)).*(sint-send(i))<=0);
    indsel=[indsel indtmp];
    % interpolate
    betaxtmp=interp1(s2,betax2,sint(indtmp),'spline');
    betaytmp=interp1(s2,betay2,sint(indtmp),'spline');
    betaxint=[betaxint betaxtmp];
    betayint=[betayint betaytmp];
    % computes the means on that section
    avbetaxi(i)=trapz(sint(indtmp),betaxtmp);
    avbetayi(i)=trapz(sint(indtmp),betaytmp);
end
avbetax=sum(avbetaxi)/sum(send-sbeg);
avbetay=sum(avbetayi)/sum(send-sbeg);
disp(['av. betax and betay with trapz method and interpolation: ',num2str(avbetax),' ',num2str(avbetay)]);

% plot
%figure;plot(sint(indsel),betaxint,'xb','LineWidth',2);hold on;
%plot(sint(indsel),betayint,'xr','LineWidth',2);
