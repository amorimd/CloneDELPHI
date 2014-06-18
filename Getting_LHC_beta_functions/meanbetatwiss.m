function [s,betax,betay,avbetax,avbetay,varargout]=meanbetatwiss(twissname,varargin);

% computes mean beta functions, from a twiss file (madx)
% and plots the beta functions
% in output: s coordinate along the ring, beta functions and their averages
%
% name of the twiss file is in twissname
%
% varargin: optional arguments:
% - varargin{1} : name of the element (or part of their name) to extract
% then in output you get s, betax and betay for those elements (average)
% if not present or '', it gives all the beta along the ring.
% then varargout{1} contains the names of those elements, varargout{2}
% their lengths, varargout{3}=sbeg (longitudinal coordinate at the begining
% of each element), varargout{4}=send (longitudinal coordinate at the endg
% of each element)
% -varargin{2} : flag : if present and 1, varargin{2} is a regular
% expression with *, . and things like that (see help of regexp function)

nlineheader=45; % number of lines in header section

% find number of columns
data=dlmread(twissname,'',nlineheader+2,2);
ncol=length(data(1,:))+2;disp(['ncol=',num2str(ncol)]);
format_name_mon='%s';format_name='%s'; % one more "%s" (due to * character)
for i=1:ncol
    format_name=[format_name,format_name_mon];
end
format_col_mon='%f';format_col='%s%s'; % first two columns are strings
for i=1:ncol-2
    format_col=[format_col,format_col_mon];
end

fid = fopen(twissname);
header = textscan(fid,'%s%s%s%f','HeaderLines',4);
fclose(fid);
fid = fopen(twissname);
name = textscan(fid,format_name,'HeaderLines',nlineheader);
fclose(fid);
fid = fopen(twissname);
C = textscan(fid,format_col,'HeaderLines',nlineheader+2);
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
  elseif strcmp(colname{i},'MUX')
      indmux=i-1;
  elseif strcmp(colname{i},'MUY')
      indmuy=i-1;
  elseif strcmp(colname{i},'DX')
      inddx=i-1;
  end
end
nameselem=strrep(C{1},'"','');
if (length(indl)~=0)
    lengt=C{indl};
end
if (inds==2)
    fid = fopen(twissname);
    C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
    fclose(fid);
end
s=C{inds};
betax=C{indbetax};
betay=C{indbetay};

avbetax1=ringlength/(2*pi*Qx);
avbetay1=ringlength/(2*pi*Qy);
disp(['av. betax and betay with smooth approx: ',num2str(avbetax1),' ',num2str(avbetay1)]);
% delete duplicates
ind=find(diff(s)==0);s2=s;betax2=betax;betay2=betay;
s2(ind)=[];betax2(ind)=[];betay2(ind)=[];

% the above is a bad way to compute the mean beta function. better to do:
%if (length(indl)~=0)
%    avbetax2=sum(betax.*lengt)/ringlength;
%    avbetay2=sum(betay.*lengt)/ringlength;
%    disp(['av. betax and betay without trapz method: ',num2str(avbetax2),' ',num2str(avbetay2)]);
%end
% finally, use method of trapz
%avbetax=trapz(s,betax)/ringlength;
%avbetay=trapz(s,betay)/ringlength;
%disp(['av. betax and betay with trapz method: ',num2str(avbetax),' ',num2str(avbetay)]);

% spline interpolation
sint=sort([[0:0.01:ringlength] s2.']);sint(find(diff(sint==0)))=[];
% interpolate
betaxint=interp1(s2,betax2,sint,'spline');
betayint=interp1(s2,betay2,sint,'spline');
%recompute the mean
avbetax=trapz(sint,betaxint)/ringlength;
avbetay=trapz(sint,betayint)/ringlength;
disp(['av. betax and betay with trapz method and interpolation: ',num2str(avbetax),' ',num2str(avbetay)]);

% note: we should find back (more or less) the tunes when doing:
disp(['real tune x: ',num2str(Qx)]);
disp(['real tune y: ',num2str(Qy)]);
%if (length(indl)~=0)
%    disp(['tune x from beta x without trapz method: ',num2str(sum(lengt./betax)/2/pi)]);
%    disp(['tune y from beta y without trapz method: ',num2str(sum(lengt./betay)/2/pi)]);
%end
%disp(['tune x from beta x with trapz method: ',num2str(trapz(s,1./betax)/2/pi)]);
%disp(['tune y from beta y with trapz method: ',num2str(trapz(s,1./betay)/2/pi)]);
disp(['tune x from beta x with trapz method and interpolation: ',num2str(trapz(sint,1./betaxint)/2/pi)]);
disp(['tune y from beta y with trapz method and interpolation: ',num2str(trapz(sint,1./betayint)/2/pi)]);

% plots
%figure;plot(s,betax,'b','LineWidth',2);hold on;plot(s,betay,'r','LineWidth',2);

% select elements
if (length(varargin)>0 ) && (length(varargin{1})>0)
    sel=varargin{1}; % element(s) to select
    if (length(varargin)>1) && (varargin{2})
        % case when sel is a "regular expression"
        a=regexp(nameselem,sel,'once');k=1;
        for j=1:length(nameselem)
          if (length(a{j})>0)&&(a{j}==1);
              indsel(k,1)=j;k=k+1;
          end
        end
    else
        indsel=strmatch(sel,nameselem);
    end
    %s=(s(indsel)+s(indsel-1))/2;
    %betax=(betax(indsel)+betax(indsel-1))/2;
    %betay=(betay(indsel)+betay(indsel-1))/2;
    varargout{1}=nameselem(indsel);
    if (length(indl)~=0)
        len=lengt(indsel);
        varargout{2}=len;
    end
    clear betax3 betay3 sbeg send;
    for i=indsel.'
        if (s(i-1)<s(i))
            % spline interpolation
            ind1=find(sint==s(i-1));ind2=find(sint==s(i));
            sint2=sint(ind1:ind2);
            betaxint2=betaxint(ind1:ind2);
            betayint2=betayint(ind1:ind2);
            % compute the mean using the spline interpolation
            betax3(find(indsel==i))=trapz(sint2,betaxint2)/lengt(i);
            betay3(find(indsel==i))=trapz(sint2,betayint2)/lengt(i);
        else
            betax3(find(indsel==i))=(betax(i)+betax(i-1))/2;
            betay3(find(indsel==i))=(betay(i)+betay(i-1))/2;
        end
        if (lengt(i)~=len(find(indsel==i)))
            disp (['Big Problem in the element length, ',num2str(i),', ',num2str(lengt(i)),', ',num2str(len(i))]);
            lengt(i),len(find(indsel==i))
        end
        if (abs(lengt(i)-(s(i)-s(i-1)))>1e-6)
            disp (['Problem in the element length, ',num2str(i),', ',num2str(lengt(i)),', ',num2str(s(i)-s(i-1))]);
            lengt(i),s(i)-s(i-1)
        end
    end
    sbeg=s(indsel-1).';
    varargout{3}=sbeg; % beginning of the selected elements
    send=s(indsel).';
    s=(s(indsel)+s(indsel-1))/2;
    varargout{4}=send; % end of the selected elements
    betax=betax3.';betay=betay3.';
    avbetax=sum(betax.*len)/sum(len);
    avbetay=sum(betay.*len)/sum(len);
    disp(['av. betax and betay in the selected elements, with trapz method and interpolation: ',num2str(avbetax),' ',num2str(avbetay)]);
end

if (nargout>7 )
    % extract phases as well
    mux=C{indmux};
    muy=C{indmuy};
    varargout{5}=mux;
    varargout{6}=muy;
end

if (nargout>11 )
    % extract dispersion as well
    dx=C{inddx};
    varargout{7}=dx;
end
