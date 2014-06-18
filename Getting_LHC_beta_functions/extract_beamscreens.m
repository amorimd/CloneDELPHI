function [sbegin,send,varargout]=extract_beamscreens(twissname,position,BSlength);

% extract beginning and end positions of the beam
% screens, given their position in terms of element around which they are
% around the ring, and their length. we assume the beam screen are centered
% around the element given in position.
% In input:
% - twissname: name of the twiss file used,
% - position: string giving a regular expression to find the elements
% around which the beam screens are. If position is a cell array of strings
% (for instance {'MQM','MQTLI'}) then the center of the beam screen section
% is at the center of the two elements given.
% - length : length (in m) of each individual beam screen.
%
% In output: s coordinate along the ring for the beginnings of each beam
% screen section (sbegin), and for their end (send).
% When present, varargout{1} contains the names of the elements selected
% and varargout{2} their lengths.
%

ncol=12; % number of columns in twiss file. if changed, change also lines 28, 31 and 50 accordingly
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
  end
end
nameselem=strrep(C{1},'"','');
if (inds==2)
    fid = fopen(twissname);
    C = textscan(fid,'%s%f%f%f%f%f%f%f%f%f%f%f','HeaderLines',nlineheader+2);
    fclose(fid);
end
s=C{inds};
if (length(indl)~=0)
    lengt=C{indl};
else
    lengt=diff([0;C{inds}]);
end


% select elements
if not(iscell(position))
    position2{1}=position;
    clear position;
    position=position2;
end
for i=1:length(position)
    a=regexp(nameselem,position{i},'once');k=1;
    for j=1:length(nameselem)
        if (length(a{j})>0)&&(a{j}==1);
            indsel{i}(k)=j;k=k+1;
        end
    end
end
% if there are several elements, sort according to their positions
if (length(position)>1)
    clear tmp
    for j=1:length(indsel{1})
        for i=1:length(position)
            tmp(i)=indsel{i}(j);
        end
        [tmp2,ind2]=sort(tmp);
        for i=1:length(position)
            indsel2{i}(j)=indsel{ind2(i)}(j);
        end        
    end
    indsel=indsel2;
end
scenter=(s(indsel{1}-1)+s(indsel{end}))/2;
sbegin=scenter.'-BSlength/2;
send=scenter.'+BSlength/2;
% test that scenter is not far away from the elements (could happen if
% there are several strings in position and if the matching did not find the elements
% around the same Q or in the same IR)
if ( abs(max(-scenter+s(indsel{end})))>BSlength )||( abs(max(scenter-s(indsel{1}-1)))>BSlength )
    disp('Pb with s');
end

varargout{1}=[];varargout{2}=[];
for j=1:length(indsel{1})
    tmp=[];tmpl=[];
    for i=1:length(position)
        tmp=[tmp nameselem(indsel{i}(j))];
        tmpl=[tmpl lengt(indsel{i}(j))];
    end
    varargout{1}=[varargout{1} ; tmp];
    varargout{2}=[varargout{2} ; tmpl];
end

% overlapping test
if (length(sbegin)>1)
    if ( min(sbegin(2:end)-send(1:end-1))<0 ) || ( (send(end)-sbegin(1)-ringlength)>0 )
        disp('Aie... Overlapping beam screens');
    end
end
