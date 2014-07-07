function [eigenval,eigenvect]=laclare(Z,taub,Nb,circum,gamma,tune,tunes, ...
    distribution,particle,chroma,alphap,M,mmax,kmax,varargin);

% solve Laclare's eigenvalue problem  in transerse, in the case of low 
% intensity perturbations (no mode coupling), for different kinds of 
% longitudinal distributions.

% This version is actually encapsulating a C++ code.

% It gives in output:
% - eigenval: the eigenvalues Deltaomega_cm (complex betatron 
% angular frequency shift of the modes m), array of size M*(2*mmax+1)*(2*plmax+1)
% (pth eigenvalue of the mth mode of the Mth multibunch mode),
% - eigenvect: the corresponding eigenvectors sigma_(x,m). This is an array
% of size M*(2*mmax+1)*(2*plmax+1)*(2*plmax+1) (lth component of the pth 
% eigenvector of the mth mode of the Mth multibunch mode) (columns are the
% eigenvectors).
%
% Input parameters:
% - Z: transverse impedance (array with three columns: frequency - NOT
% ANGULAR, real part of the impedance and imaginary part). Only the
% positive frequencies are needed,
% - taub: total bunch length in seconds,
% - Nb: number of particles per bunch,
% - circum: circumference of the machine,
% - gamma: relativistic mass factor,
% - tune: transverse betatron tune (integer part + fractional part),
% - tunes: synchrotron tune,
% - distribution: type of longitudinal distribution:
% 'waterbag' (water bag bunch),'parabolicline' (parabolic line density),
% 'parabolicamp' (parabolic amplitude density) (see Laclare)
% - particle: 'proton' or 'electron',
% - chroma: chromaticity (DeltaQ*p/Q*Deltap),
% - alphap: momentum compaction factor,
% - M: number of bunches,
% - mmax: modes considered are from -mmax to mmax,
% - kmax: number of eigenvalues computed for each coupled-bunch (nx) and
% headtail (m) mode, i.e. only the kmax largest are taken (with the
% corresponding eigenvectors)
% - varargin{1}: nxmin= minimum coupled-bunch mode number to consider 
% (0 by default)
% - varargin{2}: nxmax= maximum coupled-bunch mode number to consider 
% (-1 by default). In the end we compute the modes [nxmin:nxmax+M] modulo M

% see Elias Metral's USPAS 2009 course : Bunched beams transverse coherent
% instabilities, and J.L. Laclare lectures at CAS (1987, p.264)

path=evalc('!pwd');

if (length(strfind(path,'AFS_sync'))>=1)
    pathLacl='/home/nicolasmounet/Documents/Desktop_HD_sync/soft/My_prog_C/Laclare';
else
    pathLacl='/home/nmounet/Documents/soft/My_prog_C/Laclare';
end

eigenval=zeros(M,2*mmax+1,kmax);

if (length(varargin)>0)&&(varargin{1}>0)
	nxmin=varargin{1};
else
	nxmin=0;
end
if (length(varargin)>1)
	nxmax=varargin{2};
else
	nxmax=-1;
end

fid=fopen('input.dat','wt');
fprintf(fid,'Impedance filename\tZ.dat\n');
fprintf(fid,'Output filename\tout\n');
fprintf(fid,'Total bunch length [seconds]\t%13.8e\n',taub);
fprintf(fid,'Number of particles per bunch\t%13.8e\n',Nb);
fprintf(fid,'Machine circumference [m]\t%13.8e\n',circum);
fprintf(fid,'Relativistic gamma\t%13.8e\n',gamma);
fprintf(fid,'Transverse tune\t%13.8e\n',tune);
fprintf(fid,'Synchrotron tune\t%13.8e\n',tunes);
fprintf(fid,'Longitudinal distribution\t%s\n',distribution);
fprintf(fid,'Type of particle\t%s\n',particle);
fprintf(fid,'Chromaticity (DeltaQ*p/Q*Deltap)\t%13.8e\n',chroma);
fprintf(fid,'Momentum compaction factor\t%13.8e\n',alphap);
fprintf(fid,'Number of bunches\t%d\n',M);
fprintf(fid,'Maximum headtail mode considered\t%d\n',mmax);
fprintf(fid,'Maximum number of eigenvalues\t%d\n',kmax);
fprintf(fid,'Minimum coupled-bunch mode number to consider\t%d\n',nxmin);
fprintf(fid,'Maximum coupled-bunch mode number to consider\t%d\n',nxmax);


fclose(fid);

dlmwrite('Z.dat',Z,'precision','%20.15g','Delimiter',' ');
eval(['!cp input.dat ',pathLacl]);
eval(['!mv Z.dat ',pathLacl]);
eval(['cd ',pathLacl]);
!./laclare.x < input.dat > out
eval(['!rm input.dat Z.dat']);

plmax0=-1;
fid1=fopen('out_val.dat','r');
fid2=fopen('out_vec.dat','r');
for nxi=nxmin:M+nxmax
    nx=mod(nxi,M);
    for m=-mmax:mmax
        for k=1:kmax
            x=fscanf(fid1,'%lf',1);y=fscanf(fid1,'%lf\n',1);
            eigenval(nx+1,mmax+m+1,k)=x+1j*y;
            plmax=fscanf(fid2,'%d\n',1);
            plmax0=max(plmax,plmax0);
            for l=-plmax:plmax
                x=fscanf(fid2,'%lf ',1);y=fscanf(fid2,'%lf\n',1);
            end
            fscanf(fid2,'\n',1);
        end
        fscanf(fid1,'\n',1);fscanf(fid2,'\n',1);
    end
    fscanf(fid1,'\n',1);fscanf(fid2,'\n',1);
end

fclose(fid1);
fclose(fid2);

eigenvect=zeros(M,2*mmax+1,2*plmax0+1,kmax);
fid2=fopen('out_vec.dat','r');

for nxi=nxmin:M+nxmax
    nx=mod(nxi,M);
    for m=-mmax:mmax
        for k=1:kmax
            plmax=fscanf(fid2,'%d\n',1);
            for l=-plmax:plmax
                x=fscanf(fid2,'%lf ',1);y=fscanf(fid2,'%lf\n',1);
                eigenvect(nx+1,mmax+m+1,plmax0+l+1,k)=x+1j*y;
            end
            fscanf(fid2,'\n',1);
        end
        fscanf(fid2,'\n',1);
    end
    fscanf(fid2,'\n',1);
end

fclose(fid2);

%eval(['!rm out out_val.dat out_vec.dat']);
try eval(['cd ''',path,'''']);
catch ME1
    path1=path;
    path1(1)=''
    try eval(['cd ''',path1,'''']);
    catch ME1
        path(end)=''
        eval(['cd ''',path,'''']);
    end
end
