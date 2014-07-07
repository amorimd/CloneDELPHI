function [tau,signal]=signallaclare(eigenval,eigenvect,taub,ntau,Nb,circum,gamma,tune,tunes,M,m,nx,nturns,varargin);

% computes and plots the signal from Laclare's eigenvalues and eigenvectors (in transerse),
% n the case of low intensity perturbations (no mode coupling).
% It gives in output:
% - signal: signal that could be observed in a machine at azimuth (along
% the ring) theta=0. It's an array of nturns*ntau values.
% - tau: ntau values linearly spaced between -taub/2 and taub/2.
%
% Input parameters:
% - eigenval (scalar) and eigenvect (column of 2*plmax+1 elements): pth eigenvalue
% and corresponding eigenvector for a certain multibunch mode nx and a certain mode m
% - taub: total bunch length in seconds,
% - ntau: number of tau points considered,
% - Nb: number of particles per bunch,
% - circum: circumference of the machine,
% - gamma: relativistic mass factor,
% - tune: transverse betatron tune (integer part + fractional part),
% - tunes: synchrotron tune,
% - M: number of bunches,
% - m: mode m considered,
% - nx: multibunch mode considered,
% - nturns: number of turns considered,
% - varargin{1}: Mmax: number of bunches considered (deactivated if 0),
% - varargin{2} and varargin{3}: if present, respectively other
% eigenvalue-eigenvector pairs, to compare on the same
% plots. It can be either scalar and column vector respectively (as the two
% first input parameters) or cell-array with several such pairs.
% NOTE: the output "signal" is not changed (only with first 
% eigenvalue-eigenvector pair). Deactivated if [],
% - varargin{4-6}: if varargin{4}=1, plots also the signal from the fit 
% functions of Sacherer. It needs then chromaticity in varargin{5}, and
% slip factor eta in varargin{6}.

% see Elias Metral's USPAS 2009 course : Bunched beams transverse coherent
% instabilities, and J.L. Laclare lectures at CAS (1987, p.264)

col={'-xb','r','g','m','k','c','y'}; % color order for plot
e=1.60218e-19; % elementary charge in C
clight=299792458;% speed of light in m/s
beta=sqrt(1-1/gamma^2); % relativistic velocity factor
f0=beta*clight/circum; % revolution frequency
omega0=2*pi*f0; % revolution angular frequency
Ib=M*e*Nb*f0; % bunch intensity
omegas=tunes*omega0; % synchrotron angular frequency

% sample tau
tau=linspace(-taub/2,taub/2,ntau);

figS=figure;
% computes signal
if (length(varargin)>0)&&(varargin{1}>0)
    Mmax=varargin{1};
else
    Mmax=M;
end


% first signal
signal=zeros(nturns,Mmax,ntau);
plmax=(length(eigenvect)-1)/2;
omegac=eigenval+m*omegas; % coherent angular frequency of oscillation (complex)
% (actually, delta from tune*omega0)

for n=0:nturns-1
    %for Mb=0:floor(M/(nx+1)):Mmax-1
    for Mb=0:Mmax-1
        % for multibunch
        t=tau+n/f0+Mb/(M*f0); % times considered for this turn and this bunch
        signaltmp=zeros(1,length(tau));
        for p=-plmax:plmax
            signaltmp=signaltmp+eigenvect(p+plmax+1)*exp(1j*(nx+M*p)*omega0*t);
        end
        signaltmp=2*pi^2*Ib*signaltmp.*exp(1j*((tune-floor(tune))*omega0+omegac)*t);
        %plot(tau,real(signaltmp),'Color',[0 n/nturns Mb/Mmax],'LineWidth',2);hold on;
        % to plot all traces on same bunch
        %plot(tau,real(signaltmp),col{1},'LineWidth',2);hold on;
        % to plot the whole bunch train
        plot(tau+Mb/(M*f0),real(signaltmp),col{1},'LineWidth',2);hold on;
        % to plot on several turns
        %plot(t,real(signaltmp),col{1},'LineWidth',2);hold on;
        signal(n+1,Mb+1,:)=signaltmp.';
    end
end
%disp(['Max. relative imaginary part of signal: ', ...
%    num2str(max( (abs(imag(signal(:)))-abs(real(signal(:))) )./abs(real(signal(:)))))]);

% handle plot
xlabel('\tau (s)','FontSize',20);ylabel('Pick up signal to be measured (A.m)','FontSize',20);
grid on;set(figS,'Color',[1 1 1],'PaperType','A4','PaperOrientation','Landscape');
set(gca,'LineWidth',2,'FontSize',16,'YMinorgrid','off');

% plot also the eigenvector(s) considered as a function of frequency
% real part
figsigR=figure;
freq=[(nx-M*plmax):M:(nx+M*plmax)]*f0+(tune-floor(tune))*f0;
plot(freq,real(eigenvect.'),col{1},'LineWidth',2);hold on;
xlabel('Frequency (Hz)','FontSize',20);ylabel('Eigenvector \sigma [m] (real part)','FontSize',20);
grid on;set(figsigR,'Color',[1 1 1],'PaperType','A4','PaperOrientation','Landscape');
set(gca,'LineWidth',2,'FontSize',16,'YMinorgrid','off');
% imag. part
figsigI=figure;plot(freq,imag(eigenvect.'),col{1},'LineWidth',2);hold on;
xlabel('Frequency (Hz)','FontSize',20);ylabel('Eigenvector \sigma [m] (imag. part)','FontSize',20);
grid on;set(figsigI,'Color',[1 1 1],'PaperType','A4','PaperOrientation','Landscape');
set(gca,'LineWidth',2,'FontSize',16,'YMinorgrid','off');

% second (or more) signals (if other eigenvalue-eigenvector pair(s) is present 
% in varargin{2-3})
neig=0;
if (length(varargin)>2)&&(length(varargin{2})>0)
    eigenvaltmp=varargin{2};
    if iscell(eigenvaltmp)
        neig=length(eigenvaltmp);
        eigenval2=eigenvaltmp;
        eigenvect2=varargin{3};       
    else
        neig=1;
        eigenval2{1}=eigenvaltmp;
        eigenvecttmp=varargin{3};
        eigenvect2{1}=eigenvecttmp;
    end
    for eig=1:neig
        plmax=(length(eigenvect2{eig})-1)/2;
        omegac=eigenval2{eig}+m*omegas; % coherent angular frequency of oscillation (complex)
        % (actually, delta from tune*omega0)
        figure(figS);
        for n=0:nturns-1
            %for Mb=0:floor(M/(nx+1)):Mmax-1
            for Mb=0:Mmax-1
                % for multibunch
                t=tau+n/f0+Mb/(M*f0); % times considered for this turn and this bunch
                signaltmp=zeros(1,length(tau));
                for p=-plmax:plmax
                    signaltmp=signaltmp+eigenvect2{eig}(p+plmax+1)*exp(1j*(nx+M*p)*omega0*t);
                end
                signaltmp=2*pi^2*Ib*signaltmp.*exp(1j*((tune-floor(tune))*omega0+omegac)*t);
                %plot(tau,real(signaltmp),'Color',[0 n/nturns Mb/Mmax],'LineWidth',2);hold on;
                plot(tau,real(signaltmp),col{eig+1},'LineWidth',2);hold on;
            end
        end
        figure(figsigR);
        freq=[-(nx+M*plmax):M:(nx+M*plmax)]*f0+(tune-floor(tune))*f0;
        plot(freq,real(eigenvect2{eig}.'),col{eig+1},'LineWidth',2);hold on;
        figure(figsigI);
        plot(freq,imag(eigenvect2{eig}.'),col{eig+1},'LineWidth',2);hold on;
    end
end

if (length(varargin)>5)&&(varargin{4}==1)
    % plots the signal from the fit functions of Sacherer
    % see Elias Metral's USPAS 2009 course on bunched beam instabilities
    figure(figS);
    chroma=varargin{5};eta=varargin{6};
    omegaksi=chroma*tune*omega0/eta;
    pmtau=((tau>=-taub/2)&(tau<=taub/2)).*(cos((abs(m)+1)*pi*tau/taub)*mod(m+1,2) + ...
        sin((abs(m)+1)*pi*tau/taub)*mod(m,2));
    for n=0:nturns-1
        %for Mb=0:floor(M/(nx+1)):Mmax-1
        for Mb=0:Mmax-1
            % for multibunch
            t=tau;%+n/f0+Mb/(M*f0); % times considered for this turn and this bunch
            signaltmp=pmtau.*exp(1j*(2*pi*(n+Mb/M)*tune+omegaksi*t));
            % point at center in first signal plotted
            sigcenter=signal(n+1,Mb+1,floor(ntau/2));
            % point at center of signaltmp here
            sigtmpcenter=signaltmp(floor(ntau/2));
            % normalizing signaltmp
            signaltmp=signaltmp*(sigcenter/sigtmpcenter);
            %signaltmp=signaltmp*(2*pi^2*Ib);
            plot(tau,real(signaltmp),col{neig+2},'LineWidth',2);hold on;
        end
    end
end
    

% plot also the bunch(es) centroid versus time
% not useful: you see the betatron oscillation basically
% figcen=figure;k=1;clear time signalcen;
% for n=0:nturns-1
%    for Mb=0:M-1
%        time(k)=n/f0+Mb/(M*f0);
%        t=tau+n/f0+Mb/(M*f0);
%        signaltmp=zeros(1,length(tau));
%        for p=-plmax:plmax
%            signaltmp=signaltmp+eigenvect(p+plmax+1)*exp(1j*(nx+M*p)*omega0*t);
%        end
%        signaltmp=2*pi^2*Ib*signaltmp.*exp(1j*(tune*omega0+omegac)*t);
%        signalcen(k)=mean(signaltmp);
%        k=k+1;
%    end
% end
% plot(time,real(signalcen),'-xb','LineWidth',2);
% xlabel('time (s)','FontSize',20);ylabel('Pick up signal (bunch centroid) (A.m)','FontSize',20);
% grid on;set(figcen,'Color',[1 1 1],'PaperType','A4','PaperOrientation','Landscape');
% set(gca,'LineWidth',2,'FontSize',16,'YMinorgrid','off');

