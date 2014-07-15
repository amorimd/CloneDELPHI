#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import numpy as np
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from particle_param import *
from Impedance import *
from DELPHI import *

def VEPP_param(E0):

    ''' generates typical VEPP parameters at injection, given the electron rest energy 
    E0 in J (from e.g. function proton_param), and generate dipolar impedance (from
    single broad-band resonator). See Karliner-Popov paper (2005).
    Outputs:
    - machine: string with machine name ('VEPP'),
    - E: beam injection energy in eV (1.8 GeV),
    - gamma: relativistic mass factor,
    - sigmaz: RMS bunch length in m,
    - taub: total bunch length in s (4*RMS),
    - R: machine pysical radius (circumference/(2 pi)),
    - Qx: total horizontal tune (integer + fractional parts),
    - Qxfrac: fractional horizontal tune,
    - Qs: synchrotron tune,
    - eta: slippage factor (alpha_p-1/gamma^2),
    - M: number of bunches (1),
    - f0: revolution frequency,
    - omega0: revolution angular frequency=2pi*f0,
    - omegas: synchrotron angular frequency=Qs*omega0,
    - dphase: phase of damper w.r.t. "normal" purely resistive damper,
    - nx: coupled-bunch mode number (0),
    - R1: shunt impedance of broad-band resonator (Ohm/m) in the impedance model
    (including beta function weight),
    - Zx: horizontal dipolar impedance (funcion of frequencies),
    - f: frequencies corresponding to impedance,
    - model: name of impedance model.
    '''

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='VEPP';
    E=1.8e9; # injection energy=1.8 GeV
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    sigmaz=7.5e-2; # RMS bunch length (m)
    taub=4.*sigmaz/(beta*c); # full length in s
    circ=beta*c/8.377e5; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    Qx=7.62;Qxfrac=Qx-np.floor(Qx);
    Qs=0.025;
    alphap=0.01645; # momentum compaction factor
    M=1; # number of bunches
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    dphase=0.; # additional damper phase
    nx=0;

    # impedance definition: one broad band model (NB: shunt impedances to be multiplied by beta(location)/(R/Q) )
    beta1=15;betaav=R/Qx
    model='_1BB_Karliner-Popov';
    f=np.concatenate((10.**np.arange(-1,7),np.arange(1.e7,1.001e10,1.e7),10.**np.arange(10.1,13.1,0.1),
    	10.**np.arange(14,16)));
    R1=2.5e6*beta1/betaav;f1=R*f0*0.795/sigmaz;Q1=1.;
    Zx=resonator_impedance(R1,f1,Q1,f);

    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model;
    

def VEPP_damper(R,Qx,f0):

    ''' VEPP-4 damper model
    - R: machine (physical) radius, 
    - Qx: total tune (integer+fractional parts),
    - f0: revolution frequency,
    Outputs: "damper impedance" Zd, vs frequencies fd
    '''
    
    c=299792458;
    fd=np.concatenate((np.arange(2.e4,2.002e7,2.e4),np.arange(2.1e7,1.001e9,1e6),np.arange(1.1e9,1.01e10,1.e8)));
    fd=np.concatenate((np.flipud(-fd),np.array([0.]),fd));
    Zd=np.zeros((len(fd),2));
    #print fd,len(fd)
    L0=0.4;L1=0.1;L2=0.2;tauf=800*L1/c; # set of parameters with bare tune (including integer part)
    #print "Damper model: tauf*omega0=",tauf*2*np.pi*f0;
    dimp=damper_imp_Karliner_Popov(L0,L1,L2,tauf,R,Qx,f0,fd);
    Zd[:,0]=np.real(dimp);Zd[:,1]=np.imag(dimp);
    
    return Zd,fd;



if __name__ == "__main__":

    e,m0,c,E0=electron_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model=VEPP_param(E0);
    
    # VEPP-4 damper model
    Zd,fd=VEPP_damper(R,Qx,f0);
    strfreqflag=['','_freqdepgain'];
    strfreq=['',', freq. dependent gain'];
    strnorm=['','_norm_current_chroma'];
    
    lmax=2;nmax=1;
    flagsave=0; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=1; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    
    Ib=4*np.pi*E*Qx*Qs*sigmaz/(R**2*R1) # some normalization factor
    print "Ib=",Ib;

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # normalization factor for damper
    dnormfactor=compute_damper_matrix(0,0,nx,M,0.,omega0,Qxfrac,a,b,taub,g,
    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0]; # this is equal to 1

    # scan definition
    Qpxscan=np.array([0.,25,-7.5]);
    
    for iQp,Qpx in enumerate(Qpxscan):
    
    	omegaksi=Qpx*omega0/eta;
	
	if flagnorm:
	    # normalization factor for damper
	    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
	    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];
	    
	print "dnormfactor=",dnormfactor
	
	matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
		flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
	
	matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
		Zx,f,flag_trapz=1,abseps=1e-4);

	if (Qpx==0):
	    fdampscan=np.array([0.,2.5]);
	    Iscan=np.arange(0.2,20.2,0.2); # intensity scan in mA
	elif (Qpx==25):
	    fdampscan=np.array([0.]);
	    Iscan=np.arange(0.5,50.5,0.5); # intensity scan in mA
	elif (Qpx==-7.5):
	    fdampscan=np.array([2.5]);
	    Iscan=np.arange(0.5,50.5,0.5); # intensity scan in mA
	
     	Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch

    
    	for idamp,fdamp in enumerate(fdampscan):
	
	    # Karliner-Popov 2005 results
	    filename=path_here+'Karliner_Popov_figures/Karliner_Popov_'+machine+'_Qp'+float_to_str(Qpx)+'_f'+float_to_str(fdamp)+'_curve';

	    # output file name
	    os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
	    fileout=path_here+'../../../DELPHI_results/'+machine+'/plot_TMCI_DELPHI_vs_Karliner-Popov_'+machine+'_'+float_to_str(round(E/1e8)/10)+'GeV'+model+'_'+str(M)+'b_f'+float_to_str(fdamp)+'_'+str(nmax+1)+'radialmodes_maxazim'+str(lmax)+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];

	    dscan=Iscan*1.e-3*fdamp*2*np.pi*Qs/Ib; # damper gain scan (gain depends on intensity here)
	    print dscan
	    kmax=(2*lmax+1)*(nmax+1);
	    lambdax=np.zeros((len(Nbscan),kmax),dtype=complex);
	    
	    # computation
	    for iNb,Nb in enumerate(Nbscan):
	    
	    	coefdamper,coefZ=computes_coef(f0,dscan[iNb],b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='electron');
		
		freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas)
		
		lambdax[iNb,:]=freqshift/omegas;		
    	    
	    
	    # TMCI plots
	    strpart=['Re','Im'];
	    for ir,r in enumerate(['real','imag']):
		
		fig,ax=init_figure();
		plot_TMCI(Iscan,lambdax,ax,part=r,leg='DELPHI',patcol='xb',xlab='Intensity [mA]',
			title=machine+", $ Q^' = $ "+str(Qpx)+', $ f= $'+str(fdamp)+', '+str(nmax+1)+' radial modes, '+str(2*lmax+1)+' azimuthal modes'+strfreq[flagdamperimp]);
		
		# Karliner-Popov results
		maxs=0.;
		for k in range(1,11):
	    	    s=read_ncol_file(filename+str(k)+'_'+r+'.csv',ignored_rows=1);
		    if (k==1): leg='Karliner-Popov';
		    else: leg='';
	    	    plot(s[:,0],s[:,1],leg,'r',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",ax,0,xlab='Intensity [mA]');
		    maxs=max(maxs,np.max(s[:,1]));
		
		if ir==1: ax.set_ylim([np.floor(np.min(np.imag(lambdax.flatten()))/0.2)*0.2,np.ceil(maxs/0.2)*0.2]);

		if flagsave: end_figure(fig,ax,save=fileout+'_'+r)
		else: end_figure(fig,ax);


    if not(flagsave): pylab.show();
