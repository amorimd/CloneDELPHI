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
from VEPP import VEPP_damper

def LEP_param(E0):

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='LEP';
    E=22e9; # injection energy=22 GeV
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    sigmaz=1.3e-2; # RMS bunch length (m)
    taub=4.*sigmaz/(beta*c); # full length in s
    circ=26658.883; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    Qx=76.194;Qxfrac=Qx-np.floor(Qx);
    Qs=0.108;
    alphap=1.855e-4; # momentum compaction factor
    M=1; # number of bunches
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    dphase=0.; # additional damper phase
    nx=0;

    # impedance definition: one broad band model (NB: shunt impedances to be multiplied by beta(location)/(R/Q) )
    beta1=40.6;betaav=R/Qx
    model='_2BB_Karliner-Popov';

    f=np.concatenate((10.**np.arange(-1,7),np.arange(5.e7,1.0005e11,5.e7),10.**np.arange(11.1,13.1,0.1),
    	10.**np.arange(14,16)));
    R1=1.51e6*beta1/betaav;f1=R*f0*0.536/sigmaz;Q1=1;f1=2e9;
    R2=0.322e6*beta1/betaav;f2=R*f0*3.22/sigmaz;Q2=1;f2=12e9;
    Zx=resonator_impedance(R1,f1,Q1,f)+resonator_impedance(R2,f2,Q2,f);

    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model;
    


if __name__ == "__main__":

    e,m0,c,E0=electron_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model=LEP_param(E0);
    
    # VEPP-4 damper model
    Zd,fd=VEPP_damper(R,Qx,f0);
    strfreqflag=['','_freqdepgain'];
    strfreq=['',', freq. dependent gain'];
    strnorm=['','_norm_current_chroma'];
    
    lmax=2;nmax=1;
    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    
    Ib=4*np.pi*E*Qx*Qs*sigmaz/(R**2*R1) # some normalization factor

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # normalization factor for damper
    dnormfactor=compute_damper_matrix(0,0,nx,M,0.,omega0,Qxfrac,a,b,taub,g,
    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

    # scan definition
    Qpxscan=np.array([0.,-22.]);
    
    for iQp,Qpx in enumerate(Qpxscan):
    
    	omegaksi=Qpx*omega0/eta;
	
	if flagnorm:
	    # normalization factor for damper
	    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
	    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];
	
	matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
		flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
	
	matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
		Zx,f,flag_trapz=1,abseps=1e-4);

	if (Qpx==0):
	    fdampscan=np.array([0.,2.]);
	    Iscan=np.arange(0.025,1.525,0.025); # intensity scan in mA
	elif (Qpx==-22):
	    fdampscan=np.array([2.]);
	    Iscan=np.arange(0.025,5.025,0.025); # intensity scan in mA
	
     	Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch

    
    	for idamp,fdamp in enumerate(fdampscan):
	
	    # Karliner-Popov results
	    filename=path_here+'Karliner_Popov_figures/Karliner_Popov_Qp'+float_to_str(Qpx)+'_f'+float_to_str(fdamp);

	    # output file name
	    os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
	    fileout=path_here+'../../../DELPHI_results/'+machine+'/plot_TMCI_DELPHI_vs_Karliner-Popov_'+machine+'_'+str(round(E/1e9))+'GeV'+model+'_'+str(M)+'b_f'+float_to_str(fdamp)+'_'+str(nmax+1)+'radialmodes_maxazim'+str(lmax)+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];

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
		plot_TMCI(Iscan,lambdax,ax,part=r,leg='DELPHI',patcol='-b',xlab='Intensity [mA]',
			title=machine+", $ Q^' = $ "+str(Qpx)+', $ f= $'+str(fdamp)+', '+str(nmax+1)+' radial modes, '+str(2*lmax+1)+' azimuthal modes'+strfreq[flagdamperimp]);
		
		# Karliner-Popov results
	    	s=read_ncol_file(filename+'_'+r+'.csv',ignored_rows=1);
	    	plot(s[:,0],s[:,1],'Karliner-Popov','xr',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",ax,0,xlab='Intensity [mA]');
		maxs=np.max(s[:,1]);
		
		if ir==1: ax.set_ylim([np.floor(np.min(np.imag(lambdax.flatten()))/0.2)*0.2,np.ceil(maxs/0.2)*0.2]);

		if flagsave: end_figure(fig,ax,save=fileout+'_'+r)
		else: end_figure(fig,ax);


    if not(flagsave): pylab.show();
