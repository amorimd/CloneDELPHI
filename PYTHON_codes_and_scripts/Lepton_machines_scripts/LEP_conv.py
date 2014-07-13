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
from io_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *
from LEP import *
from VEPP import VEPP_damper

if __name__ == "__main__":

    e,m0,c,E0=electron_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model=LEP_param(E0);

    # VEPP-4 damper model
    Zd,fd=VEPP_damper(R,Qx,f0);
    strfreqflag=['','_freqdepgain'];
    strfreq=['',', freq. dependent gain'];
    strnorm=['','_norm_current_chroma'];
     
    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    
    kmax=5; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=15; # number of eigenvalues used for the plot

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
	
	if (Qpx==0):
	    fdampscan=np.array([0.,2.]);
	    Iscan=np.arange(0.025,1.525,0.025); # intensity scan in mA
	elif (Qpx==-22):
	    fdampscan=np.array([2.]);
	    Iscan=np.arange(0.025,5.025,0.025); # intensity scan in mA
	
     	Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch
    	lmax=-1;nmax=-1;matZ=None;matdamper=None;
    
    	for idamp,fdamp in enumerate(fdampscan):
	
	    print "Qpx=",Qpx,", f=",fdamp;
	    
	    # Karliner-Popov results
	    filename=path_here+'Karliner_Popov_figures/Karliner_Popov_Qp'+float_to_str(Qpx)+'_f'+float_to_str(fdamp);

	    # output file name
	    os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
	    fileout=path_here+'../../../DELPHI_results/'+machine+'/plot_TMCI_DELPHI_vs_Karliner-Popov_'+machine+'_'+str(round(E/1e9))+'GeV'+model+'_'+str(M)+'b_f'+float_to_str(fdamp)+'_converged'+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];

	    dscan=Iscan*1.e-3*fdamp*2*np.pi*Qs/Ib; # damper gain scan (gain depends on intensity here)
	    #print dscan
	    lambdax=np.zeros((len(Nbscan),kmaxplot),dtype=complex);
	    
	    # computation
	    for iNb,Nb in enumerate(Nbscan):
	    
	    	#print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;
		
		coefdamper,coefZ=computes_coef(f0,dscan[iNb],b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='electron');
		
		freqshift,v,lmax,nmax,matdamper,matZ=eigenmodesDELPHI_converged(nx,
			M,omegaksi,omega0,Qxfrac,a,b,taub,g,Zx,f,coefdamper,
			coefZ,omegas,flag_trapz=1,flagdamperimp=flagdamperimp,d=Zd,
			freqd=fd,kmax=kmax,crit=5.e-2,abseps=1.e-3,lmaxold=lmax,
			nmaxold=nmax,matdamperold=matdamper,matZold=matZ);
		
		lambdax[iNb,:]=freqshift[:kmaxplot]/omegas;		
    	    
	    
	    # TMCI plots
	    strpart=['Re','Im'];
	    for ir,r in enumerate(['real','imag']):
		
		fig,ax=init_figure();
		plot_TMCI(Iscan,lambdax,ax,part=r,leg='DELPHI',patcol='-b',xlab='Intensity [mA]',
			title=machine+", $ Q^' = $ "+str(Qpx)+', $ f= $'+str(fdamp)+', converged calculation'+strfreq[flagdamperimp]);
		
		# Karliner-Popov results
		maxs=0.;
		for k in range(1,11):
	    	    s=read_ncol_file(filename+'_'+r+'.csv',ignored_rows=1);
		    if (k==1): leg='Karliner-Popov';
		    else: leg='';
	    	    plot(s[:,0],s[:,1],leg,'xr',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",ax,0,xlab='Intensity [mA]');
		    maxs=max(maxs,np.max(s[:,1]));
		
		if ir==1: ax.set_ylim([np.floor(np.min(np.imag(lambdax.flatten()))/0.2)*0.2,np.ceil(maxs/0.2)*0.2]);

		if flagsave: end_figure(fig,ax,save=fileout+'_'+r)
		else: end_figure(fig,ax);



    if not(flagsave): pylab.show();
