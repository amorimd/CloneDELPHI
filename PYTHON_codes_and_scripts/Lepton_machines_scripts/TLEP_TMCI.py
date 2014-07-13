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
from TLEP_conv_2Dplots import TLEPZ_param


if __name__ == "__main__":

    e,m0,c,E0=electron_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model,syncdamp=TLEPZ_param(E0,Zoption=0,flagplot=False,flagsave=False);

    Zd=[];fd=[];
    strfreqflag=['','_freqdepgain'];
    strfreq=['',', freq. dependent gain'];
    strnorm=['','_norm_current_chroma'];
     
    lmax=3;nmax=2;
    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    
    dphasei=0.; # damper phase (0 -> resistive; pi/2 -> reactive)
    dphasestr='';

    os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
    fileoutroot=path_here+'../../../DELPHI_results/plot_TMCI_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV';

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # normalization factor for damper
    dnormfactor=compute_damper_matrix(0,0,nx,M,0.,omega0,Qxfrac,a,b,taub,g,
    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0]; # is actually equal to 1

    # color
    col=['b','r','g','m','k','c','y'];
    
    # scan definition
    Qxscan=np.array([640.1,640.5,640.9]);
    Zoptionscan=np.array([4]);
    Zoptionstr=['Total imp.','RF cav. from BB model','RF cav. from R. Calaga model','Res. wall impedance','Total imp. with BB model'];
    Qpxscan=np.arange(0,1);
    Iscan=np.arange(0.025,5.025,0.025); # intensity scan in mA
    Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch
    dscan=np.arange(0,1); # damper gain scan (inverse of number of damping turns)

    freqshift=np.zeros((len(Qpxscan),len(dscan),len(Nbscan)),dtype=complex);
    Ithres=np.zeros((len(dscan),len(Qpxscan)));
    
    fig=[];ax=[];
    for ir in range(2):
    	fig1,ax1=init_figure();fig.append(fig1);ax.append(ax1);
    
    for iZ,Zoption in enumerate(Zoptionscan):

	# fixed parameters
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model,syncdamp=TLEPZ_param(E0,Zoption=Zoption,flagplot=True,flagsave=True);
	dphase=dphasei;
	print "Impedance:",Zoptionstr[Zoption];
	
	for iQ,Qx in enumerate(Qxscan):

            Qxfrac=Qx-np.floor(Qx);

	
	    for iQp,Qpx in enumerate(Qpxscan):

    		omegaksi=Qpx*omega0/eta;
		headtail_phase=omegaksi*taub;
		print "Qpx=",Qpx,", headtail phase shift=",headtail_phase;

		if flagnorm:
		    # normalization factor for damper at current chromaticity
		    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    		flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
		    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

		matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
			flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);

		matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
			Zx,f,flag_trapz=1,abseps=1e-4);

    		for idamp,damp in enumerate(dscan):

		    #dscan=Iscan*1.e-3*fdamp*2*np.pi*Qs/Ib; # damper gain scan (gain depends on intensity here)
		    #print dscan

		    kmax=(2*lmax+1)*(nmax+1);
		    lambdax=np.zeros((len(Nbscan),kmax),dtype=complex);
		    # computation
		    for iNb,Nb in enumerate(Nbscan):

			coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='electron');

			freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas)

			lambdax[iNb,:]=freqshift/omegas;		

		    # TMCI plots
		    fileout=fileoutroot+model+'_'+str(M)+'b_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_'+str(nmax+1)+'radialmodes_maxazim'+str(lmax)+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];
		    #fileout=fileoutroot+'_model_scan'+'_'+str(M)+'b_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];
		    for ir,r in enumerate(['real','imag']):
    			#fig,ax=init_figure();
			plot_TMCI(Iscan,lambdax,ax[ir],part=r,leg=Zoptionstr[Zoption]+', Q='+str(Qx),patcol=col[(Zoption+len(Zoptionscan)*iQ)%7],xlab='Intensity [mA]',
			    title=machine+", $ Q^' = $ "+str(Qpx)+', $ d= $'+str(damp)+strfreq[flagdamperimp]+dphasestr,ms=6.);

			#if flagsave: end_figure(fig,ax,save=fileout+'_'+r,fontsize=25);
			#else: end_figure(fig,ax,fontsize=25);

	    
    for ir,r in enumerate(['real','imag']):
	if flagsave: end_figure(fig[ir],ax[ir],save=fileout+'_'+r,fontsize=25);
	else: end_figure(fig[ir],ax[ir],fontsize=25);

    if not(flagsave): pylab.show();
