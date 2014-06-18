#!/usr/bin/python2.6

import sys
import numpy as np
import pylab,os,re
sys.path.append("../PYTHON/")
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
    
    # synchrotron damping time in seconds (if instability slower than this -> stable)
    syncdamp=0.57; # from formula
    syncdamp=0.3; # from J. Jowett (1987) with wigglers
    sigmaz=0.0094; # from J. Jowett (1987) with wigglers
    dphase=-np.pi/2.; # damper phase (0 -> resistive; pi/2 -> reactive)
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of eigenvalues used for the plot

    if dphase==0: dphasestr=', resistive damper';
    elif dphase==-np.pi/2: dphasestr=', reactive damper';

    Ib=4*np.pi*E*Qx*Qs*sigmaz/(R**2*R1) # some normalization factor

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # normalization factor for damper
    dnormfactor=compute_damper_matrix(0,0,nx,M,0.,omega0,Qxfrac,a,b,taub,g,
    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0]; # is actually equal to 1

    # scan definition
    Qpxscan=np.arange(-20,21);
    Iscan=np.arange(0.025,5.025,0.025); # intensity scan in mA
    Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch
    dscan=np.arange(0,0.1,0.005); # damper gain scan (inverse of number of damping turns)

    freqshift=np.zeros((len(Qpxscan),len(dscan),len(Nbscan)),dtype=complex);
    Ithres=np.zeros((len(dscan),len(Qpxscan)));
    
    for iQp,Qpx in enumerate(Qpxscan):
    
    	omegaksi=Qpx*omega0/eta;
	headtail_phase=omegaksi*taub;
	print "Qpx=",Qpx,", headtail phase shift=",headtail_phase;
	
	if flagnorm:
	    # normalization factor for damper at current chromaticity
	    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
	    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];
	
    	lmax=-1;nmax=-1;matZ=None;matdamper=None;
    
    	for idamp,damp in enumerate(dscan):
	
	    #dscan=Iscan*1.e-3*fdamp*2*np.pi*Qs/Ib; # damper gain scan (gain depends on intensity here)
	    #print dscan
	    
	    lambdax=np.zeros((len(Nbscan),kmaxplot),dtype=complex);
	    # computation
	    for iNb,Nb in enumerate(Nbscan):
	    
		coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='electron');
		
		freqshifttmp,v,lmax,nmax,matdamper,matZ=eigenmodesDELPHI_converged(nx,
			M,omegaksi,omega0,Qxfrac,a,b,taub,g,Zx,f,coefdamper,
			coefZ,omegas,flag_trapz=1,flagdamperimp=flagdamperimp,d=Zd,
			freqd=fd,kmax=kmax,crit=5.e-2,abseps=1.e-3,lmaxold=lmax,
			nmaxold=nmax,matdamperold=matdamper,matZold=matZ);
			
	    	#print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;
		freqshift[iQp,idamp,iNb]=freqshifttmp[0];
		lambdax[iNb,:]=freqshifttmp[:kmaxplot]/omegas;		
		
	    # find intensity threshold
	    Ithres[idamp,iQp]=find_intensity_threshold(Iscan,freqshift[iQp,idamp,:],thresgrowth=1./syncdamp);
	    print "Qp=",Qpx,", damp=",damp,',lmax=',lmax,', nmax=',nmax,"Ithres=",Ithres[idamp,iQp];

	    # TMCI plots
	    if (iQp%10==0)and(idamp%2==0):
		fileout='../DELPHI_results/plot_TMCI_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+model+'_'+str(M)+'b_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];
	        for ir,r in enumerate(['real','imag']):
    		    fig,ax=init_figure();
		    plot_TMCI(Iscan,lambdax,ax,part=r,leg='DELPHI',patcol='b',xlab='Intensity [mA]',
			title=machine+", $ Q^' = $ "+str(Qpx)+', $ d= $'+str(damp)+strfreq[flagdamperimp]+dphasestr);

		    if flagsave: end_figure(fig,ax,save=fileout+'_'+r,fontsize=25);
		    else: end_figure(fig,ax,fontsize=25);
 	    	    
	    
    # 2D plot

    # output file name (root)
    fileout='../DELPHI_results/plot2D_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+model+'_'+str(M)+'b_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+strfreqflag[flagdamperimp]+strnorm[flagnorm];

    fig,ax=init_figure();
    
    plot2D(Ithres,Qpxscan[0],Qpxscan[-1],dscan[0],dscan[-1]," $ Q^' $ ","Damping rate (1/nb turns)",machine+dphasestr+", single-bunch instability threshold vs.$ Q^' $ and damping rate",ax,colorlabel='I [mA]');
   
    if flagsave: end_figure(fig,ax,save=fileout,fontsize=20);
    else: end_figure(fig,ax,fontsize=20);




    if not(flagsave): pylab.show();
