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
from VEPP import VEPP_damper

def TLEPZ_param(E0,Zoption=0,flagplot=False,flagsave=False):

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='TLEPZ';
    E=45.5e9; # injection energy=22 GeV
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    sigmaz=0.19e-2; # RMS bunch length (m) from F. Zimmermann
    taub=4.*sigmaz/(beta*c); # full length in s
    circ=78996.66; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    Qx=640.9;Qxfrac=Qx-np.floor(Qx); # 640 is from F. Zimmermann (?)
    alphap=9e-5; # momentum compaction factor from F. Zimmermann
    M=1; # number of bunches
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    Qs=1.29e3/f0;print "Qs=",Qs; # from F. Zimmermann
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    # synchrotron damping time in seconds (if instability slower than this -> stable)
    syncdamp=0.57; # from formula
    syncdamp=0.6; # from F. Zimmermann

    dphase=0.; # additional damper phase
    nx=0;

    # impedance definition: one broad band model (NB: shunt impedances to be multiplied by beta(location)/(R/Q) )
    beta1=50;betaav=R/Qx # B. Holzer (?)
    print "betaav=",betaav;

    fRF=np.concatenate((10.**np.arange(-1,7),np.arange(5.e7,1.0005e11,5.e7),10.**np.arange(11.1,13.1,0.1),
    	10.**np.arange(14,16)));
    R1=1.5e3*600/1.07*beta1/betaav;f1=5e9;Q1=1; # from Rama Calaga: ~3 kOhm/cavity (700 MHz)
    ZRF=resonator_impedance(R1,f1,Q1,fRF);
    
    # real impedance of RF cavities (file from Rama)
    filename=path_here+'RF_700MHz_from_Rama_Calaga/fft.bnl2EqualSmooth.dLW';
    sr=read_ncol_file(filename,ignored_rows=14);
    si=read_ncol_file(filename,ignored_rows=12500);
    fRF_Rama=sr[:,0]*1e9;ZRF_Rama=np.zeros((len(fRF_Rama),2));
    ZRF_Rama[:,0]=sr[:,1]*1e3*600/1.07;
    ZRF_Rama[:,1]=si[:,1]*1e3*600/1.07;
    # extend with zeros to avoid extrapolation
    fRF_Rama=np.concatenate((np.array([0.1,fRF_Rama[0]-1]),fRF_Rama,np.array([fRF_Rama[-1]+1,1e15])))
    ZRF_Rama=np.concatenate((np.array([[0.,ZRF_Rama[0,1]],ZRF_Rama[0,:]]),ZRF_Rama,np.array([[0.,0.],[0.,0.]])))
     
    # resistive-wall impedance
    radius=20; # beam pipe radius
    filename='../../ImpedanceWake2D/TLEP/ZxdipWTLEP_1layers'+str(radius)+'.00mm_TLEP_Al_beampipe.dat';
    s=read_ncol_file(filename,ignored_rows=1);
    fRes=s[:,0];ZRes=np.zeros((len(fRes),2));
    ZRes[:,0]=s[:,1];
    ZRes[:,1]=s[:,2];
    
    # total impedance
    # frequencies
    f=sort_and_delete_duplicates(np.concatenate((fRF_Rama,fRes)),tolerance=0.1);
    Zx=np.zeros((len(f),2));
    Zx[:,0]=np.interp(f,fRF_Rama,ZRF_Rama[:,0])+np.interp(f,fRes,ZRes[:,0]);
    Zx[:,1]=np.interp(f,fRF_Rama,ZRF_Rama[:,1])+np.interp(f,fRes,ZRes[:,1]);
    
    f2=sort_and_delete_duplicates(np.concatenate((fRF,fRes)),tolerance=0.1);
    Zx2=np.zeros((len(f2),2));
    Zx2[:,0]=np.interp(f2,fRF,ZRF[:,0])+np.interp(f2,fRes,ZRes[:,0]);
    Zx2[:,1]=np.interp(f2,fRF,ZRF[:,1])+np.interp(f2,fRes,ZRes[:,1]);
    
    # plots
    if (flagplot):
	fig,ax=init_figure();
	plot(fRF,ZRF[:,0],'Re[Z], RF cav., broadband model','--b',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(fRF,np.abs(ZRF[:,1]),'|Im[Z]|, RF cav., broadband model','--r',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(fRF_Rama,ZRF_Rama[:,0],'Re[Z], RF cav., R. Calaga model','-b',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(fRF_Rama,ZRF_Rama[:,1],'Im[Z], RF cav., R. Calaga model','-r',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(fRes,ZRes[:,0],'Re[Z], res. wall (Al), radius='+str(radius)+' mm','-g',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(fRes,np.abs(ZRes[:,1]),'|Im[Z]|, res. wall (Al), radius='+str(radius)+' mm','-m',r" $Z [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')

	plot(f,Zx[:,0],'Re[Z], total imp.','-c',r" $Z  [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	plot(f,np.abs(Zx[:,1]),'|Im[Z]|, total imp.','-k',r" $Z  [\Omega/m] $ ",ax,3,xlab='Frequency [Hz]')
	ax.set_xlim([1e3,1e13]);ax.set_ylim([10,1e9]);

	if not(flagsave): end_figure(fig,ax,fontsize=20)
	else: end_figure(fig,ax,save=path_here+'../../../DELPHI_results/TLEPZ/TLEPZ_impedance_model_comp',fontsize=20)
    
    if Zoption==0:
    	# total realistic impedance
    	return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,'_total_imp',syncdamp;

    elif Zoption==1:
    	# RF broadband only
    	return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,ZRF,fRF,'_RF_BB',syncdamp;

    elif Zoption==2:
    	# RF from Rama Calaga only
    	return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,ZRF_Rama,fRF_Rama,'_RF_Rama',syncdamp;

    elif Zoption==3:
    	# RF resistive wall only
    	return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,ZRes,fRes,'_res_wall_Al_'+str(radius)+'mm',syncdamp;

    elif Zoption==4:
    	# total impedance with RF BB
    	return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx2,f2,'_total_imp_with_BB',syncdamp;


if __name__ == "__main__":

    e,m0,c,E0=electron_param();
    Zoption=4;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model,syncdamp=TLEPZ_param(E0,Zoption=Zoption,flagplot=True,flagsave=False);

    Qx=640.9;

    Zd=[];fd=[];
    strfreqflag=['','_freqdepgain'];
    strfreq=['',', freq. dependent gain'];
    strnorm=['','_norm_current_chroma'];
     
    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    
    # synchrotron damping time in seconds (if instability slower than this -> stable)
    dphase=-np.pi/2.; # damper phase (0 -> resistive; -pi/2 -> reactive)
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of eigenvalues used for the TMCI plots

    if dphase==0: dphasestr=', resistive damper';
    elif dphase==-np.pi/2: dphasestr=', reactive damper';

    os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
    fileoutroot=path_here+'../../../DELPHI_results/'+machine+'plot_TMCI_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+model+'_'+str(M)+'b_Q'+float_to_str(Qx-np.floor(Qx));
    fileoutroot2D=path_here+'../../../DELPHI_results/'+machine+'plot2D_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+model+'_'+str(M)+'b_Q'+float_to_str(Qx-np.floor(Qx));
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # normalization factor for damper
    dnormfactor=compute_damper_matrix(0,0,nx,M,0.,omega0,Qxfrac,a,b,taub,g,
    	flagdamperimp=flagdamperimp,d=Zd,freqd=fd,abseps=1e-4);
    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0]; # is actually equal to 1

    # scan definition
    Qpxscan=np.arange(-200,210,10);
    Iscan=np.arange(0.025,5.025,0.025); # intensity scan in mA
    Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch
    dscan=np.arange(0,0.1+0.005,0.005); # damper gain scan (inverse of number of damping turns)

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
			freqd=fd,kmax=kmax,crit=5.e-2,abseps=1.,lmaxold=lmax,
			nmaxold=nmax,matdamperold=matdamper,matZold=matZ);
			
	    	#print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;
		freqshift[iQp,idamp,iNb]=freqshifttmp[0];
		lambdax[iNb,:]=freqshifttmp[:kmaxplot]/omegas;		
		
	    # find intensity threshold
	    Ithres[idamp,iQp]=find_intensity_threshold(Iscan,freqshift[iQp,idamp,:],thresgrowth=1./syncdamp);
	    print "Qp=",Qpx,", damp=",damp,',lmax=',lmax,', nmax=',nmax,"Ithres=",Ithres[idamp,iQp];

	    # TMCI plots
	    if (iQp%5==0)and(idamp%10==0):
		fileout=fileoutroot+'_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qpx)+strfreqflag[flagdamperimp]+strnorm[flagnorm];
	        for ir,r in enumerate(['real','imag']):
    		    fig,ax=init_figure();
		    plot_TMCI(Iscan,lambdax,ax,part=r,leg='DELPHI',patcol='.b',xlab='Intensity [mA]',
			title=machine+", $ Q^' = $ "+str(Qpx)+', $ d= $'+str(damp)+strfreq[flagdamperimp]+dphasestr,ms=6.);

		    if flagsave: end_figure(fig,ax,save=fileout+'_'+r,fontsize=25);
		    else: end_figure(fig,ax,fontsize=25);
 	    	    
	    
    # 2D plot

    # output file name (root)
    fileout=fileoutroot2D+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+strfreqflag[flagdamperimp]+strnorm[flagnorm];

    fig,ax=init_figure();
    
    plot2D(Ithres,Qpxscan[0],Qpxscan[-1],dscan[0],dscan[-1]," $ Q^' $ ","Damping rate (1/nb turns)",machine+dphasestr+", single-bunch instability threshold vs. $ Q^' $ and damping rate",ax,colorlabel='I [mA]');
   
    if flagsave: end_figure(fig,ax,save=fileout,fontsize=20);
    else: end_figure(fig,ax,fontsize=20);




    if not(flagsave): pylab.show();
