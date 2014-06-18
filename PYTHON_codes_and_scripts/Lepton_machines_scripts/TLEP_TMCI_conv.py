#!/usr/bin/python2.6

import sys
if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print lxplusbatchImp,lxplusbatchDEL;   

import commands
out=commands.getoutput("hostname")
if out.startswith('lxplus'):
    sys.path.insert(1,'/afs/cern.ch/user/n/nmounet/private/soft/Pymodules/numpy-install/lib64/python2.6/site-packages');
    sys.path.insert(1,'/afs/cern.ch/user/n/nmounet/private/soft/Pymodules/scipy-install/lib64/python2.6/site-packages');
    sys.path.insert(1,'/afs/cern.ch/user/n/nmounet/private/soft/Pymodules/matplotlib-install/lib64/python2.6/site-packages');

import numpy as np
import pylab,os,re
from copy import deepcopy
sys.path.append("../PYTHON/")
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from string_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *
#from VEPP import VEPP_damper
#from TLEP_conv_2Dplots import TLEPZ_param
from TLEP_imp import *


if __name__ == "__main__":

    e,m0,c,E0=electron_param();
    E=45e9;option='Z';
    optionscan=['Z','Z4C','W','H','t','tB'];
    Escan=np.array([45,45,80,120,175,175])*1e9;
    optionscan=['Z'];
    Escan=np.array([45])*1e9;
    hgap=0.03; # vertical semi-axis of vacuum pipe

    for ioption,option in enumerate(optionscan):

	E=Escan[ioption];
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,syncdamp=TLEP_param(E=E,option=option);
	beta=np.sqrt(1.-1./(gamma**2))
	# older parameters for Z option
	#R=78996.66/2/np.pi;sigmaz=0.19e-2;taub=4.*sigmaz/(beta*c);
	#f0=c*beta/78996.66;omega0=2.*np.pi*f0;
	#Qs=1.29e3/f0;omegas=Qs*omega0;alphap=9e-5;eta=alphap-1./(gamma*gamma);
	#Qx=640+Qxfrac;Qy=640+Qyfrac;

	# fixed parameters
	#Zoption=4;
	#machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qs,eta,M,f0,omega0,omegas,dphase,nx,R1,Zx,f,model,syncdamp=TLEPZ_param(E0,Zoption=Zoption,flagplot=True,flagsave=False);

	Zd=[];fd=[];
	strfreqflag=['','_freqdepgain'];
	strfreq=['',', freq. dependent gain'];
	strnorm=['','_norm_current_chroma'];

	flagsamefig=1; # 1 if all TMCI plots on the same figure
	flagsave=0; # 1 to save figure instead of plotting on screen
	flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
	flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)

	# synchrotron damping time in seconds (if instability slower than this -> stable)
	dphase=0.; # damper phase (0 -> resistive; pi/2 -> reactive)
	kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
	kmaxplot=50; # number of eigenvalues used for the plot

	if dphase==0: dphasestr=', resistive damper';
	elif dphase==-np.pi/2: dphasestr=', reactive damper';

	root_result='../DELPHI_results/'+machine;#+'_old';
	#fileoutroot='../DELPHI_results/plot_TMCI_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+model+'_'+str(M)+'b';

	# longitudinal distribution initialization
	g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

	# color
	col=['b','r','g','m','k','c','y'];mark=['.','x'];

	# scan definition
	#Qxscan=np.array([640.1,640.5,640.9]);
	Zoptionscan=['all','RW','RWnoabs','BB','RF','abs'];
	Zoptionscanleg=['Total impedance','Total resistive-wall imp.','Resistive-wall from vacuum pipe','Tapers (for photon absorbers) + RF','RF only','Total absorbers imp.'];
	Zoptionscan=['all','RWnoabs','BB'];
	Zoptionscanleg=['Total impedance','Resistive-wall from vacuum pipe','Broad-band model (tapers + RF)'];
	Zoptionscan=['all','RWnoabs'];
	Zoptionscanleg=['Total impedance','Resistive-wall from vacuum pipe'];
	ldipscan=np.array([5.5]);
	#Zoptionscan=['all'];
	#Zoptionscanleg=['Total impedance'];
	#ldipscan=np.array([1.,2.,3.5,5.5,8.,11.]);
	Qpscan=np.array([0.]);
	#Iscan=np.arange(0.01,2.005,0.005); # intensity scan in mA
	#Nbscan=Iscan*1.e-3/(e*f0); # number of particles per bunch
	Nbscan=np.arange(0.01,8.01,0.01)*1e11; # number of particles per bunch
	Iscan=Nbscan*1e3*e*f0; # intensity scan in mA
	#dscan=np.arange(0,0.+0.005,0.005); # damper gain scan (inverse of number of damping turns)
	dampscan=np.array([0.]); # damper gain scan (inverse of number of damping turns)
	Mscan=np.array([1]);

	tuneshiftQp=np.zeros((len(Zoptionscan),len(ldipscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
	Ithres=np.zeros((len(Zoptionscan),len(ldipscan),2,len(Mscan),len(dampscan),len(Qpscan)));
	imp_mod_list=[];wake_mod_list=[];

	if flagsamefig:
    	    fig=[];ax=[];
	    for ir in range(2):
    		fig1,ax1=init_figure();fig.append(fig1);ax.append(ax1);

	ldipscanleg=[];
	for iZ,Zoption in enumerate(Zoptionscan):
	    
    	    for ildip,ldip in enumerate(ldipscan):

		ldipscanleg.append('length dip.='+str(ldip)+'m');

		imp_mod,wake_mod=TLEP_imp(gamma,2*np.pi*R,R/Qx,R/Qy,length_dip=ldip,
    		    length_absorber=0.5,bending_radius=11000.,b=hgap,
		    dire='TLEP/',wake_calc=True,lxplusbatch=lxplusbatchImp,option=Zoption)

		# number of components for HEADTAIL wake file
		if (Zoption=='BB')or(Zoption=='RF'):
		    # add quadrupolar impedance filled with zeros
		    imp_mod.append(impedance_wake(a=0,b=0,c=1,d=0,plane='x',
	    		var=imp_mod[0].var,func=np.zeros((len(imp_mod[0].var),2))));
		    imp_mod.append(impedance_wake(a=0,b=0,c=0,d=1,plane='y',
	    		var=imp_mod[0].var,func=np.zeros((len(imp_mod[0].var),2))));
		    wake_mod.append(impedance_wake(a=0,b=0,c=1,d=0,plane='x',
	    		var=wake_mod[0].var,func=np.zeros((len(wake_mod[0].var),2))));
		    wake_mod.append(impedance_wake(a=0,b=0,c=0,d=1,plane='y',
	    		var=wake_mod[0].var,func=np.zeros((len(wake_mod[0].var),2))));
		

		if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('ret')):

		    imp_mod_list.append(imp_mod);
		    wake_mod_list.append(wake_mod);

		    # output in a HEADTAIL kind of wake file
		    wakefilename=root_result+'/wakeforHDTL_'+machine+'_'+Estr+'_'+Zoption+'_ldip'+float_to_str(ldip)+'m_wake.dat';
		    write_wake_HEADTAIL(wake_mod,wakefilename,beta=beta,ncomp=4);

		    for iplane,plane in enumerate(['x','y']):
			# select Zxdip or Zydip
			for iw in imp_mod:
			    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);

			for iM,M in enumerate(Mscan):

			    # normalization factor for damper
			    dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
    				flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
			    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

			    flag_trapz=0; # by default no trapz method

        		    if (M==1): nxscan=np.array([0]);flag_trapz=1;
			    #elif (M==1782): nxscan=np.array([0, 1, 300, 600, 880, 890, 891, 892, 900, 910, 950, 1000, 1200, 1500, 1780, 1781])
			    #elif (M==3564): nxscan=np.array([0, 1, 300, 600, 900, 1200, 1500, 1770, 1780, 1781, 1782, 1785, 1790, 1800, 1900, 2000, 2300, 2600, 2900, 3200, 3500, 3560, 3561, 3562, 3563])
			    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,100),np.arange(M/2-10,M/2+11),
	    			    np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

			    # to put loop on Qpscan outside the lxplus job
			    for iQp,Qp in enumerate(Qpscan):
				tuneshiftnx=np.zeros((1,len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
				#tuneshiftQp[iZ,ildip,iplane,iM,iQp,:,:len(Nbscan),:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
		    		#	nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
		    		#	a,b,taub,g,Z,freq,particle='electron',flagnorm=flagnorm,
		    		#	flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
		    		#	kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1.,
		    		#	lxplusbatch=lxplusbatchDEL,comment=machine+Zoption+'_ldip'+float_to_str(ldip)+'m_b'+float_to_str(hgap)+'mm_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
		    		#	queue='8nh',dire=root_result+'/');

				if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('ret')):

				    for idamp,damp in enumerate(dampscan):

					lambdax=tuneshiftQp[iZ,ildip,iplane,iM,iQp,idamp,:,0,0,:]/Qs;		

					# find intensity threshold
					Ithres[iZ,ildip,iplane,iM,idamp,iQp]=find_intensity_threshold(Iscan,tuneshiftQp[iZ,ildip,iplane,iM,iQp,idamp,:len(Nbscan),0,0,0]*omega0,thresgrowth=1./syncdamp);
					print machine,Zoption,", ldip=",ldip,plane,",M=",M,", Qp=",Qp,", damp=",damp,", Ithres=",Ithres[iZ,ildip,iplane,iM,idamp,iQp];

		    			# TMCI plots
					if not(flagsamefig): fileout=root_result+'/plot_TMCI_DELPHI_'+machine+'_'+float_to_str(E/1e9)+'GeV'+'_'+str(M)+'b_'+Zoption+'_ldip'+float_to_str(ldip)+'m_b'+float_to_str(hgap)+'mm_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qp)+'_'+plane+strfreqflag[flagdamperimp]+strnorm[flagnorm];
	        			else: fileout=root_result+'/plot_TMCI_DELPHI_'+machine+'_b'+float_to_str(hgap)+'mm_'+float_to_str(E/1e9)+'GeV'+'_'+str(M)+'b'+'_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qp)+strfreqflag[flagdamperimp]+strnorm[flagnorm];

					for ir,r in enumerate(['real','imag']):

					    if ir==0: ylim=[-5,5];
					    else: ylim=[-5,5];
					    
					    if not(flagsamefig):
						fig,ax=init_figure();patcol='.b';
						plot_TMCI(Nbscan,lambdax,ax,part=r,leg='',patcol='.b',xlab='Intensity [p+/b]',
						    title=machine+", $ Q^' = $ "+str(Qp)+', $ d= $'+str(damp)+strfreq[flagdamperimp]+dphasestr,ms=6.,ylim=ylim);
					    else: plot_TMCI(Nbscan,lambdax,ax[ir],part=r,leg=Zoptionscanleg[iZ]+', length dip.='+str(ldip)+'m, '+plane,patcol=mark[iplane]+col[iZ+ildip],xlab='Intensity [p+/b]',
						title=machine+", $ Q^' = $ "+str(Qp)+', $ d= $'+str(damp)+strfreq[flagdamperimp]+dphasestr,ms=6.,ylim=ylim);

					    if not(flagsamefig):
						end_figure(fig,ax,save=flagsave*(fileout+'_'+r),fontsize=25);


	if (flagsamefig)and((lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('ret'))):

	    for ir,r in enumerate(['real','imag']):
		end_figure(fig[ir],ax[ir],save=flagsave*(fileout+'_'+r),fontsize=25);


	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('ret')):

	    if len(ldipscan)>1: ldipstr='';legimp=ldipscanleg;
	    else: ldipstr='ldip'+float_to_str(ldipscan[0])+'m_';legimp=Zoptionscanleg;
	    
	    # compare the impedance models
	    plot_compare_imp_model(imp_mod_list,legimp,listcomp=['Zxdip','Zydip','Zxquad','Zyquad'],
		saveimp=root_result+'/plot_imp_'+machine+'_b'+float_to_str(hgap)+'mm_'+ldipstr+'all_options',
		saveratio=root_result+'/plot_imp_ratio_'+machine+'_b'+float_to_str(hgap)+'mm_'+ldipstr+'all_options',
		xlim=[1e3,1e12],ylim=[1e3,1e11],bounds=[20e6,3./taub],beta=beta)

            plot_compare_imp_model(wake_mod_list,legimp,listcomp=['Wxdip','Wydip','Wxquad','Wyquad'],
                saveimp=root_result+'/plot_wake_'+machine+'_b'+float_to_str(hgap)+'mm_'+ldipstr+'all_options',
                saveratio=root_result+'/plot_wake_ratio_'+machine+'_b'+float_to_str(hgap)+'mm_'+ldipstr+'all_options',
                xlim=[1e-6,100],ylim=[1e14,1e19],bounds=[20e6,3./taub],beta=beta)


	if (len(ldipscan)>1)and((lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('ret'))):

	    for iZ,Zoption in enumerate(Zoptionscan):

		for iplane,plane in enumerate(['x','y']):

	    	    for iM,M in enumerate(Mscan):

			for idamp,damp in enumerate(dampscan):

		    	    for iQp,Qp in enumerate(Qpscan):

				# plot intensity threshold vs length dipoles
				Nbthres=Ithres[iZ,:,iplane,iM,idamp,iQp]*1e-3/(e*f0);

		    		fileoutldip=root_result+'/plot_thres_vs_ldip_DELPHI_'+machine+'_b'+float_to_str(hgap)+'mm_'+float_to_str(E/1e9)+'GeV'+'_'+str(M)+'b_'+Zoption+'_d'+float_to_str(damp)+'_dphase'+float_to_str(round(dphase*100)/100)+'_converged'+'_Qp'+float_to_str(Qp)+'_'+plane+strfreqflag[flagdamperimp]+strnorm[flagnorm];

				figd,axd=init_figure();
				plot(ldipscan,Nbthres,Zoption,'b','TMCI threshold [p+/b]',axd,0,xlab='Dipole length [m]');
				end_figure(figd,axd,flagsave*fileoutldip);

    if not(flagsave): pylab.show();
