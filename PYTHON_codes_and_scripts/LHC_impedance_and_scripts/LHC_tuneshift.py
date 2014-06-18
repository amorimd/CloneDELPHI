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

from string import *
import time
import numpy as np
from copy import deepcopy
import pylab,os,re
sys.path.append("../PYTHON/")
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from tables_lib import select_in_table
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_conv import LHC_param
from LHC_imp import *
from HLLHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    E=4e12;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);
    avbetax=R/Qx;avbetay=R/Qy;
    sigmaz=0.0847; # to match Daria's computations
    beta=np.sqrt(1.-1./(gamma**2))
    taub=4.*sigmaz/(beta*c);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/LHC_tuneshift';
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=True; # to plot impedances
    nevery=1; # downsampling of the impedance (take less points than in the full model)

    wake_calc=True; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition (including measured parameters)
    beam='1';
    Qpscan=np.arange(-10.,11);
    Nbscan=np.array([100.])*1e9; # intensity scan
    
    damp=0;M=1;scenario='_2012_v2';
    squeeze='0p6m_3m_0p6m_3m';
   
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((2,len(Qpscan),1,len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((2,len(Qpscan),1,len(Nbscan),1,1),dtype=complex);

    param_filename_coll='../Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
    settings_filename_coll='../Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';

    # compute imp. model
    imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
	    settings_filename_coll,dire="../LHC_elements/",commentcoll=scenario,
	    direcoll='Coll'+scenario+'/',
	    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeeze,
	    wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario)
	
	
    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# dump into a file
	#filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	#pick.dump(imp_mod,filemodel);
	#filemodel.close();

	# write Ascii files with each component
	#write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	#	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
	#	dire=root_result+'/')

	if (wake_calc):
	    # write Ascii files with each component
	    #write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    #    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
	    #    dire=root_result+'/')

	    # write HEADTAIL input file
	    # first do the longitudinal matching
	    delta_lin,eps0_lin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=True);
	    print "Longitudinal linear matching: delta_p/p=",delta_lin,", eps0=",eps0_lin,", Qs=",Qs;
	    delta_nlin,eps0_nlin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=False);
	    print "Longitudinal non-linear matching: delta_p/p=",delta_nlin,", eps0=",eps0_nlin,", Qs=",Qs;
	    # check that Qs is correct
	    if ((Qs-Qs_from_RF_param(V,h,gamma,eta))>1e-8): print " Pb in Qs !!";sys.exit();

	    # write .cfg files (linear and non-linear matching)
	    cfgnamelin=root_result+'/'+machine+"_"+Estr+scenario+'_lin.cfg';
	    write_CFG_HEADTAIL(cfgnamelin,Nb=Nbscan[0],betax=R/Qx,betay=R/Qy,
	    	sigmaz=sigmaz,delta=delta_lin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		gamma=gamma,nturns=200000,Qx=Qx,Qy=Qy,isyn=1,start_turn=199000,end_turn=199100,VRF=V,
		dratex=damp,dratey=damp);

	    cfgnamenlin=root_result+'/'+machine+"_"+Estr+scenario+'_nlin.cfg';
	    write_CFG_HEADTAIL(cfgnamenlin,Nb=Nbscan[0],betax=R/Qx,betay=R/Qy,
	    	sigmaz=sigmaz,delta=delta_nlin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		gamma=gamma,nturns=200000,Qx=Qx,Qy=Qy,isyn=4,start_turn=199000,end_turn=199100,VRF=V,
		dratex=damp,dratey=damp);
	    
	    # wake for HEADTAIL
	    write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		

    
	# DELPHI computation
	for iplane,plane in enumerate(['x','y']):


	    # total impedance

	    # select Zxdip+Zxquad or Zydip+Zyquad
	    for iw in imp_mod:
		if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Zd=deepcopy(iw.func[::nevery,:]);freqd=deepcopy(iw.var[::nevery]);
	    for iw in imp_mod:
		if test_impedance_wake_comp(iw,0,0,1-iplane,iplane,plane): Zq=deepcopy(iw.func[::nevery,:]);freqq=deepcopy(iw.var[::nevery]);
	    # sum the two
	    freq=sort_and_delete_duplicates(np.concatenate((freqd,freqq)));
	    Z=np.zeros((len(freq),2),dtype=float);
	    for icol in range(2): Z[:,icol]=np.interp(freq,freqd,Zd[:,icol],right=0.)+np.interp(freq,freqq,Zq[:,icol],right=0.);
		
	    for iQp,Qp in enumerate(Qpscan):

		tuneshiftnx=np.zeros((1,1,1,kmaxplot),dtype=complex);
		
	    	tuneshiftQp[iplane,iQp,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[iplane,iQp,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    [0],[damp],Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			    flag_trapz=1,flagdamperimp=0,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5e-2,abseps=1e-4,flagm0=True,
			    lxplusbatch=lxplusbatchDEL,comment=scenario+'_'+plane+'_Qp'+float_to_str(Qp),
			    queue='8nh',dire=root_result+'/');
	    
	    
    # now the plots (outside loop on scenarios)
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	for iplane,plane in enumerate(['x','y']):

	    # plots vs Q'
	    for Nb in Nbscan:

		# initialize plots vs Q'
		figQp=[];axQp=[];figQpm0=[];axQpm0=[];
		for ir in range(2):
		    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
		    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

		# output files name for plots vs Q'
		fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+'_'+plane;
		fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+'_'+plane;

		strpart=['Re','Im'];
		for ir,r in enumerate(['real','imag']):

		    # output files name for data vs Qp
		    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+'_'+plane;
		    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+'_'+plane;

		    #if flagcompute:
		    ts=getattr(tuneshiftQp[iplane,:,0,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
		    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
		    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")

		    tsm0=getattr(tuneshiftm0Qp[iplane,:,0,pylab.mlab.find(Nbscan==Nb),0,0],r);
		    data=np.hstack((Qpscan.reshape((-1,1)),tsm0.reshape((-1,1))));
		    write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
		    #else:
		    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
		    #    Qpscan=s[:,0];ts=s[:,1];

		    sgn=1;sgnstr='';
		    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
		    ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
		    if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate

		    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))],'DELPHI','-b',ylab,axQp[ir],0,xlab=" $ Q^' $ ");
		    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))],'DELPHI','-b',ylab+" (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");


		for ir,r in enumerate(['real','imag']):
		    end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
		    end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))


    if not(flagsave): pylab.show();
