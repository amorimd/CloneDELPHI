#!/usr/bin/python

import sys
if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print lxplusbatchImp,lxplusbatchDEL;   

from string import *
import numpy as np
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure,cmap
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param
from LHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
       
    kmax=5; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=50; # number of converged eigenvalues (kmax most unstable ones are converged)
    root_result=path_here+'../../../DELPHI_results/'+machine+'/somedeviceIR3';
    os.system("mkdir -p "+root_result);
    suffix='';#suffix='_only_TCSG_IR7' # suffix for output files 
    
    model='_somedeviceIR3';
    name='somedeviceIR3_ss304L';material='ss304L';radius=0.04;
    angle=0;length=0.2;betaratio=3.; # last is ratio w.r.t. average beta functions

    # scan definitions
    coatingscan=1e-6*np.array([0,1,2,3,5,8,10,12,15,20,35,50,100]);
    Qpscan=np.arange(-10,21,1);
    dampscan=np.array([0,0.02]); # damper gain scan
    #Nbscan=np.arange(1.e10,5.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.5e11]); # intensity scan
    Mscan=np.array([1,1782]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
        
    # compute the total impedance
    scenarioscan=['_2012','_relaxed','_mm_kept','_nominal','_sigma_kept'];
    iscenario=1;scenario=scenarioscan[iscenario];
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_7TeV_B1_postLS1_veryrelaxed.dat',
    	'_Allthemachine_7TeV_B1_postLS1_relaxed.dat','_Allthemachine_7TeV_B1_postLS1_baseline.dat','_Allthemachine_7TeV_B1_postLS1_baseline.dat']
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    # broad-band model
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());

    # compute model for collimators
    param_filename_coll=path_here+"Coll_settings/collgaps_fromRoderik_modifNico_materialnames"+scenario+".dat";
    beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;
    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,
	comment=scenario,dire='Coll'+scenario+'/');

    # compute the rest
    param_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_param.dat"
    beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat"
    imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
	    wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,comment='_'+Estr,dire='Rest_'+Estr+'/');

    # add up
    imp_mod_tot=[];wake_mod_tot=[];
    add_impedance_wake(imp_mod_tot,imp_mod_coll,1,1);
    add_impedance_wake(wake_mod_tot,wake_mod_coll,1,1);
    add_impedance_wake(imp_mod_tot,imp_mod_rest,1,1);
    add_impedance_wake(wake_mod_tot,wake_mod_rest,1,1);
    add_impedance_wake(imp_mod_tot,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod_tot,wake_mod_BB,1,1);

    # plot and compare with zbase
    compare_imp_vs_zbase(imp_mod_tot,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

    tuneshiftQp=np.zeros((len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    powloss=np.zeros((len(coatingscan),len(Mscan),len(Nbscan)));
    kickfactor=np.zeros(len(coatingscan));
   
    for icoat,coat in enumerate(coatingscan):
    
	# compute model for TCDQ
	strcoat='_coatCu_'+float_to_str(1e6*coat)+'mum';
	imp_mod_dev=[];wake_mod_dev=[];
	
	if (icoat==0):
	    imp_mod,wake_mod=LHC_element_iw_model(name,[material],radius,radius,E,length,'H',
	    	thickness=[np.inf],wake_calc=wake_calc,lxplusbatch=lxplusbatchImp,
		comment='_coatCu_0mum',dire='test'+model+'/')
	else:
	    imp_mod,wake_mod=LHC_element_iw_model(name,['Cu300K',material],radius,radius,E,length,'H',
	    	thickness=[coat,np.inf],wake_calc=wake_calc,lxplusbatch=lxplusbatchImp,
		comment=strcoat,dire='test'+model+'/')
	
    	# take into account beta functions
	add_impedance_wake(imp_mod_dev,imp_mod,betaratio,betaratio);
    	add_impedance_wake(wake_mod_dev,wake_mod,betaratio,betaratio);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    filemodel=open('impedances'+model+strcoat+'.txt','w');
	    pick.dump(imp_mod_dev,filemodel);
	    filemodel.close();

	    # plot and compare with post-LS1 total impedance
	    if (np.mod(icoat,3)==0):
		maxratio=plot_compare_imp_model([imp_mod_tot,imp_mod_dev],
	    	    ['full post-LS1 impedance ('+scenario+')','DUT, '+material+', radius='+str(radius*1e3)+' mm, L='+str(length)+' m, Cu coating '+str(coat*1e6)+' $\mu $m'],
	    	    listcomp=['Zlong','Zxdip','Zydip'],
		    saveimp=root_result+'/plot_imp_'+machine+"postLS1_"+scenario+'_'+name+strcoat,
	            saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_"+scenario+'_'+name+strcoat);

	    # compute transverse kick factor and power loss 
	    kickfactor[icoat]=transverse_kick_factor(imp_mod_dev,sigmaz,powerspectrum='gaussian');
	    for iM,M in enumerate(Mscan):
		for iNb,Nb in enumerate(Nbscan):
		    powloss[icoat,iM,iNb]=power_loss(imp_mod_dev,sigmaz,gamma,Nb,M,2*np.pi*R,powerspectrum='gaussian');


	    # DELPHI computation
	    for iplane,plane in enumerate(['x','y']):
	        # select Zxdip or Zydip
		for iw in imp_mod_dev:
		    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);

		for iM,M in enumerate(Mscan):

		    flag_trapz=0; # by default no trapz method for DELPHI

        	    if (M==1): nxscan=np.array([0]);flag_trapz=1;
		    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,100),np.arange(M/2-10,M/2+11),
	    		np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

		    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);

		    tuneshiftQp[icoat,iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[icoat,iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			    nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1.e-7,flagm0=True,
			    lxplusbatch=lxplusbatchDEL,comment=machine+strcoat+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
			    dire=root_result+'/');
		    


    # now the plots (outside loop on coating)
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	for iplane,plane in enumerate(['x','y']):

	    for iM,M in enumerate(Mscan):

		for idamp,damp in enumerate(dampscan):

		    for iNb,Nb in enumerate(Nbscan):

			# initialize plots vs Qp
			figQp=[];axQp=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

			# output file name for plots vs Qp
			fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			strpart=['Re','Im'];
			for ir,r in enumerate(['real','imag']):

			    for icoat,coat in enumerate(coatingscan):

				# output file name for data vs Qp
				strcoat='_coatCu_'+str(1e6*coat)+'mum';
				fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strcoat+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    #    if flagcompute:
				ts=getattr(tuneshiftQp[icoat,iplane,iM,:,idamp,iNb,0,0,0],r);
				data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				#else:
			    #	s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
			    #	Qpscan=s[:,0];ts=s[:,1];

				sgn=1;sgnstr='';
				if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				plot(Qpscan,sgn*ts,"CFC with Cu coating "+str(coat*1e6)+" $ \mu $ m",'',"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			# finish plot vs Qp
			for ir,r in enumerate(['real','imag']):
			    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r,fontsize=20)
			    else: end_figure(figQp[ir],axQp[ir]);

	
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	# plot transverse kick factor vs coating thickness
	fig,ax=init_figure();
	plot(coatingscan*1e6,-kickfactor,'Stainless steel with Cu coating','-b',r" $ \kappa_\perp $ [V/(C.m)]",ax,3,xlab="Cu thickness [$ \mu $ m]");
	plot(coatingscan*1e6,-np.ones(len(coatingscan))*kickfactor[0],'Bare stainless steel','--r',r" $ \kappa_\perp $ [V/(C.m)]",ax,3,xlab="Cu thickness [$ \mu $ m]");
	if flagsave: end_figure(fig,ax,save=root_result+'/plot_kickfactor_vs_coat_'+machine)
	else: end_figure(fig,ax);

	for iM,M in enumerate(Mscan):

	    for iNb,Nb in enumerate(Nbscan):
	    
		# plot power loss vs coating thickness
		fileoutplotpow=root_result+'/plot_power_loss_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11';

		figpow,axpow=init_figure();
		plot(coatingscan*1e6,np.abs(powloss[:,iM,iNb]),"Stainless steel with Cu coating",'-b',"Power loss [W]",axpow,3,xlab="Cu thickness [ $ \mu $ m]");
		plot(coatingscan*1e6,np.ones(len(coatingscan))*np.abs(powloss[0,iM,iNb]),"Bare stainless steel",'--r',"Power loss [W]",axpow,3,xlab="Cu thickness [ $ \mu $ m]");
		if flagsave: end_figure(figpow,axpow,save=fileoutplotpow)
		else: end_figure(figpow,axpow);

		for iplane,plane in enumerate(['x','y']):

		    for idamp,damp in enumerate(dampscan):

			# plot tuneshift of m=0 mode at Q'=0 vs coating thickness

			# output file name
			fileoutplottu=root_result+'/plot_tuneshift_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

		    	figtu,axtu=init_figure();
			#print np.squeeze(np.real(tuneshiftm0Qp[pylab.mlab.find(Qpscan==0),:,iplane,iM,idamp,iNb]))
			iQp=pylab.mlab.find(Qpscan==0);#print np.real(tuneshiftm0Qp[:,iplane,iM,iQp,idamp,iNb]);
			plot(coatingscan*1e6,np.real(tuneshiftm0Qp[:,iplane,iM,iQp,idamp,iNb,0,0]),"Stainless steel with Cu coating",'-b',"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			plot(coatingscan*1e6,np.ones(len(coatingscan))*np.real(tuneshiftm0Qp[0,iplane,iM,iQp,idamp,iNb,0,0]),"Bare stainless steel",'--r',"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			if flagsave: end_figure(figtu,axtu,save=fileoutplottu)
			else: end_figure(figtu,axtu);


    if not(flagsave): pylab.show();
