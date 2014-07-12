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

    E=450e9;#E=6500e9;
    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
    flagSach=(lxplusbatchDEL=='launch'); # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    #flagSach=0; # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    
    kmax=2; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of converged eigenvalues (kmax most unstable ones are converged)
    suffix='';#suffix='_only_TCSG_IR7' # suffix for output files 
    lmaxSach=3; # max headtail mode for Sacherer  
    
    model='_TDI_new_coating';
    name='TDI_new_coating';
    ftypescan=2;fminrefine=8.e8;fmaxrefine=2.e12;nrefine=40000;nflog=100;nrefine=2000;

    # scan definitions
    coatingmat_scan=[['Cu300K_in_TDI','Ti_in_TDI'],
    	['NEG_in_TDI','Cu300K_in_TDI','NEG_in_TDI','Ti_in_TDI'],
    	['NEG','Cu300K_in_TDI','NEG','Ti_in_TDI'],
	['NEG_in_TDI','Cu300K_in_TDI','NEG_in_TDI']];
    coatingthickness_scan=[np.array([1,5])*1e-6,np.array([1,1,0.3,5])*1e-6,
    	np.array([1,1,0.3,5])*1e-6,np.array([1,1,5])*1e-6]
    leg=[r"no NEG",r"NEG with $ \rho=6\mu\Omega. $ m",
    	r"NEG with $ \rho=1\mu\Omega. $ m",r"no Ti, NEG with $ \rho=6\mu\Omega. $ m"];
    strmat=["no_NEG","NEG_in_TDI","NEG","no_Ti_NEG_in_TDI"];
    coatingscan=1e-6*np.array([0,1,2,5,10]);
    Qpscan=np.arange(-10,21,1);
    dampscan=np.array([0,0.02]); # damper gain scan
    #Nbscan=np.arange(1.e10,5.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.15,1.7,2.2,3.5])*1e11; # intensity scan
    Mscan=np.array([1,1782,3564]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
        
    # compute the total impedance
    scenarioscan=['_2012','_nominal'];
    iscenario=int(E>=6.5e12);scenario=scenarioscan[iscenario];
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=['_Allthemachine_'+float_to_str(round(E/1e9))+'GeV_B1.dat']
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    # broad-band model
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());

    # compute model for collimators
    sqstr='';
    if E==4e12: sqstr='sq0p6_';
    elif E>=6.5e12: sqstr='sq0p55_';
    if E<=4e12: param_filename_coll=path_here+"Coll_settings/coll_ph1_beta_"+float_to_str(round(E/1e9))+"GeV_"+sqstr+"b1_2012.txt";
    elif E>=6.5e12: param_filename_coll=path_here+"Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat";
    beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;
    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,
	lxplusbatch=lxplusbatchImp,comment=scenario,dire='Coll'+scenario+'/');

    namesref,material,angle,length,halfgap,betax,betay=read_coll_files(param_filename_coll,settings_filename_coll,beta_filename_coll,namesref=['TDI']);
    if (len(material)!=1): print 'Pb: not only one TDI...';sys.exit();
    print halfgap;
    length=[2.8];
    if E==450e9: halfgap=[3.79e-3];
    else: halfgap=[55e-3];
    print namesref,material,angle,length,halfgap,betax,betay
    root_result=path_here+'../../../DELPHI_results/'+machine+'/TDI_'+float_to_str(1000*halfgap[0])+'mm_'+Estr;
    os.system("mkdir -p "+root_result);

    # compute the rest
    if E==450e9: sqstr='sq11m_10m_11m_10m';
    elif E==4e12: sqstr='sq0p6m_3m_0p6m_3m';
    elif E>=6.5e12: sqstr='sq0p55m_10m_0p55m_10m';
    param_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_param.dat"
    beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_"+sqstr+".dat"
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
    #compare_imp_vs_zbase(imp_mod_tot,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
    #	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

    tuneshiftQp=np.zeros((len(coatingmat_scan),len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(coatingmat_scan),len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    tuneshiftQpSach=np.zeros((len(coatingmat_scan),len(coatingscan),2,len(Mscan),len(Qpscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
    tuneshiftQpSachm0=np.zeros((len(coatingmat_scan),len(coatingscan),2,len(Mscan),len(Qpscan),len(Nbscan),1),dtype=complex);
    powloss=np.zeros((len(coatingmat_scan),len(coatingscan),len(Mscan),len(Nbscan)));
    kickfactor=np.zeros((len(coatingmat_scan),len(coatingscan)));
   
    # load Sacherer tuneshifts (when part of the scan was already done)
    #fileSach=open(root_result+'/Sacherer.txt','r');
    #tuneshiftQpSach=pick.load(fileSach);
    #tuneshiftQpSachm0=pick.load(fileSach);
    #fileSach.close();
    #tuneshiftQpSach=np.concatenate((tuneshiftQpSach,np.zeros((1,len(coatingscan),2,len(Mscan),len(Qpscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex)));
    #tuneshiftQpSachm0=np.concatenate((tuneshiftQpSachm0,np.zeros((1,len(coatingscan),2,len(Mscan),len(Qpscan),len(Nbscan),1),dtype=complex)));
    
    for imat,coatingmat in enumerate(coatingmat_scan):
    	#imat+=3;print "imat=",imat
	for icoat,coat in enumerate(coatingscan):

	    # compute model for hBN block of TDI with new coating
	    strcoat='_coatCu_'+float_to_str(1e6*coat)+'mum_'+strmat[imat];
	    coatingthickness=coatingthickness_scan[imat];
	    if (imat!=0): coatingthickness[1]=coat;
	    else: coatingthickness[0]=coat;
	    print coatingmat,coatingthickness
	    imp_mod_dev=[];wake_mod_dev=[];

	    if (coat==0):
		imp_mod,wake_mod=LHC_singlecoll_iw_model(name,material[0],halfgap[0],angle[0],gamma,length[0],
		    wake_calc=wake_calc,coatingmat='Ti_in_TDI',coatingthickness=5e-6,
		    fpar=freq_param(ftypescan=ftypescan,nflog=nflog,fminrefine=fminrefine,
		    fmaxrefine=fmaxrefine,nrefine=nrefine),lxplusbatch=lxplusbatchImp,
		    comment='_coatCu_0mum',dire='test'+model+'/',queue='2nd')
	    else:
		imp_mod,wake_mod=LHC_singlecoll_iw_model(name,material[0],halfgap[0],angle[0],gamma,length[0],
		    wake_calc=wake_calc,coatingmat=coatingmat,coatingthickness=coatingthickness,
		    fpar=freq_param(ftypescan=ftypescan,nflog=nflog,fminrefine=fminrefine,
		    fmaxrefine=fmaxrefine,nrefine=nrefine),lxplusbatch=lxplusbatchImp,
		    comment=strcoat,dire='test'+model+'/',queue='2nd')

    	    # take into account beta functions
	    add_impedance_wake(imp_mod_dev,imp_mod,betax[0]/avbetax,betay[0]/avbetay);
    	    add_impedance_wake(wake_mod_dev,wake_mod,betax[0]/avbetax,betay[0]/avbetay);

	    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

		# dump into a file
		filemodel=open('impedances'+model+strcoat+'.txt','w');
		pick.dump(imp_mod_dev,filemodel);
		filemodel.close();

		# plot and compare with post-LS1 total impedance
		maxratio=plot_compare_imp_model([imp_mod_tot,imp_mod_dev],
	    	    ['full post-LS1 impedance ('+scenario+')',"DELPHI, TDI with new Cu coating "+str(coat*1e6)+" $ \mu $ m, "+leg[imat]],
	    	    listcomp=['Zlong','Zxdip','Zydip'],
		    saveimp=root_result+'/plot_imp_'+machine+"postLS1_"+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+name+strcoat,
	            saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_"+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+name+strcoat);

		# compute transverse kick factor and power loss 
		kickfactor[imat,icoat]=transverse_kick_factor(imp_mod_dev,sigmaz,powerspectrum='gaussian');
		for iM,M in enumerate(Mscan):
		    for iNb,Nb in enumerate(Nbscan):
			powloss[imat,icoat,iM,iNb]=power_loss(imp_mod_dev,sigmaz,gamma,Nb,M,2*np.pi*R,powerspectrum='gaussian');
			print coatingmat,", Cu",coat," mum, M=",M,", Nb=",Nb,", powerloss=",powloss[imat,icoat,iM,iNb],"W"

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
			tuneshiftnxSach=np.zeros((len(Qpscan),len(nxscan),len(Nbscan),1,2*lmaxSach+1),dtype=complex);
			ZeffSach=np.zeros((len(Qpscan),len(nxscan),1,2*lmaxSach+1),dtype=complex);

			tuneshiftQp[imat,icoat,iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[imat,icoat,iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
				nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
				a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
				flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
				kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1.e-7,flagm0=True,
				lxplusbatch=lxplusbatchDEL,comment=machine+strcoat+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
				queue='1nd',dire=root_result+'/');

			if flagSach:
			    # Sacherer (no damper)
			    tuneshiftQpSach[imat,icoat,iplane,iM,:,:,:,:],tuneshiftnxSach,tuneshiftQpSachm0[imat,icoat,iplane,iM,:,:,:],ZeffSach=sacherer(imp_mod,
				    Qpscan,nxscan,Nbscan,[omegas],M,omega0,eval('Q'+plane),gamma,eta,taub,lmaxSach,
				    particle='proton',modetype='sinusoidal',compname='Z'+plane+'dip');

    if flagSach:
	# save Sacherer tuneshifts
	fileSach=open(root_result+'/Sacherer.txt','w');
	pick.dump(tuneshiftQpSach,fileSach);
	pick.dump(tuneshiftQpSachm0,fileSach);
	fileSach.close();
    else:
	# load Sacherer tuneshifts
	fileSach=open(root_result+'/Sacherer.txt','r');
	tuneshiftQpSach=pick.load(fileSach);
	tuneshiftQpSachm0=pick.load(fileSach);
	fileSach.close();


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

			    for imat,coatingmat in enumerate(coatingmat_scan):

				for icoat,coat in enumerate(coatingscan):

				    # output file name for data vs Qp
				    #strcoat='_coatCu_'+str(1e6*coat)+'mum';
	    			    strcoat='_coatCu_'+float_to_str(1e6*coat)+'mum_'+strmat[imat];
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strcoat+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				#    if flagcompute:
				    ts=getattr(tuneshiftQp[imat,icoat,iplane,iM,:,idamp,iNb,0,0,0],r);
				    Sachstr='';
				    if damp==0:
					# compare with Sacherer most unstable mode
					tsSach=getattr(tuneshiftQpSach[imat,icoat,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
					Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));

				    else:
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));

				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift"+Sachstr)
				    #else:
				#	s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				#	Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    if (coat!=0):
				    	plot(Qpscan,sgn*ts,"DELPHI, TDI with new Cu coating "+str(coat*1e6)+" $ \mu $ m, "+leg[imat],col[icoat]+mark[imat],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
				    	if damp==0: plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, with Cu coating '+str(coat*1e6)+" $ \mu $ m, "+leg[imat],'--'+col[icoat]+mark[imat],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
				    elif (imat==0):
				    	plot(Qpscan,sgn*ts,"DELPHI, old TDI",col[icoat]+mark[imat],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
				    	if damp==0: plot(Qpscan,np.squeeze(sgn*tsSach),"Sacherer, old TDI",'--'+col[icoat]+mark[imat],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			# finish plot vs Qp
			for ir,r in enumerate(['real','imag']):
			    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r,fontsize=20)
			    else: end_figure(figQp[ir],axQp[ir]);

	
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	# plot transverse kick factor vs coating thickness
	fig,ax=init_figure();
	for imat,coatingmat in enumerate(coatingmat_scan):
	    plot(coatingscan*1e6,-kickfactor[imat,:],'TDI with new Cu coating, '+leg[imat],'-'+col[imat]+mark[imat],r" $ \kappa_\perp $ [V/(C.m)]",ax,3,xlab="Cu thickness [$ \mu $ m]");
	    if (imat==0): plot(coatingscan*1e6,-np.ones(len(coatingscan))*kickfactor[0,0],'Old TDI','--k',r" $ \kappa_\perp $ [V/(C.m)]",ax,3,xlab="Cu thickness [$ \mu $ m]");
	if flagsave: end_figure(fig,ax,save=root_result+'/plot_kickfactor_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV')
	else: end_figure(fig,ax);

	for iM,M in enumerate(Mscan):

	    for iNb,Nb in enumerate(Nbscan):
	    
		# plot power loss vs coating thickness
		fileoutplotpow=root_result+'/plot_power_loss_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11';

		figpow,axpow=init_figure();
		for imat,coatingmat in enumerate(coatingmat_scan):
		    plot(coatingscan*1e6,np.abs(powloss[imat,:,iM,iNb]),"TDI with new Cu coating, "+leg[imat],'-'+col[imat]+mark[imat],"Power loss [W]",axpow,3,xlab="Cu thickness [ $ \mu $ m]");
		    if (imat==0): plot(coatingscan*1e6,np.ones(len(coatingscan))*np.abs(powloss[0,0,iM,iNb]),"Old TDI",'--k',"Power loss [W]",axpow,3,xlab="Cu thickness [ $ \mu $ m]");
		if flagsave: end_figure(figpow,axpow,save=fileoutplotpow)
		else: end_figure(figpow,axpow);

		for iplane,plane in enumerate(['x','y']):

		    for idamp,damp in enumerate(dampscan):

			# plot tuneshift of m=0 mode at Q'=0 vs coating thickness

			for imat,coatingmat in enumerate(coatingmat_scan):

			    # output file name
			    fileoutplottu=root_result+'/plot_tuneshiftm0_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+strmat[imat]+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutdatatu=root_result+'/data_tuneshiftm0_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strmat[imat]+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    
			    iQp=pylab.mlab.find(Qpscan==0)
			    
			    ts=np.real(tuneshiftm0Qp[imat,:,iplane,iM,iQp,idamp,iNb,0,0])
			    Sachstr='';
			    if damp==0:
				# compare with Sacherer m=0 mode
				tsSach=np.real(tuneshiftQpSachm0[imat,:,iplane,iM,iQp,iNb,0])
				Sachstr=" Sacherer_real_tuneshift_atQp0"
				data=np.hstack((coatingscan.reshape((-1,1))*1e6,ts.reshape((-1,1)),tsSach.reshape((-1,1))));

			    else:
				data=np.hstack((coatingscan.reshape((-1,1))*1e6,ts.reshape((-1,1))));
			    
			    write_ncol_file(fileoutdatatu+'_real.dat',data,header="Coating_thickness[mum]\treal_tuneshift_atQp0"+Sachstr)

			    #for icoat,coat in enumerate(coatingscan):
			    #	print coatingmat,", Cu",coat," mum, M=",M,", Nb=",Nb,", plane",plane,", damp=",damp,", tuneshift (m=0, Qp=0) (DELPHI)=",np.real(tuneshiftm0Qp[imat,icoat,iplane,iM,iQp,idamp,iNb,0,0]);
		    	    
			    figtu,axtu=init_figure();
			    #print np.squeeze(np.real(tuneshiftm0Qp[pylab.mlab.find(Qpscan==0),:,iplane,iM,idamp,iNb]))
			    iQp=pylab.mlab.find(Qpscan==0);#print np.real(tuneshiftm0Qp[:,iplane,iM,iQp,idamp,iNb]);
			    plot(coatingscan*1e6,np.squeeze(ts),"DELPHI, TDI with new Cu coating, "+leg[imat],'-'+col[imat]+mark[imat],"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			    if (imat==0): plot(coatingscan*1e6,np.ones(len(coatingscan))*np.real(tuneshiftm0Qp[0,0,iplane,iM,iQp,idamp,iNb,0,0]),"Old TDI",'-k',"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			    
			    if damp==0:
				plot(coatingscan*1e6,np.squeeze(tsSach),'Sacherer, with Cu coating, '+leg[imat],'--'+col[imat]+mark[imat],"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
				if (imat==0): plot(coatingscan*1e6,np.ones(len(coatingscan))*np.real(tuneshiftQpSachm0[0,0,iplane,iM,iQp,iNb,0]),"Old TDI",'--k',"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			    	#for icoat,coat in enumerate(coatingscan):
				#    print coatingmat,", Cu",coat," mum, M=",M,", Nb=",Nb,", plane",plane,", damp=",damp,", tuneshift (m=0, Qp=0) (Sacherer)=",np.real(tuneshiftQpSachm0[imat,icoat,iplane,iM,iQp,iNb,0]);
			    
			    if flagsave: end_figure(figtu,axtu,save=fileoutplottu)

			    else: end_figure(figtu,axtu);


    if not(flagsave): pylab.show();
