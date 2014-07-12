#!/usr/bin/python

import sys
if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print lxplusbatchImp,lxplusbatchDEL;   

from string import *
import time# as ti
import numpy as np
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from tables_lib import select_in_table
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param
from LHC_imp import *
from HLLHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_new_postLS1';
    #root_result='../../scratch0/LHC_new_postLS1';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=False; # to plot impedances
    nevery=2; # downsampling of the impedance (take less points than in the full model)

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)

    
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition
    subscan=np.array([0,1,2,3,4]);strsubscan='_LHC_v2_all_postLS1_options';margin_factor=1;
    subscan=np.array([0,1,2,4]);strsubscan='_LHC_v2_all_postLS1_options_wo_nominal';margin_factor=1;
    subscan=np.array([2,4,5,6,7]);strsubscan='_LHC_v2_postLS1_Roderik_options';margin_factor=1;
    subscan=np.array([0,3,5,6]);strsubscan='_LHC_v2_postLS1_Roderik_options_nominal';margin_factor=1;
    #subscan=np.array([2]);strsubscan='_LHC_v2_mm_kept';margin_factor=1;
    scenarioscan=np.array(['_2012_v2','_feb2014_relaxed','_feb2014_mm_kept','_feb2014_nominal',
    	'_feb2014_sigma_kept','_may2014_mm_kept_IR3_15_sig_TCL6_open','_may2014_sigma_kept_TCL6_open',
	'_may2014_sigma_kept_TCL6_open_TCSG_12_sigma']);
    dircollscan=np.array(['_2012_v2','_relaxed','_mm_kept','_nominal','_sigma_kept',
    	'_mm_kept','_sigma_kept','_sigma_kept'])
    legscen=np.array(["2012","relaxed settings","tight settings (mm kept)","nominal settings",
    	"tight settings ($ \sigma $ kept)","2012 mm kept",
	"2 sigma retraction","2 sigma retraction with 12 $ \sigma $ TCS IR7"]);
    squeezescan=np.array(['0p6m_3m_0p6m_3m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m',
    	'0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m']);
    Escan=np.array([4e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12]);
    
    param_filename_coll_root=path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames';
   
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_4TeV_B1_physics_fill_3265.dat'])
    #Qpscan=np.arange(-20,31);
    Qpscan=np.array([14,15,16]);
    Qpaver=np.array([14,15,16]);
    #Qpscanplot=np.array([0,15]);
    Qpscanplot=np.array([15]);Qpscanplot=[];
    iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver

    dampscan=np.array([0, 0.001, 0.005, 0.01, 0.02]); # damper gain scan
    #dampscan=np.array([0, 0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.3e11,1.5e11,1.7e11,3.e11,4.e11]); # intensity scan for plot vs Qp
    Nbscanplot=[]; # intensity scan for plot vs Qp
    NbscanHEADTAILv1=np.array([1.7e11]);
    NbscanHEADTAILv2=np.array([1.5e11]);
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((len(subscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(subscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    #t1=ti.clock()
    for iscenario,scenario in enumerate(scenarioscan[subscan]):
    	
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	print "scenario: ",scenario
	
	if scenario.startswith('_2012'):
	    param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
	    settings_filename_coll=path_here+'Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';
	else:
	    param_filename_coll=param_filename_coll_root+scenario+'.dat';
	    settings_filename_coll=param_filename_coll;
	
	# compute imp. model
	scenariobis=scenario.replace('_feb2014','');
	imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
		settings_filename_coll,dire=path_here+"LHC_elements/",commentcoll=scenariobis,
		direcoll='Coll'+dircollscan[subscan[iscenario]]+'/',
		lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
		wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario)
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    #filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    #pick.dump(imp_mod,filemodel);
	    #filemodel.close();
	    
	    # write Ascii files with each component
	    write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
	    	dire=root_result+'/')
	    
	    # plot and compare with zbase
	    #if subscan[iscenario]<=1: compare_imp_vs_zbase(imp_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[subscan[iscenario]],
	    #	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

	    if (wake_calc):
		# write Ascii files with each component
		#write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	#    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
		#    dire=root_result+'/')
	    
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# dip only
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)		

    
    #t2=ti.clock();
    #print "time for impedance computation [seconds]: ",t2-t1;
    
    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	#t1=ti.clock()
	# plot and compare all scenarios
	maxratio_sb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	    ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3);
	
	
	maxratio_mb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	    ylim=[1e5,1e10],bounds=[8e3,40e6],legpos=3);

	# plot and compare all scenarios with 2012 wake
	if wake_calc:
	    maxratio_w=plot_compare_imp_model(wake_mod_list,legscen[subscan],listcomp=['Wlong','Wxdip','Wydip'],
		saveimp=root_result+'/plot_wake_'+machine+"_scenarios"+strsubscan,
		saveratio=root_result+'/plot_wake_ratio_'+machine+"_scenarios"+strsubscan,
		xlim=[1e-5,1e6],ylim=[1e8,1e19],yliml=[1e6,1e19],bounds=[8e3,5e10],legpos=0);

	#t2=ti.clock();
	#print "time for impedance plots [seconds]: ",t2-t1;
        
	#t1=ti.clock()
	# DELPHI scans now
	for iscenario,scenario in enumerate(scenarioscan[subscan]):

	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
	    
	    # write HEADTAIL input file
	    # first do the longitudinal matching
	    delta_lin,eps0_lin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=True);
	    print "Longitudinal linear matching: delta_p/p=",delta_lin,", eps0=",eps0_lin,", Qs=",Qs;
	    delta_nlin,eps0_nlin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=False);
	    print "Longitudinal non-linear matching: delta_p/p=",delta_nlin,", eps0=",eps0_nlin,", Qs=",Qs;
	    # check that Qs is correct
	    if ((Qs-Qs_from_RF_param(V,h,gamma,eta))>1e-8): print " Pb in Qs !!";sys.exit();
	    
	    # write .cfg files (linear and non-linear matching)
	    if len(Nbscanplot)>0: 
		cfgnamelin=root_result+'/'+machine+"_"+Estr+scenario+'_lin.cfg';
		write_CFG_HEADTAIL(cfgnamelin,Nb=Nbscanplot[0],betax=R/Qx,betay=R/Qy,
	    	    sigmaz=sigmaz,delta=delta_lin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		    gamma=gamma,nturns=200000,Qx=Qx,Qy=Qy,isyn=1,start_turn=199000,end_turn=199100,VRF=V);

		cfgnamenlin=root_result+'/'+machine+"_"+Estr+scenario+'_nlin.cfg';
		write_CFG_HEADTAIL(cfgnamenlin,Nb=Nbscanplot[0],betax=R/Qx,betay=R/Qy,
	    	    sigmaz=sigmaz,delta=delta_nlin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		    gamma=gamma,nturns=200000,Qx=Qx,Qy=Qy,isyn=4,start_turn=199000,end_turn=199100,VRF=V);
	    
	    # DELPHI computation
	    tuneshiftQp[iscenario,:,:,:,:,:,:,:,:],tuneshiftm0Qp[iscenario,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_list[iscenario],
	    	Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-4,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV',
		queue='2nd',dire=root_result+'/',flagQpscan_outside=True);

	#t2=ti.clock();
	#print "time for DELPHI computations [seconds]: ",t2-t1;
	
	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    #t1=ti.clock()
		    # plots vs Q'
		    for idamp,damp in enumerate(dampscan):

			for Nb in Nbscanplot:

			    # initialize plots vs Q'
			    #t1bis=ti.clock();
			    figQp=[];axQp=[];figQpm0=[];axQpm0=[];
			    for ir in range(2):
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);
			    #t2bis=ti.clock();
			    #print "M=",M,plane,"damp=",damp,",Nb=",Nb,", time for plots vs Qp, init_figure [seconds]: ",t2bis-t1bis;

			    # output files name for plots vs Q'
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    #t1bis=ti.clock();
			    strpart=['Re','Im'];
			    for ir,r in enumerate(['real','imag']):

				for iscenario,scenario in enumerate(scenarioscan[subscan]):

				    # output files name for data vs Qp
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    
				    tsm0=getattr(tuneshiftm0Qp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
				    if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate
				    
				    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
				    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],ylab+" (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");

				    if (M==1)and(damp in [0,0.02]):
					# compare with HEADTAIL
				    	EstrHEAD=float_to_str(round(Escan[subscan[iscenario]]/1e11)/10.)+'TeV';
					fact=1;
					#if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0 (conversion growth rate to tune shift)
					
					if (scenario=='_2012')and(Nb in NbscanHEADTAILv1):
					    nsl=500;npr=1000000;nlin=1; # HEADTAIL parameters for comparison
					    if (ir==0): fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_damper_dip_1b_ntwake20_nkick1_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_drate"+float_to_str(damp)+"_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin);
					    sufHEADTAIL="_aver_Sussix_most_tau_finer.txt";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
					
					elif ((scenario=='_2012_v2')or(scenario=='_feb2014_mm_kept'))and(Nb in NbscanHEADTAILv2):
					    nsl=200;npr=500000;nlin=1; # HEADTAIL parameters for comparison
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
					    if (scenario=='_feb2014_mm_kept'): scenariobis='_mm_kept';
					    else: scenariobis=scenario;
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenariobis+"_dip_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_pre1_nsig2_drate"+float_to_str(damp);
					    sufHEADTAIL="_aver_most_tau_finer.txt";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
			    #t2bis=ti.clock();
			    #print "M=",M,plane,"damp=",damp,",Nb=",Nb,", time for plots vs Qp, plot & write data [seconds]: ",t2bis-t1bis;


			    #t1bis=ti.clock()
			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
			        axQp[ir].set_xlim([-10,25]);#axQp[ir].set_ylim([0,1.6e-5]);
				end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
				end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))
			    #t2bis=ti.clock();
			    #print "M=",M,plane,"damp=",damp,",Nb=",Nb,", time for plots vs Qp, end_figure [seconds]: ",t2bis-t1bis;

		    #t2=ti.clock();
		    #print "M=",M,plane,", time for plots vs Qp [seconds]: ",t2bis-t1bis;
			
		    #t1=ti.clock()
		    # plots vs Q' & damper
		    for Nb in Nbscanplot:
		    	
			for iscenario,scenario in enumerate(scenarioscan[subscan]):

			    # initialize plots vs Q'
			    figQp=[];axQp=[];figQpm0=[];axQpm0=[];
			    for ir in range(2):
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

			    # output files name for plots vs Q'
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+scenario+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+scenario+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];
			    for idamp,damp in enumerate(dampscan):

			    	for ir,r in enumerate(['real','imag']):

				    # output files name for data vs Qp
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    
				    tsm0=getattr(tuneshiftm0Qp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),tsm0.reshape((-1,1))));
				    write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    #if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
				    if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate
				    
				    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+',d='+str(damp),col[idamp],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
				    plot(Qpscan[~np.isnan(np.squeeze(tsm0))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(tsm0))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+',d='+str(damp),col[idamp],ylab+" (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");

			    # finish plot vs Qp & damper
			    for ir,r in enumerate(['real','imag']):
			        axQp[ir].set_xlim([-15,25]);
				end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
				end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))

		    #t2=ti.clock();
		    #print "M=",M,plane,", time for plots vs Qp & damper [seconds]: ",t2-t1;
			
		    #t1=ti.clock()
		    # plots vs Nb, and TMCI plots
		    for idamp,damp in enumerate(dampscan):

			for Qp in Qpscanplot:

			    # initialize plots vs Nb
			    figNb=[];axNb=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figNb.append(fig);axNb.append(ax);
			    figNbm0=[];axNbm0=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figNbm0.append(fig);axNbm0.append(ax);
			    # initialize TMCI plots
			    figTMCI=[];axTMCI=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figTMCI.append(fig);axTMCI.append(ax);
			    
			    pat=['.',''];

			    # output file name for plots vs Q'
			    fileoutplotNb=root_result+'/plot_vs_Nb_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotNbm0=root_result+'/plot_vs_Nb_m0_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];ylim=([-5,3],[-0.01,1]);
			    for ir,r in enumerate(['real','imag']):

				for iscenario,scenario in enumerate(scenarioscan[subscan]):

				    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
				    # output file name for data vs Nb
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataNb=root_result+'/data_vs_Nb_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutdataNbm0=root_result+'/data_vs_Nb_m0_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

				    ts_most=getattr(tuneshiftQp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,0],r);
				    data=np.hstack((Nbscan.reshape((-1,1)),ts_most.reshape((-1,1))));
				    #write_ncol_file(fileoutdataNb+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_most_unstable")

				    tsm0=getattr(tuneshiftm0Qp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0],r);
				    data=np.hstack((Nbscan.reshape((-1,1)),tsm0.reshape((-1,1))));
				    #write_ncol_file(fileoutdataNbm0+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_m0")

				    sgn=1;sgnstr='';
				    #if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
				    if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate
				    
				    plot(Nbscan[~np.isnan(np.squeeze(ts_most))],np.squeeze(sgn*ts_most)[~np.isnan(np.squeeze(ts_most))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],ylab,axNb[ir],0,xlab="Intensity [p+/bunch]");
				    plot(Nbscan[~np.isnan(np.squeeze(tsm0))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(tsm0))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],ylab+" (mode 0)",axNbm0[ir],0,xlab="Intensity [p+/bunch]");
				    
				    # set plot axes
				    axNb[ir].set_xlim([0,8e11]);
				    maxy=np.ceil(np.max(np.abs(sgn*ts_most))*1e5)/1e5;
				    if ir==0: axNb[ir].set_ylim([-maxy,0]);
				    else: ylim[ir][1]=np.ceil(maxy/omegas*5)/5.;axNb[ir].set_ylim([0,maxy]);

				    if (False)and(((scenario=='_2012_v2')or(scenario=='_feb2014_mm_kept'))and(M==1)):
					# compare with HEADTAIL
				    	EstrHEAD=float_to_str(round(Escan[subscan[iscenario]]/1e11)/10.)+'TeV';
					fact=1;
					#if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0 (conversion growth rate to tune shift)
					sufHEADTAIL="_aver_most_tau_finer.txt";
					sufHEADTAILm0="_m0_finer.txt";
					if (scenario=='_feb2014_mm_kept'): scenariobis='_mm_kept';
					else: scenariobis=scenario;
					
					if (damp==0)and(Qp==0):
					    # HEADTAIL parameters for comparison
					    nsl=200;npr=500000;
					    nlin=1;nlinstr='_lin';dipstr='_dip';
					    #nlin=4;nlinstr='_nlin';dipstr='';
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotNb=fileoutplotNb+'_vs_HEADTAIL'+nlinstr+dipstr;
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenariobis+dipstr+"_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_csi"+float_to_str(Qp)+"_qsec0_oct0_drate"+float_to_str(damp)+"_pre1_nsig2";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0]*1e11,fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],ylab,axNb[ir],0,xlab="Intensity [p+/bunch]");
					    
					    # tuneshift of mode 0
					    sm0=read_ncol_file(rootHEADTAIL+sufHEADTAILm0,ignored_rows=1);
                                	    plot(sm0[:,0]*1e11,fact*sm0[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],ylab+" (mode 0)",axNbm0[ir],0,xlab="Intensity [p+/bunch]");
					
					elif (damp==0.02)and(Qp==15):
					    # HEADTAIL parameters for comparison
					    nsl=200;npr=500000;
					    nlin=1;nlinstr='_lin';dipstr='_dip';
					    #nlin=4;nlinstr='_nlin';dipstr='';
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotNb=fileoutplotNb+'_vs_HEADTAIL'+nlinstr+dipstr;
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenariobis+dipstr+"_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_csi"+float_to_str(Qp)+"_qsec0_oct0_drate"+float_to_str(damp)+"_pre1_nsig2";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
					    plot(s[:,0]*1e11,fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],ylab,axNb[ir],0,xlab="Intensity [p+/bunch]");
				    
				    # TMCI plot
				    ts=tuneshiftQp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,:];
				    plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,patcol=pat[ir]+col[iscenario],xlab='Nb [p+/b]',
					title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);
	    

			    # finish plots vs Nb and TMCI plots
			    for ir,r in enumerate(['real','imag']):
				end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
				end_figure(figNbm0[ir],axNbm0[ir],save=flagsave*(fileoutplotNbm0+'_'+r))
				end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r));

		    #t2=ti.clock();
		    #print "M=",M,plane,", time for plots vs Nb & TMCI plots [seconds]: ",t2-t1;


    #####################################################
    # plots stabilizing emittance vs Nb for certain Qp
    #####################################################

    # post-LS1 scenarios
    #legbeam=['Classical 25 ns','BCMS 25 ns', 'Classical 50 ns','BCMS 50 ns'];
    #emitbeam=[3.75,1.9,2.2,1.6];
    #intbeam=[1.35,1.15,1.65,1.6];
    #Mbeam=[3564,3564,1782,1782];
    #colbeam=['ok','xk','oc','xc'];
    # from Giovanni Rumolo LBOC talk, 8/4/2014, with 0.6 mm.mrad emittance blow-up in LHC
    # except emittance standard 25ns (3.75 -> nominal)
    legbeam=['25 ns, standard','25 ns, BCMS','25ns, standard, 8b+4e','50 ns, standard (2012)'];
    emitbeam=[3.75,1.9,2.9,2.2];
    intbeam=[1.3,1.3,1.8,1.7];
    Mbeam=[3564,3564,3564,1782];
    colbeam=['ok','xk','+k','dc'];
    # HL-LHC scenarios
    #legbeam=['PIC low emit.','PIC high emit.', 'US1',"US2 low int. & low emit.","US2 high int. & high emit."];
    #emitbeam=[1.8,2.22,2.62,2.26,2.5];
    #intbeam=[1.38,1.38,1.9,1.9,2.2];
    #colbeam=['or','+m','xg','dk','vb'];

    #t1=ti.clock()
    # OLD VERSION
    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    #if flagcompute:
	    # test that first scenario is 2012 updated model, and that any other is an HL-LHC one
	    if (len(subscan)>1)and((subscan[0]==0)and(all(subscan[1:]>0))):
	    
		# parameters
		oct2012signscan=['_pos_oct','_neg_oct'];
		oct2012scan=[510.,250];Nb2012=1.5e11;emit2012=2.5;en2012=4.;
		octLS1=550.;enLS1=6.5; 
		#octHL=550.;enHL=7.;

		colscen=['b','m','k','c'];

		for isign,octsign in enumerate(oct2012signscan):

		    for iM,M in enumerate(Mscan):

			for idamp,damp in enumerate(dampscan):

			    # initialize plot
			    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.7]);

			    # output file name for plots vs emittance
			    fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

			    # plot postLS1 beam parameters
			    for ib,legb in enumerate(legbeam):
	  			plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)",ms=16);

			    # with assumption of linearity between intensity and growth rate:
			    # oct*emit/en = impfact * factlin *int;
			    factor2012lin=oct2012scan[isign]*emit2012/(en2012*Nb2012/1e11);
			    # better : oct*emit/en^2 =  fact *imag(-tuneshift(int));
			    #factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
			    # take 50ns case in 2012
			    iM0=pylab.mlab.find(Mscan==1782);
			    imtu0=np.min(getattr(tuneshiftQp[0,:,iM0,iQpaver,idamp,pylab.mlab.find(Nbscan==Nb2012),0,0,0],'imag'));
			    #imtu0=np.max(np.abs(tuneshiftQp[0,:,iM0,iQpaver,idamp,pylab.mlab.find(Nbscan==Nb2012),0,0,0]));
			    print imtu0
			    factor2012=oct2012scan[isign]*emit2012/(en2012**2*abs(imtu0));
			    # emit=fact*en^2/oct *imag(-tuneshift)

			    for iscenario,scenario in enumerate(scenarioscan[subscan[1:]]):

				# output file name for data int. vs emitttance
				Estr=float_to_str(round(Escan[subscan[iscenario+1]]/1e9))+'GeV';
				fileoutdataemit=root_result+'/data_int_vs_emit_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

				imtu=np.zeros(len(Nbscan));
				for iNb,Nb in enumerate(Nbscan):
				    ts=getattr(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0],'imag');
		    		    imtu[iNb]=np.min(ts[~np.isnan(ts)]);
		    		    #imtu[iNb]=np.max(np.abs(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0]));
				print imtu

				#emitLS1=factor2012*enHL**2/octHL * np.abs(imtu);
				emitLS1=factor2012*enLS1**2/octLS1 * np.abs(imtu);
				plot(emitLS1,Nbscan/1e11,"Limit with "+legscen[subscan[iscenario+1]],colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
				# with linear assumption
				#if damp==0: emitLS1lin=np.max(np.abs(maxratio_mb[:,iscenario,0]+0*maxratio_mb[:,iscenario,1]))*factor2012lin*enLS1/octLS1 * (Nbscan/1e11);
				#else: emitLS1lin=np.max(np.abs(maxratio_sb[:,iscenario,0]+0*maxratio_sb[:,iscenario,1]))*factor2012lin*enLS1/octLS1 * (Nbscan/1e11);
				#plot(emitLS1lin,Nbscan/1e11,"Limit with "+legscen[subscan[iscenario+1]]+' '+Estr+', linear approx.','--'+colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
				ax.set_xlim([0,6]);ax.set_ylim([0,6]);

				data=np.hstack((emitLS1.reshape((-1,1)),Nbscan.reshape((-1,1))));
				write_ncol_file(fileoutdataemit+'_'+r+'.dat',data,header="emit\t"+strpart[ir]+"Nb")


			    # finish plot
			    end_figure(fig,ax,save=flagsave*fileoutplotemit,legpos=(0.3,0.5),legfontsize=20)

    #t2=ti.clock();
    #print "time for int. vs. emit plots (OLD) [seconds]: ",t2-t1;


    # NEW VERSION (14/05/2014)
    colscen=['b','r','g','m','k','c','y'];
    colscen2=np.array([(0,0,1),(1,0,0),(0,1,0),(1,0,1),(0,0,0),(0,1,1),(1,1,0)]);
    hatch=np.array(["/","\\","|","X","/","\\","|","X"]);

    oct2012signscan=['_neg_oct','_pos_oct'];
    legsignscan=["negative polarity (\"old\")","positive polarity (\"new\")"];
    octLS1=570.;

    # relevant 2012 instability data - based on files ../Mesures_LHC/instability_table_B2V_pos_oct_flattop.csv &
    # ../Mesures_LHC/instability_table_B2H_pos_oct_flattop.csv (Q' > 10)
    # and file ../Mesures_LHC/instability_table_B2H_neg_oct_squeeze.csv (Q'>5)
    # There are also a 3 flat top instabilities taken from an LMC talk
    # by G. Arduini in August 2012 (slides in ../Docs/lmc_145c_talk_instabilities2012_LHC_Gianluigi.pdf,
    # slide 14&15, fills 2919, 2920 & 2932 - damper gain from Timber & trim editor)
    en2012=4;M_2012=1782;
    # negative octupole polarity
    dataoct_neg2012=np.array([58.,200.,200.]);
    meanoct_neg2012=np.average(dataoct_neg2012);erroct_neg2012=np.sqrt(np.var(dataoct_neg2012));
    dataemit_neg2012=np.array([2.3,2.5,2.5]);
    meanemit_neg2012=np.average(dataemit_neg2012);erremit_neg2012=0.5;
    dataNb_neg2012=np.array([1.4,1.47,1.46])*1e11;
    dataQp_neg2012=np.array([5.9,7.,7.]);
    dataplane_neg2012=np.array(['x','x','x']);
    datadamp_neg2012=1./np.array([50.,100.,100.]);
    # positive octupole polarity
    dataoct_pos2012=np.array([209.,487.,510.,35.,35.,510.]);
    meanoct_pos2012=np.average(dataoct_pos2012);erroct_pos2012=np.sqrt(np.var(dataoct_pos2012));
    dataemit_pos2012=np.array([2.3,2.3,2.3,2.2,2.2,2.5]);
    meanemit_pos2012=np.average(dataemit_pos2012);erremit_pos2012=0.5;
    dataNb_pos2012=np.array([1.64,1.64,1.63,1.64,1.64,1.44])*1e11;
    dataQp_pos2012=np.array([15.4,18.3,10.3,13.6,17.8,9.]);
    dataplane_pos2012=np.array(['x','x','x','x','y','x']);
    datadamp_pos2012=1./np.array([100.,100.,100.,100.,100.,100.]);

    # scan parameters for DELPHI computation at 4 TeV for those experimental cases
    oct_2012=np.hstack((dataoct_neg2012,dataoct_pos2012));
    emit_2012=np.hstack((dataemit_neg2012,dataemit_pos2012));
    Nb_2012=np.hstack((dataNb_neg2012,dataNb_pos2012));
    Qp_2012=np.hstack((dataQp_neg2012,dataQp_pos2012));
    damp_2012=np.hstack((datadamp_neg2012,datadamp_pos2012));
    plane_2012=np.hstack((dataplane_neg2012,dataplane_pos2012));
    ind_neg2012=np.arange(3);ind_pos2012=np.arange(3,9); # indices in the above tables, for resp. neg. and oct. polarities

    # compute first with 2012 imp. model
    param_filename_coll_2012=path_here+'Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
    settings_filename_coll_2012=path_here+'Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';
    # fixed parameters
    machine,E_2012,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=en2012*1e12);
    avbetax_2012=R/Qx;avbetay_2012=R/Qy; # average beta functions used
    # model
    imp_mod_2012,wake_mod_2012=LHC_imp_model_v2(E_2012,avbetax_2012,avbetay_2012,param_filename_coll_2012,
	    settings_filename_coll_2012,TDIcoating='preLS1',dire=path_here+"LHC_elements/",commentcoll='_2012_v2',direcoll='Coll_2012_v2/',
	    lxplusbatch=lxplusbatchImp,BPM=False,beam=beam,squeeze='0p6m_3m_0p6m_3m',
	    wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave='_2012_v2')
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    #t1=ti.clock()
    tuneshiftQp_2012=np.zeros((len(Qp_2012),1,1,1,1,kmaxplot),dtype=complex);
    factor_2012=np.zeros(len(Qp_2012));
    for iQp,Qp in enumerate(Qp_2012):

	# select parameters
	plane=plane_2012[iQp];iplane=int(plane=='y');print plane,iplane
	Nb=Nb_2012[iQp];damp=damp_2012[iQp];
	# select Zxdip or Zydip
	for iw in imp_mod_2012:
	    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);

	flag_trapz=0; # by default no trapz method
	if (M_2012==1): nxscan=np.array([0]);flag_trapz=1;
	else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M_2012,M_2012/20),np.arange(M_2012/2-10,M_2012/2+11),
		np.arange(M_2012-10,M_2012),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

    	tuneshiftnx=np.zeros((1,len(nxscan),1,1,1,1,kmaxplot),dtype=complex);

	tuneshiftQp_2012[iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
	    nxscan,[damp],[Nb],[omegas],[dphase],M_2012,omega0,eval('Q'+plane),
	    gamma,eta,a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,flag_trapz=flag_trapz,
	    flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
	    abseps=1.e-4,flagm0=False,lxplusbatch=lxplusbatchDEL,
	    comment=machine+'_2012_v2_'+float_to_str(round(E_2012/1e9))+'GeV_'+str(M_2012)+'b_Qp'+str(Qp)+'_'+plane,
	    queue='2nd',dire=root_result+'/');

	# "stability factor" for each case
	factor_2012[iQp]=-oct_2012[iQp]*emit_2012[iQp]/(en2012**2*np.imag(tuneshiftQp_2012[iQp,0,0,0,0,0]));
	print "all factors 2012:",factor_2012;
	factor_neg_oct_2012_mean=np.average(factor_2012[ind_neg2012]);
	factor_pos_oct_2012_mean=np.average(factor_2012[ind_pos2012]);
	factor_neg_oct_2012_sig=np.sqrt(np.var(factor_2012[ind_neg2012]));
	factor_pos_oct_2012_sig=np.sqrt(np.var(factor_2012[ind_pos2012]));
	print "averages & sigmas:",factor_neg_oct_2012_mean,factor_neg_oct_2012_sig,factor_pos_oct_2012_mean,factor_pos_oct_2012_sig;

    # end of computations with 2012 model
    #t2=ti.clock();
    #print "time for 2012 computations [seconds]: ",t2-t1;

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    #t1=ti.clock()
	    # plots of averages, without "error bars", all scenarios on same plot
	    for isign,octsign in enumerate(oct2012signscan):

		# assumes oct*emit/en^2 =  fact *imag(-tuneshift(int));
		#factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
		fact=eval('factor'+octsign+'_2012_mean');
		sig_fact=eval('factor'+octsign+'_2012_sig');
		# emit=fact*en^2/oct *imag(-tuneshift)

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			# initialize plot
			figall,axall=init_figure(axes=[0.15,0.1,0.8,0.7]);

			# output file name for plots vs emittance
			fileoutplotemitall=root_result+'/plot_int_vs_emit_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';

			# plot postLS1 beam parameters
			for ib,legb in enumerate(legbeam):
	  		    if (M==1)or(M==Mbeam[ib]):
				plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",axall,0,xlab="Norm. emittance ($ \mu $ m)",ms=16);

			for iscenario,scenario in enumerate(scenarioscan[subscan]):

			    # output file name for data int. vs emitttance
			    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
			    fileoutdataemit=root_result+'/data_int_vs_emit_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];
			    fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

			    # initialize plot with only one scenario and its "error bar" (spread)
			    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.75])
			    # plot postLS1 beam parameters
			    for ib,legb in enumerate(legbeam):
	  			if (M==1)or(M==Mbeam[ib]):
				    plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)",ms=16);

			    imtu=np.zeros(len(Nbscan));
			    for iNb,Nb in enumerate(Nbscan):
				ts=getattr(tuneshiftQp[iscenario,:,iM,iQpaver,idamp,iNb,0,0,0],'imag');
		    		if len(np.isnan(ts))==0: imtu[iNb]=np.min(ts);
				else: imtu[iNb]=np.min(ts[~np.isnan(ts)]);
		    		#imtu[iNb]=np.max(np.abs(tuneshiftQp[iscenario,:,iM,iQpaver,idamp,iNb,0,0,0]));
			    print imtu

			    emit=fact*(Escan[subscan[iscenario]]/1.e12)**2/octLS1 * np.abs(imtu);
			    sig_emit=sig_fact*(Escan[subscan[iscenario]]/1.e12)**2/octLS1 * np.abs(imtu)
			    #ax.errorbar(emit,Nbscan/1e11,xerr=sig_emit,fmt=colscen[iscenario],label="Stab. limit "+legscen[subscan[iscenario]],lw=3);
			    plot(emit,Nbscan/1e11,"Stab. limit "+legscen[subscan[iscenario]],colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",axall,0,lw=3.,xlab="Norm. emittance ($ \mu $ m)");
			    axall.set_xlim([0,6]);axall.set_ylim([0,6]);
			    #i+=1;

			    data=np.hstack((emit.reshape((-1,1)),Nbscan.reshape((-1,1)),sig_emit.reshape((-1,1))));
			    write_ncol_file(fileoutdataemit+'.dat',data,header="emit\tsig_emit\tNb")
			    			    
			    # plot of each individual line with "error bars"
			    plot(emit,Nbscan/1e11,"Stab. limit with "+legsignscan[isign],colscen[isign],"Intensity ($ 10^{11} $p+/b)",ax,0,lw=3.,xlab="Norm. emittance ($ \mu $ m)");
			    ax.fill_betweenx(Nbscan/1e11,emit-sig_emit,emit+sig_emit,color=colscen[isign],lw=2,alpha=0.3)
			    ax.set_xlim([0,6]);ax.set_ylim([0,6]);

			    # finish plot with individual scenario
			    end_figure(fig,ax,save=flagsave*fileoutplotemit,legpos=(0.3,0.7))

			# finish plot
			end_figure(figall,axall,save=flagsave*fileoutplotemitall,legpos=(0.3,0.5),legfontsize=20)

	    
	    # plots of individual scenarios, average & "error bars", both polarities on same plot
	    for iM,M in enumerate(Mscan):

		for idamp,damp in enumerate(dampscan):

		    for iscenario,scenario in enumerate(scenarioscan[subscan]):

			# output file name for data int. vs emitttance
			Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
			fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+'_converged'+strnorm[flagnorm];

			imtu=np.zeros(len(Nbscan));
			for iNb,Nb in enumerate(Nbscan):
			    ts=getattr(tuneshiftQp[iscenario,:,iM,iQpaver,idamp,iNb,0,0,0],'imag');
		    	    if len(np.isnan(ts))==0: imtu[iNb]=np.min(ts);
			    else: imtu[iNb]=np.min(ts[~np.isnan(ts)]);
		    	    #imtu[iNb]=np.max(np.abs(tuneshiftQp[iscenario,:,iM,iQpaver,idamp,iNb,0,0,0]));

			# initialize plot with only one scenario and its "error bar" (spread)
			fig,ax=init_figure(axes=[0.15,0.1,0.8,0.75])
			# plot postLS1 beam parameters
			for ib,legb in enumerate(legbeam):
	  		    if (M==1)or(M==Mbeam[ib]):
				plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)",ms=16);

			for isign,octsign in enumerate(oct2012signscan):

			    # assumes oct*emit/en^2 =  fact *imag(-tuneshift(int));
			    #factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
			    fact=eval('factor'+octsign+'_2012_mean');
			    sig_fact=eval('factor'+octsign+'_2012_sig');
			    # emit=fact*en^2/oct *imag(-tuneshift)

			    emit=fact*(Escan[subscan[iscenario]]/1.e12)**2/octLS1 * np.abs(imtu);
			    sig_emit=sig_fact*(Escan[subscan[iscenario]]/1.e12)**2/octLS1 * np.abs(imtu)
			    
			    # plot of each individual line with "error bars"
			    plot(emit,Nbscan/1e11,"Stab. limit with "+legsignscan[isign],colscen[isign],"Intensity ($ 10^{11} $p+/b)",ax,0,lw=3.,xlab="Norm. emittance ($ \mu $ m)");
			    ax.fill_betweenx(Nbscan/1e11,emit-sig_emit,emit+sig_emit,color=colscen[isign],lw=2,alpha=0.3)
			    #ax.set_xlabel("Norm. emittance ($ \mu $ m)");
			    #ax.set_ylabel("Intensity ($ 10^{11} $p+/b)");
			    ax.set_xlim([0,6]);ax.set_ylim([0,6]);


			# finish plot
			end_figure(fig,ax,save=flagsave*fileoutplotemit,legpos=(0.3,0.7))

	    #t2=ti.clock();
	    #print "time for int. vs. emit plots (NEW) [seconds]: ",t2-t1;

    if not(flagsave): pylab.show();
