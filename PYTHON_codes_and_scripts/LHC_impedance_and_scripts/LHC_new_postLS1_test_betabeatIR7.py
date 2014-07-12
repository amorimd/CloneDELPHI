#!/usr/bin/python

import sys
if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print lxplusbatchImp,lxplusbatchDEL;   


from string import *
import time
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
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_new_postLS1_test_beatbeatIR7';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=True; # to plot impedances
    nevery=2; # downsampling of the impedance (take less points than in the full model)

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    coef='0.2_verynice_beatings';
    collscenario='_prematch'+coef;
    # scan definition
    subscan=np.array([0,1,2,3]);strsubscan='_LHC_v1_v2_old_new_optics'+collscenario;margin_factor=1;
    scenarioscan=np.array(['_nominal','_nominal_new_'+coef,'_nominal_v2','_nominal_v2_new_'+coef]);
    legscen=np.array(['Nominal settings, old imp. model, old optics',
    	'Nominal settings, old imp. model, new optics',
    	'Nominal settings, updated imp. model, old optics',
	'Nominal settings, updated imp. model, new optics']);
    squeezescan=np.array(['0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m']);
    Escan=np.array([7e12,7e12,7e12,7e12]);
    
    direcoll_scan=['Coll_nominal/','Coll_nominal/','Coll_nominal/','Coll_nominal/'];
    param_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat']);
    settings_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/collgaps_settings_nominal_betabetaIR7'+collscenario+'_b'+beam+'.dat',
    	path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/collgaps_settings_nominal_betabetaIR7'+collscenario+'_b'+beam+'.dat']);
    beta_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/coll_ph1_beta_7000GeV_betabeatIR7'+collscenario+'_b'+beam+'.txt',
    	path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_nominal.dat',
    	path_here+'Coll_settings/coll_ph1_beta_7000GeV_betabeatIR7'+collscenario+'_b'+beam+'.txt']);
   
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_4TeV_B1_physics_fill_3265.dat'])
    Qpscan=np.arange(-20,31);
    Qpaver=np.array([14,15,16]);
    Qpscanplot=np.array([0,15]);
    iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver

    dampscan=np.array([0,0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.5e11,1.7e11]); # intensity scan for plot vs Qp
    NbscanHEADTAILv1=np.array([1.7e11]);
    NbscanHEADTAILv2=np.array([1.5e11]);
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    for iscenario,scenario in enumerate(scenarioscan[subscan]):
    	
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	print "scenario: ",scenario
	
	# compute imp. model
	if subscan[iscenario]<2:
	    imp_mod,wake_mod=LHC_imp_model_v1(E,avbetax,avbetay,param_filename_coll_scan[subscan[iscenario]],
		    settings_filename_coll_scan[subscan[iscenario]],beta_filename_coll=beta_filename_coll_scan[subscan[iscenario]],
		    dire=path_here+"LHC_elements/",commentcoll='_nominal',direcoll=direcoll_scan[subscan[iscenario]],
		    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
		    wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario)
	else:
	    imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll_scan[subscan[iscenario]],
		    settings_filename_coll_scan[subscan[iscenario]],beta_filename_coll=beta_filename_coll_scan[subscan[iscenario]],
		    dire=path_here+"LHC_elements/",commentcoll='_nominal',direcoll=direcoll_scan[subscan[iscenario]],
		    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
		    wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario)
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    pick.dump(imp_mod,filemodel);
	    filemodel.close();
	    
	    # write Ascii files with each component
	    write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
		dire=root_result+'/')
	    
	    # plot and compare with zbase
	    #if subscan[iscenario]<=1: compare_imp_vs_zbase(imp_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[subscan[iscenario]],
	    #	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

	    if (wake_calc):
		# write Ascii files with each component
		write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
		    dire=root_result+'/')
	    
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# dip only
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)		
		# plot and compare with zbase
		if subscan[iscenario]<=1: compare_imp_vs_zbase(wake_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[subscan[iscenario]],
	    	    save=root_result+'/plot_wake_vs_zbase_'+machine+scenario,
		    xlim=[1e-5,1e6],ylim=[1e8,1e19],wake_flag=True);

    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios
	plot_compare_imp_model(imp_mod_list[:2],legscen[subscan[:2]],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios_v1"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios_v1"+strsubscan,
	    ylim=[1e5,1e10],bounds=[8e3,2e9],legpos=3);
	
	plot_compare_imp_model(imp_mod_list[2:],legscen[subscan[2:]],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios_v2"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios_v2"+strsubscan,
	    ylim=[1e5,1e10],bounds=[8e3,2e9],legpos=3);
	
	#maxratio_mb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	#    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	#    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	#    ylim=[1e5,1e10],bounds=[8e3,40e6],legpos=3);

	# plot and compare all scenarios with 2012 wake
	if wake_calc:
	    maxratio_w=plot_compare_imp_model(wake_mod_list,legscen[subscan],listcomp=['Wlong','Wxdip','Wydip'],
		saveimp=root_result+'/plot_wake_'+machine+"_scenarios"+strsubscan,
		saveratio=root_result+'/plot_wake_ratio_'+machine+"_scenarios"+strsubscan,
		xlim=[1e-5,1e6],ylim=[1e8,1e19],yliml=[1e6,1e19],bounds=[8e3,5e10],legpos=0);

        
    if False:
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
	    	Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],	omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-4,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV',
		queue='2nd',dire=root_result+'/',flagQpscan_outside=True);
	    

	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			# plots vs Q'
			for Nb in Nbscanplot:

			    # initialize plots vs Q'
			    figQp=[];axQp=[];figQpm0=[];axQpm0=[];
			    for ir in range(2):
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			    	fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

			    # output files name for plots vs Q'
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

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
				    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
				    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");

				    if (M==1):
					# compare with HEADTAIL
				    	EstrHEAD=float_to_str(round(Escan[subscan[iscenario]]/1e12))+'TeV';
					fact=1;
					if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0 (conversion growth rate to tune shift)
					
					if (scenario=='_2012')and(Nb in NbscanHEADTAILv1):
					    nsl=500;npr=1000000;nlin=1; # HEADTAIL parameters for comparison
					    if (ir==0): fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_damper_dip_1b_ntwake20_nkick1_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_drate"+float_to_str(damp)+"_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin);
					    sufHEADTAIL="_aver_Sussix_most_tau_finer.txt";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
					
					elif ((scenario=='_2012_v2')or(scenario=='_HLLHC_round_Crab_wire_TCT'))and(Nb in NbscanHEADTAILv2):
					    nsl=200;npr=500000;nlin=1; # HEADTAIL parameters for comparison
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenario+"_dip_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_pre1_nsig2_drate"+float_to_str(damp);
					    sufHEADTAIL="_aver_most_tau_finer.txt";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");


			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
			        axQp[ir].set_xlim([-15,25]);
				end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
				end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))

			
			# plots vs Nb, and TMCI plots
			for Qp in Qpscanplot:

			    # initialize plots vs Nb
			    figNb=[];axNb=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figNb.append(fig);axNb.append(ax);
			    # initialize TMCI plots
			    figTMCI=[];axTMCI=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figTMCI.append(fig);axTMCI.append(ax);
			    
			    pat=['.',''];

			    # output file name for plots vs Q'
			    fileoutplotNb=root_result+'/plot_vs_Nb_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];ylim=([-5,3],[-0.01,1]);
			    for ir,r in enumerate(['real','imag']):

				for iscenario,scenario in enumerate(scenarioscan[subscan]):

				    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
				    # output file name for data vs Nb
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataNb=root_result+'/data_vs_Nb_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

				    ts_most=getattr(tuneshiftQp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,0],r);
				    data=np.hstack((Nbscan.reshape((-1,1)),ts_most.reshape((-1,1))));
				    write_ncol_file(fileoutdataNb+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_most_unstable")

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Nbscan[~np.isnan(np.squeeze(ts_most))],np.squeeze(sgn*ts_most)[~np.isnan(np.squeeze(ts_most))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axNb[ir],0,xlab="Intensity [p+/bunch]");
				    
				    # set plot axes
				    axNb[ir].set_xlim([0,8e11]);
				    maxy=np.ceil(np.max(np.abs(ts_most))*1e5)/1e5;
				    if ir==0: axNb[ir].set_ylim([-maxy,0]);
				    else: ylim[ir][1]=np.ceil(maxy/Qs*5)/5.;axNb[ir].set_ylim([0,maxy]);

				    if ((scenario=='_2012_v2')or(scenario=='_HLLHC_round_Crab_wire_TCT'))and(M==1):
					# compare with HEADTAIL
				    	EstrHEAD=float_to_str(round(Escan[subscan[iscenario]]/1e12))+'TeV';
					fact=1;
					if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0 (conversion growth rate to tune shift)
					sufHEADTAIL="_aver_most_tau_finer.txt";
					
					if (damp==0)and(Qp==0):
					    # HEADTAIL parameters for comparison
					    nsl=200;npr=500000;
					    nlin=1;nlinstr='_lin';dipstr='_dip';
					    #nlin=4;nlinstr='_nlin';dipstr='';
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotNb=fileoutplotNb+'_vs_HEADTAIL'+nlinstr+dipstr;
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenario+dipstr+"_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_csi"+float_to_str(Qp)+"_qsec0_oct0_drate"+float_to_str(damp)+"_pre1_nsig2";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
                                	    plot(s[:,0]*1e11,fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axNb[ir],0,xlab="Intensity [p+/bunch]");
					
					elif (damp==0.02)and(Qp==15):
					    # HEADTAIL parameters for comparison
					    nsl=200;npr=500000;
					    nlin=1;nlinstr='_lin';dipstr='_dip';
					    #nlin=4;nlinstr='_nlin';dipstr='';
					    if (ir==0)and(scenario=='_2012_v2'): fileoutplotNb=fileoutplotNb+'_vs_HEADTAIL'+nlinstr+dipstr;
					    rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_HLLHC_v2/fcutoff_50GHz/LHC_1b_"+EstrHEAD+"_B"+beam+scenario+dipstr+"_ntwake10_nkick1_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_csi"+float_to_str(Qp)+"_qsec0_oct0_drate"+float_to_str(damp)+"_pre1_nsig2";
					    # plot
					    s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
					    plot(s[:,0]*1e11,fact*s[:,3*iplane+ir+1],'HEADTAIL, '+legscen[subscan[iscenario]]+' '+Estr,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axNb[ir],0,xlab="Intensity [p+/bunch]");
				    
				    # TMCI plot
				    ts=tuneshiftQp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,:];
				    plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,patcol=pat[ir]+col[iscenario],xlab='Nb [p+/b]',
					title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);
	    

			    # finish plots vs Nb and TMCI plots
			    for ir,r in enumerate(['real','imag']):
				end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
				end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r));


	    # plots stabilizing emittance vs Nb for certain Qp

	    #if flagcompute:
	    # test that first scenario is 2012 updated model, and that any other is an HL-LHC one
	    if (len(subscan)>1)and((subscan[0]==1)and(all(subscan[1:]>1))):
	    

		# parameters
		oct2012signscan=['_pos_oct','_neg_oct'];
		oct2012scan=[510.,250];Nb2012=1.5e11;emit2012=2.5;en2012=4.;
		octLS1=550.;enLS1=6.5; 
		octHL=550.*2.;enHL=7.; # 2 factor due to highest oct. efficiency in HL-LHC (from their beta) (rough estimate from 15cm round optics)
		#impfact=[0.77, 1.13, 1.53];
		# post-LS1 scenarios
		#legbeam=['Classical 25 ns','BCMS 25 ns', 'Classical 50 ns','BCMS 50 ns'];
		#emitbeam=[3.75,1.9,2.2,1.6];
		#intbeam=[1.35,1.15,1.65,1.6];
		#colbeam=['or','+r','or','+g'];
		# HL-LHC scenarios
		legbeam=['PIC low emit.','PIC high emit.', 'US1',"US2 low int. & low emit.","US2 high int. & high emit."];
		emitbeam=[1.8,2.22,2.62,2.26,2.5];
		intbeam=[1.38,1.38,1.9,1.9,2.2];
		colbeam=['or','+m','xg','dk','vb'];

		colscen=['b','m','k','c'];

		for isign,octsign in enumerate(oct2012signscan):

		    for iM,M in enumerate(Mscan):

			for idamp,damp in enumerate(dampscan):

			    # initialize plot
			    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);

			    # output file name for plots vs emittance
			    fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+strsubscan+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

			    # plot postLS1 beam parameters
			    for ib,legb in enumerate(legbeam):
	  			plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");

			    # with assumption of linearity between intensity and growth rate:
			    # oct*emit/en = impfact * factlin *int;
			    factor2012lin=oct2012scan[isign]*emit2012/(en2012*Nb2012/1e11);
			    # better : oct*emit/en^2 =  fact *imag(-tuneshift(int));
			    #factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
			    # take iM=0 (50ns case in 2012)
			    imtu0=np.min(getattr(tuneshiftQp[0,:,0,iQpaver,idamp,pylab.mlab.find(Nbscan==Nb2012),0,0,0],'imag'));
			    #imtu0=np.max(np.abs(tuneshiftQp[0,:,0,iQpaver,idamp,pylab.mlab.find(Nbscan==Nb2012),0,0,0]));
			    print imtu0
			    factor2012=oct2012scan[isign]*emit2012/(en2012**2*abs(imtu0));
			    # emit=fact*en^2/oct *imag(-tuneshift)

			    for iscenario,scenario in enumerate(scenarioscan[subscan[1:]]):

				# output file name for data int. vs emitttance
				Estr=float_to_str(round(Escan[subscan[iscenario+1]]/1e9))+'GeV';
				fileoutdataemit=root_result+'/data_int_vs_emit_'+machine+'_'+Estr+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

				imtu=np.zeros(len(Nbscan));
				for iNb,Nb in enumerate(Nbscan):
		    		    imtu[iNb]=np.min(getattr(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0],'imag'));
		    		    #imtu[iNb]=np.max(np.abs(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0]));
				print imtu

				emitLS1=factor2012*enHL**2/octHL * np.abs(imtu);
				plot(emitLS1,Nbscan/1e11,"Limit "+legscen[subscan[iscenario+1]]+' '+Estr,colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
				# with linear assumption
				if damp==0: emitLS1lin=np.max(np.abs(maxratio_mb[:,iscenario,0]+0*maxratio_mb[:,iscenario,1]))*factor2012lin*enHL/octHL * (Nbscan/1e11);
				else: emitLS1lin=np.max(np.abs(maxratio_sb[:,iscenario,0]+0*maxratio_sb[:,iscenario,1]))*factor2012lin*enHL/octHL * (Nbscan/1e11);
				plot(emitLS1lin,Nbscan/1e11,"Limit with "+legscen[subscan[iscenario+1]]+' '+Estr+', linear approx.','--'+colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
				ax.set_xlim([0,6]);ax.set_ylim([0,6]);

				data=np.hstack((emitLS1.reshape((-1,1)),Nbscan.reshape((-1,1))));
				write_ncol_file(fileoutdataemit+'_'+r+'.dat',data,header="emit\t"+strpart[ir]+"Nb")


			    # finish plot
			    end_figure(fig,ax,save=flagsave*fileoutplotemit,legpos=1,fontsize=35)


    if not(flagsave): pylab.show();
