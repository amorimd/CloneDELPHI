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
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from tables_lib import select_in_table
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

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/postLS1';
    #root_result='/home/nmounet/Documents/DELPHI_results/'+machine+'/postLS1';
   
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition
    scenarioscan=['_2012','_relaxed','_mm_kept','_nominal','_sigma_kept'];
    legscen=['2012','relaxed settings','tight settings (mm kept)','nominal settings','tight settings ($ \sigma $ kept)'];
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_7TeV_B1_postLS1_veryrelaxed.dat',
    	'_Allthemachine_7TeV_B1_postLS1_relaxed.dat','_Allthemachine_7TeV_B1_postLS1_baseline.dat','_Allthemachine_7TeV_B1_postLS1_baseline.dat']
    Qpscan=np.arange(-20,31,1);
    Qpaver=np.array([14,15,16]);
    iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver

    dampscan=np.array([0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.7e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # compute non changing part of the impedance
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    # rest (i.e. not coll.) of the machine wall impedance
    param_filename_rest="../LHC_elements/beam_screens_warm_pipe_LHC_param.dat"
    #beta_filename_rest="../LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p6m_3m_0p6m_3m.dat"
    #beta_filename_rest="../LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat"
    
    #imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
    #	wake_calc=False,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatch);
    
    # broad-band model
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=False,
    	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());
    

    tuneshiftQp=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    for iscenario,scenario in enumerate(scenarioscan):
    
	# scenario is between 0 and 4:
	# 0 -> 4 TeV 2012 tight settings,
	# 1 -> 6.5 TeV relaxed settings,
	# 2 -> 6.5 TeV same settings as 2012 in mm,
	# 3 -> 6.5 TeV nominal settings (design report),
	# 4 -> 6.5 TeV same settings as 2012 in sigma,
	# settings for cases 1 to 4 are from Roderik Bruce.

	# postLS1 files
	if iscenario==0:
	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=4e12);
	    param_filename_coll="../Coll_settings/coll_ph1_beta_"+str(int(E/1e9))+"GeV_sq0p6_b1_2012.txt";
	    beta_filename_coll=param_filename_coll;
	    settings_filename_coll="../Coll_settings/coll_settings_physics_fill_3265_B1.txt";
	else:
    	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    	    param_filename_coll="../Coll_settings/collgaps_fromRoderik_modifNico_materialnames"+scenario+".dat";
            beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;
	
	# compute model for collimators
	imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,
	    comment=scenario,dire='Coll'+scenario+'/');

	# compute the rest (once for 4 TeV and once for 6.5TeV - depends on energy)
	if iscenario==0: beta_filename_rest="../LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p6m_3m_0p6m_3m.dat"
	else: beta_filename_rest="../LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat"
    	
	if (iscenario<2): imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
		wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,comment='_'+Estr,dire='Rest_'+Estr+'/');
    
	# add up
	imp_mod=[];wake_mod=[];
	add_impedance_wake(imp_mod,imp_mod_coll,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll,1,1);
	add_impedance_wake(imp_mod,imp_mod_rest,1,1);
	add_impedance_wake(wake_mod,wake_mod_rest,1,1);
	add_impedance_wake(imp_mod,imp_mod_BB,1,1);
	add_impedance_wake(wake_mod,wake_mod_BB,1,1);
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    pick.dump(imp_mod,filemodel);
	    filemodel.close();
	    
	    # plot and compare with zbase
	    compare_imp_vs_zbase(imp_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
	    	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

	    if (wake_calc):
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"wakeforhdtl_PyZbase+Allthemachine_"+Estr+"_B"+str(beam)+"_"+scenario,beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# plot and compare with zbase
		compare_wake_vs_zbase(wake_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
	    	    save=root_result+'/plot_wake_vs_zbase_'+machine+scenario,
		    xlim=[1e-5,1e6],ylim=[1e8,1e19],wake_flag=True);

    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios with 2012 impedance
	maxratio_sb=plot_compare_imp_model(imp_mod_list,legscen,listcomp=['Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"postLS1_scenarios",
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_scenarios",
	    ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=1);

	maxratio_mb=plot_compare_imp_model(imp_mod_list,legscen,listcomp=['Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"postLS1_scenarios",
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_scenarios",
	    ylim=[1e5,1e10],bounds=[8e3,40e6],legpos=1);

	# plot and compare all scenarios with 2012 wake
	if wake_calc:
	    maxratio_w=plot_compare_imp_model(wake_mod_list,legscen,listcomp=['Wxdip','Wydip'],
		saveimp=root_result+'/plot_wake_'+machine+"postLS1_scenarios",
		saveratio=root_result+'/plot_wake_ratio_'+machine+"postLS1_scenarios",
		xlim=[1e-5,1e6],ylim=[1e8,1e19]);

        # DELPHI loops now
	for iscenario,scenario in enumerate(scenarioscan):

	    if iscenario==0:
		machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=4e12);
	    else:
    		machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    	    
	    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

	    for iplane,plane in enumerate(['x','y']):
	        # select Zxdip or Zydip
		for iw in imp_mod_list[iscenario]:
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
		    	tuneshiftQp[iscenario,iplane,iM,iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
		    		nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
		    		a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
		    		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
		    		kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
		    		lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
		    		queue='1nd',dire=root_result+'/');
		    #print "Waiting 20 minutes...";
		    #time.sleep(1200);
		    
		    # case with Qpscan inside each lxplus job
		    #tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
		    #
		    #tuneshiftQp[iscenario,iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
		#		nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
		#		a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
		#		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
		#		kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
		#		lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
		#		queue='1nw',dire=root_result+'/');


	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    # plots vs Q'
	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			for Nb in Nbscanplot:

			    # initialize plots vs Qp
			    figQp=[];axQp=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

			    # output file name for plots vs Qp
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];
			    for ir,r in enumerate(['real','imag']):

				for iscenario,scenario in enumerate(scenarioscan):

				    # output file name for data vs Qp
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

				    if (M==1)and(scenario=='_2012'):
					# compare with HEADTAIL
					nsl=500;npr=1000000;nlin=1; # HEADTAIL parameters for comparison
					fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
					rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_damper_dip_1b_ntwake20_nkick1_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_drate"+float_to_str(damp)+"_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin);
					#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
					#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
					sufHEADTAIL="_aver_Sussix_most_tau_finer.txt";
					s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
					fact=1;
                                	if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0
                                	plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
				if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
				else: end_figure(figQp[ir],axQp[ir]);


	    # plots stabilizing emittance vs Nb for certain Qp

	    #if flagcompute:

	    # parameters
	    oct2012signscan=['_pos_oct','_neg_oct'];
	    oct2012scan=[510.,250];Nb2012=1.5e11;emit2012=2.5;en2012=4.;
	    octLS1=550;enLS1=7;
	    #impfact=[0.77, 1.13, 1.53];
	    legbeam=['classical 25 ns','BCMS 25 ns', 'classical 50 ns','BCMS 50 ns'];
	    emitbeam=[3.75,1.9,2.2,1.6];
	    intbeam=[1.35,1.15,1.65,1.6];
	    colbeam=['or','+r','og','+g'];

	    colscen=['b','m','k','c'];

	    for isign,octsign in enumerate(oct2012signscan):
	    
		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			# initialize plot
			fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);

			# output file name for plots vs emittance
			fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

			# plot postLS1 beam parameters
			for ib,legb in enumerate(legbeam):
	  		    plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");

			# with assumption of linearity between intensity and growth rate:
			# oct*emit/en = impfact * factlin *int;
			factor2012lin=oct2012scan[isign]*emit2012/(en2012*Nb2012/1e11);
			# better : oct*emit/en^2 =  fact *imag(-tuneshift(int));
			#factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
			# take 50ns case in 2012
			imtu0=np.min(getattr(tuneshiftQp[0,:,pylab.mlab.find(Mscan==1782),iQpaver,idamp,pylab.mlab.find(Nbscan==Nb2012),0,0,0],'imag'));
			print imtu0
			factor2012=oct2012scan[isign]*emit2012/(-en2012**2*imtu0);
			# emit=fact*en^2/oct *imag(-tuneshift)

			for iscenario,scenario in enumerate(scenarioscan[1:]):

			    # output file name for data int. vs emitttance
			    fileoutdataemit=root_result+'/data_int_vs_emit_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged'+strnorm[flagnorm];

			    imtu=np.zeros(len(Nbscan));
			    for iNb,Nb in enumerate(Nbscan):
		    		imtu[iNb]=np.min(getattr(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0],'imag'));
			    print imtu

			    emitLS1=factor2012*enLS1**2/octLS1 * (-imtu);
			    plot(emitLS1,Nbscan/1e11,"Limit with "+legscen[iscenario+1],colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
			    ax.set_xlim([0,6]);
			    # with linear assumption
			    if damp==0: emitLS1lin=np.max(np.abs(maxratio_mb[:,iscenario,0]+0*maxratio_mb[:,iscenario,1]))*factor2012lin*enLS1/octLS1 * (Nbscan/1e11);
			    else: emitLS1lin=np.max(np.abs(maxratio_sb[:,iscenario,0]+0*maxratio_sb[:,iscenario,1]))*factor2012lin*enLS1/octLS1 * (Nbscan/1e11);
	 		    
			    #plot(emitLS1lin,Nbscan/1e11,"Limit with "+legscen[iscenario+1]+', linear','--'+colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");

			    data=np.hstack((emitLS1.reshape((-1,1)),Nbscan.reshape((-1,1))));
			    write_ncol_file(fileoutdataemit+'_'+r+'.dat',data,header="emit\t"+strpart[ir]+"Nb")


			# finish plot
			if flagsave: end_figure(fig,ax,save=fileoutplotemit)
			else: end_figure(fig,ax);


    if not(flagsave): pylab.show();
