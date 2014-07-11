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
from HLLHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/LHC_BS_TMCI';
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=True; # to plot impedances
    nevery=1; # downsampling of the impedance (take less points than in the full model)

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition
    subscan=np.array([0,1,2]);strsubscan='_LHC_inj_4TeV_7TeV';
    #subscan=np.array([0]);strsubscan='_LHC_inj';
    scenarioscan=np.array(['_BS_inj','_BS_flat_top_2012','_BS_nominal']);
    legscen=np.array(['LHC BS RW injection','LHC BS RW flat top 2012','LHC BS RW nominal']);
    squeezescan=np.array(['11m_10m_11m_10m','0p6m_3m_0p6m_3m','0p55m_10m_0p55m_10m']);
    Escan=np.array([450e9,4e12,7e12]);
    
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_BeamScreens_450GeV_B1.dat','_BeamScreens_4TeV_B1.dat'])
    Qpscan=np.array([0]);
    Qpscanplot=np.array([0]);

    dampscan=np.array([0]); # damper gain scan
    Nbscan=np.arange(1.e11,100.2e11,2.e10); # intensity scan
    #Mscan=np.array([1,1782,3564]); # scan on number of bunches
    Mscan=np.array([1]); # scan on number of bunches

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
	
	# compute imp. model (RW impedance of BS only, without weld)
	
	# beta functions
	beta_filename_rest="../LHC_elements/LHC_beta_length_B"+str(beam)+"_sq"+squeezescan[subscan[iscenario]]+".dat"
	# parameters
	param_filename_RW="../LHC_elements/LHC_RW_param.dat"
	# names of the BS
	namesRW=read_ncol_file_identify_header(param_filename_RW,'name');
	namesBS=select_LHC_names(namesRW,pattern='BS'); # only BS
	# imp. model
	imp_mod_RW_BS,wake_mod_RW_BS=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	    beta_filename_rest,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=namesBS,
	    lxplusbatch=lxplusbatchImp,comment='_'+Estr,dire='BS_v2_'+Estr+'/');

	imp_mod_list.append(imp_mod_RW_BS);
	wake_mod_list.append(wake_mod_RW_BS);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    pick.dump(imp_mod_RW_BS,filemodel);
	    filemodel.close();
	    
	    # write Ascii files with each component
	    write_imp_wake_mod(imp_mod_RW_BS,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
		dire=root_result+'/')
	    
	    # plot and compare with zbase
	    if subscan[iscenario]<=1: compare_imp_vs_zbase(imp_mod_RW_BS,root_zbase=zbaseroot,suffix_zbase=zbasesuf[subscan[iscenario]],
	    	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario,ncomp=3);
    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios with 2012 impedance
	maxratio_sb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	    ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3);
	
        # DELPHI now
	for iscenario,scenario in enumerate(scenarioscan[subscan]):

	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
	    
	    # DELPHI computation
	    tuneshiftQp[iscenario,:,:,:,:,:,:,:,:],tuneshiftm0Qp[iscenario,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_list[iscenario],
	    	Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],	omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-3,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV',
		queue='8nh',dire=root_result+'/',flagQpscan_outside=True);
	    

	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

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

			    strpart=['Re','Im'];ylim=([-5,3],[-0.001,0.02]);
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
				    #axNb[ir].set_xlim([0,8e11]);
				    maxy=np.ceil(np.max(np.abs(ts_most))*1e5)/1e5;
				    if ir==0: axNb[ir].set_ylim([-maxy,0]);
				    else: axNb[ir].set_ylim([0,maxy]);

				    # TMCI plot
				    ts=tuneshiftQp[iscenario,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,:];
				    plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,patcol=pat[ir]+col[iscenario],xlab='Nb [p+/b]',
					title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);
	    

			    # finish plots vs Nb and TMCI plots
			    for ir,r in enumerate(['real','imag']):
				end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
				end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r));



    if not(flagsave): pylab.show();
