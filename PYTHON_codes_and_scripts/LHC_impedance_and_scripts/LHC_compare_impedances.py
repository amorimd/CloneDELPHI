#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

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
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_compare';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=False; # to plot impedances
    nevery=2; # downsampling of the impedance (take less points than in the full model)
    if lxplusbatchDEL.startswith('retrieve'): nevery=200;

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','xr','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # cases scan definition
    #subscan=np.array([0,1]);strsubscan='_LHC_v2_2011_100turns_sq1m_vs_2012_50turns';margin_factor=1;
    #subscan=np.array([0]);strsubscan='_LHC_v2_2011';margin_factor=1;
    #subscan=np.array([2,3]);strsubscan='_LHC_v2_2012_inj_vs_2015_inj';margin_factor=1;
    subscan=np.array([4,5]);strsubscan='_LHC_v2_2011_0turns_sq1m_vs_2012_0turns';margin_factor=1;
    scenarioscan=np.array(['_2011','_2012_v2','_2012_v2','_2012_v2','_2011','_2012_v2']);
    commentscenario_scan=np.array(['','','','_2015_optics_change','_no_damper','_no_damper']);
    legscen=np.array(['3.5 TeV 2011 (17/10)','4TeV 2012','Injection 2012','Injection 2015','3.5 TeV 2011 (17/10) (no damper)','4TeV 2012 (no damper)']);
    beamscan=np.array(['1','1','1','1','1','1']);
    Nbscan=np.array([145., 150., 130., 130., 145., 150.])*1e9; # intensity scan
    squeezescan=np.array(['1p0m_10m_1p0m_3m','0p6m_3m_0p6m_3m','11m_10m_11m_10m','11m_10m_11m_10m_2015','1p0m_10m_1p0m_3m','0p6m_3m_0p6m_3m']);
    Escan=np.array([3.5e12, 4e12, 450e9, 450e9, 3.5e12, 4e12]);
    dampscan=np.array([0.01, 0.02, 0.02, 0.02, 0., 0.]); # damper gain scan
    
    #zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    #zbasesuf=np.array(['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_4TeV_B1_physics_fill_3265.dat'])
    # other scans
    Qpscan=np.arange(-10,21);
    #Qpaver=np.array([14,15,16]);
    Qpscanplot=np.array([0,15]);
    #iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver

    #Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    #Nbscanplot=np.array([1.5e11,1.7e11,2.5e11,3.e11,4.e11]); # intensity scan for plot vs Qp
    #Mscan=np.array([1,1782,3564]); # scan on number of bunches
    Mscan=np.array([1782]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),1,1,1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),1,1,1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    for iscenario,scenario in enumerate(scenarioscan[subscan]):
    	
	beam=beamscan[subscan[iscenario]];
	sq=squeezescan[subscan[iscenario]];
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	commentscenario=commentscenario_scan[subscan[iscenario]];
	print "scenario: ",scenario,commentscenario;
	
	if scenario.startswith('_2012')and(E==4e12):
	    settings_filename_coll=path_here+'Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';
	    sufcoll='_2012.txt';
	    param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_'+float_to_str(round(E/1e9))+'GeV_sq'+sq[:3]+'_b'+beam+sufcoll;
	elif scenario.startswith('_2012')and(E==450e9):
	    settings_filename_coll=path_here+'Coll_settings/coll_ph1_beta_'+float_to_str(round(E/1e9))+'GeV_b'+beam+'_2012.txt';
	    param_filename_coll=settings_filename_coll;
	elif scenario.startswith('_2011'):
	    settings_filename_coll=path_here+'Coll_settings/coll_settings_B'+beam+'_'+float_to_str(round(E/1e9))+'GeV_fill1727.txt';
	    sufcoll='.txt';
	    sq=sq[:2]+'5'+sq[3:]; # very dirty (not same optics for collimators - 1.5m squeeze instead of 1m ... but should have little impact)
	    param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_'+float_to_str(round(E/1e9))+'GeV_sq'+sq[:3]+'_b'+beam+sufcoll;
	elif scenario.startswith('_2012'):
	    settings_filename_coll=path_here+'Coll_settings/coll_ph1_beta_'+float_to_str(round(E/1e9))+'GeV_b'+beam+'_2012.txt',
	    sufcoll='_2012.txt';
	    param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_'+float_to_str(round(E/1e9))+'GeV_sq'+sq[:3]+'_b'+beam+sufcoll;
	
	print settings_filename_coll,param_filename_coll
	
	# compute imp. model
	imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
		settings_filename_coll,dire=path_here+"LHC_elements/",commentcoll=scenario,direcoll='Coll_2012_v2/',
		lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
		wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario+commentscenario)
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    #filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    #pick.dump(imp_mod,filemodel);
	    #filemodel.close();
	    
	    # write Ascii files with each component
	    #write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    #	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
	    #	dire=root_result+'/')
	    
	    # plot and compare with zbase
	    #if subscan[iscenario]<=1: compare_imp_vs_zbase(imp_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[subscan[iscenario]],
	    #	save=root_result+'/plot_imp_vs_zbase_'+machine+scenario);

	    if (wake_calc):
		# write Ascii files with each component
		#write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	#    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
		#    dire=root_result+'/')
	    
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+commentscenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# dip only
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+commentscenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)		

    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios
	maxratio_sb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	    ylim=[1e5,1e9],bounds=[40e6,2e9],legpos=3,markimp=['','x','','','']);
	
	
	#maxratio_mb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
	#    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	#    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	#    ylim=[1e5,1e10],bounds=[8e3,40e6],legpos=3);

	# plot and compare all scenarios with 2012 wake
	if wake_calc:
	    maxratio_w=plot_compare_imp_model(wake_mod_list,legscen[subscan],listcomp=['Wlong','Wxdip','Wydip'],
		saveimp=root_result+'/plot_wake_'+machine+"_scenarios"+strsubscan,
		saveratio=root_result+'/plot_wake_ratio_'+machine+"_scenarios"+strsubscan,
		xlim=[1e-5,1e6],ylim=[1e8,1e19],yliml=[1e6,1e19],bounds=[8e3,5e10],legpos=0,markimp=['','x','','','']);

        
	# DELPHI scans now
	for iscenario,scenario in enumerate(scenarioscan[subscan]):

	    Nb=Nbscan[subscan[iscenario]];
	    damp=dampscan[subscan[iscenario]];
	    commentscenario=commentscenario_scan[subscan[iscenario]];
	    
	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	    print scenario,", Qs=",Qs,", omegas=",omegas;
	    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
	    
	    # DELPHI computation
	    tuneshiftQp[iscenario,:,:,:,:,:,:,:,:],tuneshiftm0Qp[iscenario,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_list[iscenario],
	    	Mscan,Qpscan,[damp],[Nb],[omegas],[dphase],omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-4,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+commentscenario+'_'+float_to_str(round(E/1e9))+'GeV',
		queue='2nd',dire=root_result+'/',flagQpscan_outside=True);
	    

	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    # plots vs Q'

		    # initialize plots vs Q'
		    figQp=[];axQp=[];figQpm0=[];axQpm0=[];
		    for ir in range(2):
			fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

		    # output files name for plots vs Q'
		    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+strsubscan+'_'+str(M)+'b_converged'+strnorm[flagnorm]+'_'+plane;
		    fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+strsubscan+'_'+str(M)+'b_converged'+strnorm[flagnorm]+'_'+plane;

		    strpart=['Re','Im'];
		    for ir,r in enumerate(['real','imag']):

			for iscenario,scenario in enumerate(scenarioscan[subscan]):

			    # output files name for data vs Qp
			    Nb=Nbscan[subscan[iscenario]];
			    damp=dampscan[subscan[iscenario]];
			    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
			    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+commentscenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+commentscenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    #if flagcompute:
			    ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,0,0,0,0,0],r);
			    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
			    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")

			    tsm0=getattr(tuneshiftm0Qp[iscenario,iplane,iM,:,0,0,0,0],r);
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
			    plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))],'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr,col[iscenario],ylab + " (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");



		    # finish plot vs Qp
		    for ir,r in enumerate(['real','imag']):
			#axQp[ir].set_xlim([-10,20]);
			end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
			end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))

						

    if not(flagsave): pylab.show();
