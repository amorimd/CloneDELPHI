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
import time# as ti
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

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/LHC_new_postLS1_bunchlength_Qs_tests';
    #root_result='../../scratch0/LHC_new_postLS1';
    
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

    # scenario definition
    scenarioscan=np.array(['_2012_v2','_feb2014_relaxed','_feb2014_mm_kept','_feb2014_nominal',
    	'_feb2014_sigma_kept','_may2014_mm_kept_IR3_15_sig_TCL6_open','_may2014_sigma_kept_TCL6_open',
	'_may2014_sigma_kept_TCL6_open_TCSG_12_sigma']);
    dircollscan=np.array(['_2012_v2','_relaxed','_mm_kept','_nominal','_sigma_kept',
    	'_mm_kept','_sigma_kept','_sigma_kept'])
    legscen=np.array(["2012","relaxed settings","tight settings (mm kept)","nominal settings",
    	"tight settings ($ \sigma $ kept)","mm kept with TCL6 open & 15 $ \sigma $ TCS IR3",
	"sigma kept with TCL6 open","sigma kept with TCL6 open & 12 $ \sigma $ TCS IR7"]);
    squeezescan=np.array(['0p6m_3m_0p6m_3m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m',
    	'0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m','0p55m_10m_0p55m_10m']);
    Escan=np.array([4e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12,6.5e12]);
    # choice of the scenario
    iscenario=5;scenario=scenarioscan[iscenario];
    
    param_filename_coll_root='../Coll_settings/collgaps_fromRoderik_modifNico_materialnames';
   
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_4TeV_B1_physics_fill_3265.dat'])
    Qpscan=np.arange(-20,31);
    Qpscanplot=np.array([0,15]);

    dampscan=np.array([0, 0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,5.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.3e11]); # intensity scan for plot vs Qp
    Mscan=np.array([1,3564]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # bunchlength and Qs scans
    taub_scan=np.array([1.,1.25])*1e-9;
    Qs_scan=np.array([1.3,1.83,2.12])*1e-3;
    omegas_scan=Qs_scan*omega0;
    
    # parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[iscenario]);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    print "scenario: ",scenario

    # select coll. settings file	
    if scenario.startswith('_2012'):
	param_filename_coll='../Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
	settings_filename_coll='../Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';
    else:
	param_filename_coll=param_filename_coll_root+scenario+'.dat';
	settings_filename_coll=param_filename_coll;

    # compute imp. model
    scenariobis=scenario.replace('_feb2014','');
    imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
	    settings_filename_coll,dire="../LHC_elements/",commentcoll=scenariobis,
	    direcoll='Coll'+dircollscan[iscenario]+'/',
	    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[iscenario],
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

	# plot and compare with zbase
	#if iscenario<=1: compare_imp_vs_zbase(imp_mod,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
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

    
    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	tuneshiftQp=np.zeros((len(taub_scan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegas_scan),1,kmaxplot),dtype=complex);
	tuneshiftm0Qp=np.zeros((len(taub_scan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegas_scan),1),dtype=complex);
	
	# DELPHI scans now
	for itaub,taub in enumerate(taub_scan):

    	    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
	    
	    # DELPHI computation
	    tuneshiftQp[itaub,:,:,:,:,:,:,:,:],tuneshiftm0Qp[itaub,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod,
	    	Mscan,Qpscan,dampscan,Nbscan,omegas_scan,[dphase],omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-4,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+float_to_str(taub*1e9)+'ns',
		queue='2nd',dire=root_result+'/',flagQpscan_outside=True);

	
	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x','y']):

		for iM,M in enumerate(Mscan):

		    # plots vs Q', for each Qs
		    for idamp,damp in enumerate(dampscan):

			for iQs,Qs in enumerate(Qs_scan):
			
			    for Nb in Nbscanplot:

				# initialize plots vs Q'
				figQp=[];axQp=[];figQpm0=[];axQpm0=[];
				for ir in range(2):
			    	    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			    	    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

				# output files name for plots vs Q'
				Estr=float_to_str(round(Escan[iscenario]/1e9))+'GeV';
				fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				strpart=['Re','Im'];
				for ir,r in enumerate(['real','imag']):

				    for itaub,taub in enumerate(taub_scan):

					# output files name for data vs Qp
					fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
					fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

					#if flagcompute:
					ts=getattr(tuneshiftQp[itaub,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),iQs,0,0],r);
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")

					tsm0=getattr(tuneshiftm0Qp[itaub,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),iQs,0],r);
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
					#else:
					#    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
					#    Qpscan=s[:,0];ts=s[:,1];

					sgn=1;sgnstr='';
					if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
					if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate

					plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))],"total bunch length = "+str(taub*1e9)+" ns",col[itaub],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
					plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))],"total bunch length = "+str(taub*1e9)+" ns",col[itaub],ylab+" (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");


				# finish plot vs Qp
				for ir,r in enumerate(['real','imag']):
			            #axQp[ir].set_xlim([-15,25]);#axQp[ir].set_ylim([0,1.6e-5]);
				    end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
				    end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))

			
		    # plots vs Q', for each taub
		    for idamp,damp in enumerate(dampscan):

			for itaub,taub in enumerate(taub_scan):
			
			    for Nb in Nbscanplot:

				# initialize plots vs Q'
				figQp=[];axQp=[];figQpm0=[];axQpm0=[];
				for ir in range(2):
			    	    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);
			    	    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQpm0.append(fig);axQpm0.append(ax);

				# output files name for plots vs Q'
				Estr=float_to_str(round(Escan[iscenario]/1e9))+'GeV';
				fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+Estr+scenario+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				strpart=['Re','Im'];
				for ir,r in enumerate(['real','imag']):

				    for iQs,Qs in enumerate(Qs_scan):

					ts=getattr(tuneshiftQp[itaub,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),iQs,0,0],r);
					tsm0=getattr(tuneshiftm0Qp[itaub,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),iQs,0],r);

					sgn=1;sgnstr='';
					if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
					if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate

					plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*ts)[~np.isnan(np.squeeze(ts))]," $ Q_s= $ "+str(Qs),col[iQs],ylab,axQp[ir],0,xlab=" $ Q^' $ ");
					plot(Qpscan[~np.isnan(np.squeeze(ts))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts))]," $ Q_s= $ "+str(Qs),col[iQs],ylab+" (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");


				# finish plot vs Qp
				for ir,r in enumerate(['real','imag']):
			            #axQp[ir].set_xlim([-15,25]);#axQp[ir].set_ylim([0,1.6e-5]);
				    end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))
				    end_figure(figQpm0[ir],axQpm0[ir],save=flagsave*(fileoutplotQpm0+'_'+r))

			
		    # plots vs Nb, and TMCI plots
		    for idamp,damp in enumerate(dampscan):
		    
		    	for iQs,Qs in enumerate(Qs_scan):

		    	    for itaub,taub in enumerate(taub_scan):

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
				    Estr=float_to_str(round(Escan[iscenario]/1e9))+'GeV';
				    fileoutplotNb=root_result+'/plot_vs_Nb_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutplotNbm0=root_result+'/plot_vs_Nb_m0_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

				    strpart=['Re','Im'];ylim=([-5,3],[-0.01,1]);
				    for ir,r in enumerate(['real','imag']):

					# output file name for data vs Nb
					fileoutdataNb=root_result+'/data_vs_Nb_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
					fileoutdataNbm0=root_result+'/data_vs_Nb_m0_'+machine+'_'+Estr+scenario+'_Qs'+float_to_str(Qs)+'_'+float_to_str(taub*1e9)+'ns_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

					ts_most=getattr(tuneshiftQp[itaub,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,iQs,0,0],r);
					data=np.hstack((Nbscan.reshape((-1,1)),ts_most.reshape((-1,1))));
					#write_ncol_file(fileoutdataNb+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_most_unstable")

					tsm0=getattr(tuneshiftm0Qp[itaub,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,iQs,0],r);
					data=np.hstack((Nbscan.reshape((-1,1)),tsm0.reshape((-1,1))));
					#write_ncol_file(fileoutdataNbm0+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_m0")

					sgn=1;sgnstr='';
					#if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
					if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate

					plot(Nbscan[~np.isnan(np.squeeze(ts_most))],np.squeeze(sgn*ts_most)[~np.isnan(np.squeeze(ts_most))],'','b',ylab,axNb[ir],0,xlab="Intensity [p+/bunch]");
					plot(Nbscan[~np.isnan(np.squeeze(tsm0))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(tsm0))],'','b',ylab+" (mode 0)",axNbm0[ir],0,xlab="Intensity [p+/bunch]");

					# set plot axes
					axNb[ir].set_xlim([0,5e11]);
					maxy=np.ceil(np.max(np.abs(sgn*ts_most))*1e5)/1e5;
					if ir==0: axNb[ir].set_ylim([-maxy,0]);
					else: ylim[ir][1]=np.ceil(maxy/omegas*5)/5.;axNb[ir].set_ylim([0,maxy]);

					# TMCI plot
					ts=tuneshiftQp[itaub,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,iQs,0,:];
					plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='',patcol=pat[ir],xlab='Nb [p+/b]',
					    title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);


				    # finish plots vs Nb and TMCI plots
				    for ir,r in enumerate(['real','imag']):
					end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
					end_figure(figNbm0[ir],axNbm0[ir],save=flagsave*(fileoutplotNbm0+'_'+r))
					end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r));



    if not(flagsave): pylab.show();
