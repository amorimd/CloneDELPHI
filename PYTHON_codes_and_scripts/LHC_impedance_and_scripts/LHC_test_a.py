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
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/HLLHC_test_a';
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
    subscan=np.array([0]);strsubscan='_HLLHC_round_Crab_BBC';margin_factor=1;
    scenario='_HLLHC_round_Crab_wire_TCT';iscenario=0;
    legscen=np.array(['']);
    squeezescan=np.array(['0p15m_round']);
    Escan=np.array([7e12]);
    optionCrabscan=np.array(['']);
    optionBBCscan=np.array([2]);
    # scan of a parameter
    ascan=np.array([5.,6.,7.,8.,9.,10.,12.]);
    
    param_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat']);
    settings_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat']);
    
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_Allthemachine_4TeV_B1_physics_fill_3265.dat','_Allthemachine_4TeV_B1_physics_fill_3265.dat'])
    Qpscan=np.arange(-20,31,2);
    #Qpaver=np.array([14,15,16]);
    Qpscanplot=np.array([0,16]);
    #iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver

    dampscan=np.array([0,0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.7e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    Mscan=np.array([1]); # scan on number of bunches


    tuneshiftQp=np.zeros((len(ascan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(ascan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    # compute imp. model
    collscen='_HLLHC';
    imp_mod,wake_mod=HLLHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll_scan[subscan[iscenario]],
	    settings_filename_coll_scan[subscan[iscenario]],dire=path_here+"LHC_elements/",
	    commentcoll=collscen,direcoll='Coll'+collscen+'_v2/',
	    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
	    wake_calc=wake_calc,optionCrab=optionCrabscan[subscan[iscenario]],margin_factor=margin_factor,
	    optionBBC=optionBBCscan[subscan[iscenario]],flagplot=flagplot,root_result=root_result,commentsave=scenario)

    imp_mod_list.append(imp_mod);
    wake_mod_list.append(wake_mod);


    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

        # DELPHI now
	    
	for ia,a in enumerate(ascan):
	    
	    # DELPHI computation
	    tuneshiftQp[ia,:,:,:,:,:,:,:,:],tuneshiftm0Qp[ia,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_list[iscenario],
	    	Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],	omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-3,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV'+'_a'+float_to_str(a),
		queue='1nd',dire=root_result+'/',flagQpscan_outside=True);
	    

	# now the plots (outside loop on a)
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
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];
			    for ir,r in enumerate(['real','imag']):

				for ia,a in enumerate(ascan):

				    # output files name for data vs Qp
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+scenario+'_a'+float_to_str(a)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+scenario+'_a'+float_to_str(a)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[ia,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    
				    tsm0=getattr(tuneshiftm0Qp[ia,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshiftm0")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+', a='+str(a),col[ia],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
				    plot(Qpscan,np.squeeze(sgn*tsm0),'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+', a='+str(a),col[ia],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ (mode 0)",axQpm0[ir],0,xlab=" $ Q^' $ ");

			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
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
			    fileoutplotNb=root_result+'/plot_vs_Nb_'+machine+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];ylim=([-5,3],[-0.01,1]);
			    for ir,r in enumerate(['real','imag']):

				for ia,a in enumerate(ascan):

				    # output file name for data vs Nb
				    Estr=float_to_str(round(Escan[subscan[iscenario]]/1e9))+'GeV';
				    fileoutdataNb=root_result+'/data_vs_Nb_'+machine+'_'+Estr+scenario+'_a'+float_to_str(a)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;

				    ts_most=getattr(tuneshiftQp[ia,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,0],r);
				    data=np.hstack((Nbscan.reshape((-1,1)),ts_most.reshape((-1,1))));
				    write_ncol_file(fileoutdataNb+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_most_unstable")

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Nbscan,np.squeeze(sgn*ts_most),'DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+', a='+str(a),col[ia],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axNb[ir],0,xlab="Intensity [p+/bunch]");
				    
				    # set plot axes
				    axNb[ir].set_xlim([0,8e11]);
				    maxy=np.ceil(np.max(np.abs(ts_most))*1e5)/1e5;
				    if ir==0: axNb[ir].set_ylim([-maxy,0]);
				    else: axNb[ir].set_ylim([0,maxy]);ylim[ir][1]=np.ceil(maxy/Qs*5)/5.;

				    # TMCI plot
				    ts=tuneshiftQp[ia,iplane,iM,pylab.mlab.find(Qpscan==Qp),idamp,:,0,0,:];
				    plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='DELPHI, '+legscen[subscan[iscenario]]+' '+Estr+', a='+str(a),patcol=pat[ir]+col[ia],xlab='Nb [p+/b]',
					title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);
	    

			    # finish plots vs Nb and TMCI plots
			    for ir,r in enumerate(['real','imag']):
				end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
				end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r),legpos=3);


    if not(flagsave): pylab.show();
