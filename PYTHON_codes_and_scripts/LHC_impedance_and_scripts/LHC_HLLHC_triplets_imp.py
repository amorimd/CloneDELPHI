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
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from string_lib import *
from tables_lib import select_in_table
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param
from LHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    # post-LS1 coll. model
    model='_relaxed';
    E=6.5e12; # energy in eV (4e12 -> 2012 case, 6.5e12 -> post-LS1 case)
    model='_2012';
    E=4e12; # energy in eV (4e12 -> 2012 case, 6.5e12 -> post-LS1 case)

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
    print machine+model.replace('_',' ')+": sigmaz=",sigmaz;

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/triplets_HLLHC';
    os.system("mkdir -p "+root_result);
    suffix='_triplets'; # suffix for impedance plots
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    flagDEL=0; # 0 to avoid DELPHI computing/plots
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':','.-'];
    

    # scan definition
    BBscan=['','_noBS98','_slope0p4'];
    holescan=['','_lowthickness'];
    RWscan=['_aC_Cu_ss','_aC_ss','_Cu_ss','_Cu_ss_noBS98','_ss'];
    #nmat=[3,2,2,2,2]; # number of material layers for each case
    
    Qpscan=np.arange(0,21,1);
    Qpaver=np.array([14,15,16]);
    iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver
    dampscan=np.array([0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,5.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.5e11]); # intensity scan for plot vs Qp
    Mscan=np.array([1]); # scan on number of bunches

    # parameters for power loss computation
    NbLHC=1.5e11;NbHLLHC=3.5e11;Mp=2*1404; # 2 50ns beams
    
    # compute non changing part of the impedance
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    if (model=='_2012'):
	# 2012 files
	param_filename_coll=path_here+"Coll_settings/coll_ph1_beta_"+str(int(E/1e9))+"GeV_sq0p6_b1_2012.txt";
	beta_filename_coll=param_filename_coll;
	settings_filename_coll=path_here+"Coll_settings/coll_settings_physics_fill_3265_B1.txt";
    else:
	# postLS1 files
    	param_filename_coll=path_here+"Coll_settings/collgaps_fromRoderik_modifNico_materialnames"+model+".dat";
        beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;

    # compute model for collimators
    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,
	comment=model,dire='Coll'+model+'/');

    # rest (i.e. not coll.) of the machine wall impedance
    param_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_param.dat"
    if E==4e12:
    	beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p6m_3m_0p6m_3m.dat"
    elif E>=6.5e12:
    	beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat"
    
    imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
    	wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,comment='_'+Estr,dire='Rest_'+Estr+'/');
    
    # broad-band model
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());
    
    # add up
    imp_mod=[];wake_mod=[];
    add_impedance_wake(imp_mod,imp_mod_coll,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll,1,1);
    add_impedance_wake(imp_mod,imp_mod_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_rest,1,1);
    add_impedance_wake(imp_mod,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod,wake_mod_BB,1,1);

    # compute separately LHC triplets impedance
    if E==4e12:
    	beta_filename_triplets=path_here+"LHC_elements/triplets_LHC_beta_length_B1_sq0p6m_3m_0p6m_3m.dat";
    elif E>=6.5e12:
    	beta_filename_triplets=path_here+"LHC_elements/triplets_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat";
    param_filename_triplets_BB=path_here+"LHC_elements/triplets_LHC_BB_param.dat";
    param_filename_triplets_holes=path_here+"LHC_elements/triplets_LHC_pumpingholes_param.dat";
    param_filename_triplets_RW=path_here+"LHC_elements/triplets_LHC_RW_param.dat";
    weld_filename=path_here+"LHC_elements/from_Carlo_weld_factor/presentscreen/weld_factor_current_arcBS_from_CZannini.dat";
    
    namesBB=read_ncol_file_identify_header(param_filename_triplets_BB,'name');
    namestaper=select_LHC_names(namesBB,pattern='taper');
    namesBPMs=select_LHC_names(namesBB,pattern='BPM');
    namesRW=read_ncol_file_identify_header(param_filename_triplets_RW,'name');
    namesBS=select_LHC_names(namesRW,pattern='BS');
    nameswarm=select_LHC_names(namesRW,pattern='Warm');
    
    print machine+" broad-band:";
    # Broad-Band contributions
    imp_mod_triplets_BB_taper,wake_mod_triplets_BB_taper=LHC_manyBB_resonator(avbetax,avbetay,param_filename_triplets_BB,
    	beta_filename_triplets,fcutoff=5e9,Q=1,beta=1,wake_calc=wake_calc,namesref=namestaper,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=[sigmaz,gamma,NbLHC,Mp]);
    imp_mod_triplets_BB_BPMs,wake_mod_triplets_BB_BPMs=LHC_manyBB_resonator(avbetax,avbetay,param_filename_triplets_BB,
    	beta_filename_triplets,fcutoff=5e9,Q=1,beta=1,wake_calc=wake_calc,namesref=namesBPMs,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=[sigmaz,gamma,NbLHC,Mp]);
    
    print machine+" holes:";
    # holes contributions
    imp_mod_triplets_holes=LHC_many_holes(avbetax,avbetay,param_filename_triplets_holes,
    	beta_filename_triplets,fcutoff=5e9,namesref=None,
	fpar=freq_param(ftypescan=0,nflog=100),power_loss_param=[sigmaz,gamma,NbLHC,Mp]);
    
    print machine+" resistive-wall:";
    # Resistive-Wall contributions
    # beam screen : with weld factor
    imp_mod_triplets_RW_BS,wake_mod_triplets_RW_BS=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_triplets_RW,
    	beta_filename_triplets,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=namesBS,
	power_loss_param=[sigmaz,NbLHC,Mp],freq_dep_factor_file=weld_filename,
	lxplusbatch=lxplusbatchImp,comment='_triplets_'+Estr,dire='triplets_'+Estr+'/');
    # warm pipe: no weld
    imp_mod_triplets_RW_warm,wake_mod_triplets_RW_warm=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_triplets_RW,
    	beta_filename_triplets,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=nameswarm,
	power_loss_param=[sigmaz,NbLHC,Mp],
	lxplusbatch=lxplusbatchImp,comment='_triplets_'+Estr,dire='triplets_'+Estr+'/');
	
    # total LHC triplets BB impedance
    imp_mod_triplets_BB=[];wake_mod_triplets_BB=[];
    add_impedance_wake(imp_mod_triplets_BB,imp_mod_triplets_BB_taper,1,1);
    add_impedance_wake(wake_mod_triplets_BB,wake_mod_triplets_BB_taper,1,1);
    add_impedance_wake(imp_mod_triplets_BB,imp_mod_triplets_BB_BPMs,1,1);
    add_impedance_wake(wake_mod_triplets_BB,wake_mod_triplets_BB_BPMs,1,1);
    
    # total LHC triplets RW impedance
    imp_mod_triplets_RW=[];wake_mod_triplets_RW=[];
    add_impedance_wake(imp_mod_triplets_RW,imp_mod_triplets_RW_BS,1,1);
    add_impedance_wake(wake_mod_triplets_RW,wake_mod_triplets_RW_BS,1,1);
    add_impedance_wake(imp_mod_triplets_RW,imp_mod_triplets_RW_warm,1,1);
    add_impedance_wake(wake_mod_triplets_RW,wake_mod_triplets_RW_warm,1,1);

    # total LHC triplets impedance
    imp_mod_triplets=[];wake_mod_triplets=[];
    add_impedance_wake(imp_mod_triplets,imp_mod_triplets_BB,1,1);
    add_impedance_wake(wake_mod_triplets,wake_mod_triplets_BB,1,1);
    add_impedance_wake(imp_mod_triplets,imp_mod_triplets_holes,1,1);
    add_impedance_wake(imp_mod_triplets,imp_mod_triplets_RW,1,1);
    add_impedance_wake(wake_mod_triplets,wake_mod_triplets_RW,1,1);
    
    if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):
	# dump into a file
	filemodel=open(root_result+'/impedances'+suffix+'.txt','w');
	pick.dump(imp_mod_triplets,filemodel);
	filemodel.close();

	# plot the various imp. contributions of triplets
	plot_compare_imp_model([imp_mod_triplets,imp_mod_triplets_RW,imp_mod_triplets_BB_taper,imp_mod_triplets_holes,imp_mod_triplets_BB_BPMs],
    	    [' total triplets',' resistive-wall',' tapers broad-band',' pumping holes',' BPMs broad-band'],
	    listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_contrib_'+machine+model+suffix,saveratio=root_result+'/plot_imp_ratio_contrib_'+machine+model+suffix);
    
    # now HL-LHC
    
    # some initializations
    tuneshiftQp=np.zeros((len(RWscan),len(BBscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);
    machine='HLLHC';print machine+": sigmaz=",sigmaz;
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    beta_filename_HLLHC_triplets=path_here+"LHC_elements/triplets_HLLHC_beta_length_B1_sq0p1m_10m_0p1m_3m.dat";

    for iRW,RW in enumerate(RWscan):

	param_filename_HLLHC_triplets_RW=path_here+"LHC_elements/triplets_HLLHC_RW_param"+RW+".dat";
    	weld_filename_HLLHC=path_here+"LHC_elements/from_Carlo_weld_factor/newscreen/weld_factor_HLLHC_pos3_2mm_from_CZannini.dat";

        print machine+" resistive-wall:",RW.replace('_',' ');
	# Resistive-Wall contributions
	if (RW.find('Cu')!=-1):
	    # apply weld factor
	    imp_mod_HLLHC_triplets_RW,wake_mod_HLLHC_triplets_RW=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_HLLHC_triplets_RW,
    		beta_filename_HLLHC_triplets,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,
		power_loss_param=[sigmaz,NbHLLHC,Mp],freq_dep_factor_file=weld_filename_HLLHC,
		lxplusbatch=lxplusbatchImp,comment='_'+machine+'_triplets_'+RW+Estr,dire='triplets_'+machine+'_'+Estr+'/');
	else:
	    # no weld factor since stainless steel
	    imp_mod_HLLHC_triplets_RW,wake_mod_HLLHC_triplets_RW=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_HLLHC_triplets_RW,
    		beta_filename_HLLHC_triplets,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,
		power_loss_param=[sigmaz,NbHLLHC,Mp],freq_dep_factor_file=None,
		lxplusbatch=lxplusbatchImp,comment='_'+machine+'_triplets_'+RW+Estr,dire='triplets_'+machine+'_'+Estr+'/');

	if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):

	    for iBB,BB in enumerate(BBscan):

    		param_filename_HLLHC_triplets_BB=path_here+"LHC_elements/triplets_HLLHC_BB_param"+BB+".dat";

		namesBB=read_ncol_file_identify_header(param_filename_HLLHC_triplets_BB,'name');
		namestaper=select_LHC_names(namesBB,pattern='taper');
		namesBPMs=select_LHC_names(namesBB,pattern='BPM');

        	print machine+" broad-band:",BB.replace('_',' ');
		# Broad-Band contributions
		imp_mod_HLLHC_triplets_BB_taper,wake_mod_HLLHC_triplets_BB_taper=LHC_manyBB_resonator(avbetax,avbetay,param_filename_HLLHC_triplets_BB,
    		    beta_filename_HLLHC_triplets,fcutoff=5e9,Q=1,beta=1,wake_calc=wake_calc,namesref=namestaper,
		    fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=[sigmaz,gamma,NbHLLHC,Mp]);
		imp_mod_HLLHC_triplets_BB_BPMs,wake_mod_HLLHC_triplets_BB_BPMs=LHC_manyBB_resonator(avbetax,avbetay,param_filename_HLLHC_triplets_BB,
    		    beta_filename_HLLHC_triplets,fcutoff=5e9,Q=1,beta=1,wake_calc=wake_calc,namesref=namesBPMs,
		    fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=[sigmaz,gamma,NbHLLHC,Mp]);

		# total LHC triplets BB impedance
		imp_mod_HLLHC_triplets_BB=[];wake_mod_HLLHC_triplets_BB=[];
		add_impedance_wake(imp_mod_HLLHC_triplets_BB,imp_mod_HLLHC_triplets_BB_taper,1,1);
		add_impedance_wake(wake_mod_HLLHC_triplets_BB,wake_mod_HLLHC_triplets_BB_taper,1,1);
		add_impedance_wake(imp_mod_HLLHC_triplets_BB,imp_mod_HLLHC_triplets_BB_BPMs,1,1);
		add_impedance_wake(wake_mod_HLLHC_triplets_BB,wake_mod_HLLHC_triplets_BB_BPMs,1,1);
    
		for ihole,hole in enumerate(holescan):

		    param_filename_HLLHC_triplets_holes=path_here+"LHC_elements/triplets_HLLHC_pumpingholes_param"+hole+".dat";

        	    print machine+" holes:",hole.replace('_',' ');
		    # holes contributions
		    imp_mod_HLLHC_triplets_holes=LHC_many_holes(avbetax,avbetay,param_filename_HLLHC_triplets_holes,
    			beta_filename_HLLHC_triplets,fcutoff=5e9,namesref=None,
			fpar=freq_param(ftypescan=0,nflog=100),power_loss_param=[sigmaz,gamma,NbHLLHC,Mp]);
    
		    # total LHC triplets impedance
		    imp_mod_HLLHC_triplets=[];wake_mod_HLLHC_triplets=[];
		    add_impedance_wake(imp_mod_HLLHC_triplets,imp_mod_HLLHC_triplets_BB,1,1);
		    add_impedance_wake(wake_mod_HLLHC_triplets,wake_mod_HLLHC_triplets_BB,1,1);
		    add_impedance_wake(imp_mod_HLLHC_triplets,imp_mod_HLLHC_triplets_holes,1,1);
		    add_impedance_wake(imp_mod_HLLHC_triplets,imp_mod_HLLHC_triplets_RW,1,1);
		    add_impedance_wake(wake_mod_HLLHC_triplets,wake_mod_HLLHC_triplets_RW,1,1);

		    # dump into a file
		    filemodel=open(root_result+'/impedances'+BB+hole+RW+'.txt','w');
		    pick.dump(imp_mod,filemodel);
		    filemodel.close();

		    # plot the various imp. contributions of HL-LHC triplets
		    plot_compare_imp_model([imp_mod_HLLHC_triplets,imp_mod_HLLHC_triplets_RW,imp_mod_HLLHC_triplets_BB_taper,imp_mod_HLLHC_triplets_holes,imp_mod_HLLHC_triplets_BB_BPMs],
    			[' total HL-LHC triplets',' resistive-wall',' tapers broad-band',' pumping holes',' BPMs broad-band'],
			listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_contrib_'+machine+BB+hole+RW+suffix,
			saveratio=root_result+'/plot_imp_ratio_contrib_'+machine+BB+hole+RW+suffix);

		    # plot LHC vs HL-LHC triplets impedance and ratio
		    plot_compare_imp_model([imp_mod_triplets,imp_mod_HLLHC_triplets],[' LHC triplets',' HL-LHC triplets'+BB.replace('_',' ')+hole.replace('_',' ')+RW.replace('_',' ')],
			listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_LHC_'+model+machine+BB+hole+RW+suffix,
			saveratio=root_result+'/plot_imp_ratio_LHC_'+model+machine+BB+hole+RW+suffix);

		    # plot LHC vs total impedance and ratio
		    plot_compare_imp_model([imp_mod,imp_mod_triplets],[' total'+model.replace('_',' '),' total LHC triplets'],
			listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_LHC'+model+suffix,
			saveratio=root_result+'/plot_imp_ratio_LHC'+model+suffix);

		    # plot HL-LHC vs total impedance and ratio
		    plot_compare_imp_model([imp_mod,imp_mod_HLLHC_triplets],[' total'+model.replace('_',' '),' total HL-LHC triplets'+BB.replace('_',' ')+hole.replace('_',' ')+RW.replace('_',' ')],
			listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_'+machine+BB+hole+RW+suffix,
			saveratio=root_result+'/plot_imp_ratio_'+machine+BB+hole+RW+suffix);


		    # DELPHI loops now
		    for iplane,plane in enumerate(['x','y']):
	        	# select Zxdip or Zydip
			for iw in imp_mod:
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

			    if flagDEL:
				tuneshiftQp[iRW,iBB,iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
					nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
					a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
					flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
					kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=5.e-3,
					lxplusbatch=lxplusbatchDEL,comment=machine+BB+RW+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
					dire=root_result+'/');



    # now the plots (outside loop on scenarios)
    if (flagDEL)and((lxplusbatchDEL.startswith('retrieve'))or(lxplusbatchDEL==None)):

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

			    for iBB,BB in enumerate(BBscan):

				for iRW,RW in enumerate(RWscan):

				    # output file name for data vs Qp
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+BB+RW+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[iBB,iRW,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Qpscan,sgn*ts,'DELPHI,'+BB.replace('_',' ')+RW.replace('_',' '),linetype[iBB]+col[iRW],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			# finish plot vs Qp
			for ir,r in enumerate(['real','imag']):
			    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
			    else: end_figure(figQp[ir],axQp[ir]);



    if not(flagsave): pylab.show();
