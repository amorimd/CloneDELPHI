#!/usr/bin/python

import sys
from string import *
import numpy as np
import pickle as pick
import pylab,re
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from string_lib import *
from Impedance import *
from LHC_imp import *
from LHC_coll_imp import *


def HLLHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll=None,TDIcoating='preLS1',dire="../LHC_elements/",commentcoll='',direcoll='Coll_HLLHC_v2/',lxplusbatch=None,
	BPM=True,beam='1',squeeze='0p15m_round',wake_calc=False,ftypescan=0,nflog=100,zpar=z_param(),
	optionCrab=None,optionBBC=0,margin_factor=1.,optionMo_TCS37=None,optionTCu_triplet='',
	fcutoffBB=50e9,flagplot=False,root_result='../DELPHI_results/LHC',commentsave=''):

    ''' Total HL-LHC impedance model as of Nov. 2013 (model v1 was actually
    the LHC model v1 with HL-LHC collimators and beta functions)
    - E: particle energy in eV
    - avbetax, avbetay: average beta function from smooth approximation (R/Q) in m
    - param_filename_coll: file with collimator parameters (including materials, 
      angle, length, beta functions)
    - settings_filename_coll: collimator settings file
    - beta_filename_coll: collimator beta functions file (if present - otherwise take param_filename_coll)
    - TDIcoating: kind of coatgin for first block of TDI: can be 'preLS1'
      (5mum Ti) or 'postLS1' (1mum NEG + 2mum CU + 0.3mum NEG + 5mum Ti) or 
      directly a list of layer objects
    - dire: directory where to find all files with parameters for the rest of the machine
    - commentcoll: comment for the collimators IW2D computation
    - direcoll: subdirectory of ImpedanceWake2D where to put collimators impedance
    - lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus
    		   if 'retrieve' -> retrieve outputs
    - BPM: coll. geometry contains a BPM cavity if BPM is True, otherwise old LHC coll. geometry
    - beam: beam number ('1' or '2')
    - squeeze: suffix of filename with beta functions, for the rest of the 
    machine (everything except collimators). The first number gives the squeeze in IP1, 
    from which the broad-band impedance (from design report) is evaluated (squeezed if <2m,
    otherwise injection BB model is taken)
    - wake_calc: True to compute wake function as well
    - ftypescan, nflog and zpar: parameters for frequency and distance scans
    - optionCrab: None to put no crab cavities, or:
    	'' to put crab cavities as 2 broad-band models,
	'_BNL' to put crab cavities as lsit of HOMS, for BNL kind of cavities,
	'_SLAC' to put crab cavities as lsit of HOMS, for SLAC kind of cavities,
    - optionBBC: 0 -> no beam-beam wire compensator,
    		 1 -> beam-beam wire compensator as stand-alone wire (stripline BPM model),
		 2 -> beam-beam wire compensator embedded in tungsten collimator,
    - optionMo_TCS37: if not None, indicate Mo coating thickness on TCS collimators
    in IR3 and IR7,
    - option_TCu_triplet: if '', copper temperature in beam screen triplets is 20K, otherwise
    this temperature should be indicated (e.g. '50K' -> only case implemented now),
    - margin_factor: additional factor to take into account all the unknowns. 
    	-> We multiply the final model by this value,
    - fcutoffBB: cutoff frequency for broad-band models,
    - flagplot: if True, plot impedances and percent of each part of the model,
    - root_result: used only with flagplot: directory where to put impedance plots,
    - commentsave: used only with flagplot: additional comment for filename of impedance plots.
    '''

    imp_mod=[];wake_mod=[];
    
    # compute model for collimators
    if beta_filename_coll==None: beta_filename_coll=param_filename_coll;
    
    if optionMo_TCS37==None:
	imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,zpar=zpar,namesref=None,
	    TDIcoating=TDIcoating,BPM=BPM,fcutoffBB=fcutoffBB,lxplusbatch=lxplusbatch,comment=commentcoll,dire=direcoll);
    else:
	# case with some metallic Mo coating on the TCS in IR3 and IR7
	# select coll. names
	namestot=read_ncol_file_identify_header(param_filename_coll,'name');
	namesref=select_LHC_coll_IR(namestot,pattern='TCS',IRlist=[3,7]);print namesref
	# anti-selection
	namesanti=invert_selection(namestot,namesref)
	# model for Mo-coated coll.
	imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,zpar=zpar,namesref=namesref,coatingmat='Mo',coatingthickness=optionMo_TCS37,
	    TDIcoating=TDIcoating,BPM=BPM,fcutoffBB=fcutoffBB,lxplusbatch=lxplusbatch,comment=commentcoll+'_Mo'+float_to_str(optionMo_TCS37*1e6)+'mum',dire=direcoll);
	# model for the other coll.
	imp_mod_coll_RW2,wake_mod_coll_RW2,imp_mod_coll_geom2,wake_mod_coll_geom2=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,zpar=zpar,namesref=namesanti,
	    TDIcoating=TDIcoating,BPM=BPM,fcutoffBB=fcutoffBB,lxplusbatch=lxplusbatch,comment=commentcoll,dire=direcoll);	
	add_impedance_wake(imp_mod_coll_RW,imp_mod_coll_RW2,1,1);
	add_impedance_wake(wake_mod_coll_RW,wake_mod_coll_RW2,1,1);
	add_impedance_wake(imp_mod_coll_geom,imp_mod_coll_geom2,1,1);
	add_impedance_wake(wake_mod_coll_geom,wake_mod_coll_geom2,1,1);
    	
    
    # beta functions for all the rest
    beta_filename_rest=dire+"HLLHC_beta_length_B"+str(beam)+"_sq"+squeeze+".dat"
    
    isq=squeeze.find('m');
    if (float(squeeze[:isq].replace('p','.'))<0.5): Bfield_triplet=11.;Bfieldstr='_B11T'; # use B=11T in triplets when squeezed, for magneto-resistance effect (from E. Todesco)
    else: Bfield_triplet=None;Bfieldstr=''; # B field computed from energy with LHC bending radius (not really correct, but small effect anyway...)

    # compute model for the rest of the RW, with weld factor for beam screens
    if E>=1e12: Estr=float_to_str(E/1e12)+'TeV';
    else: Estr=float_to_str(E/1e9)+'GeV';
    param_filename_RW=dire+"HLLHC_RW_param_aC_Cu"+optionTCu_triplet+"_ss.dat"
    #weld_filename=dire+"from_Carlo_weld_factor/newscreen/weld_factor_HLLHC_pos3_2mm_from_CZannini_invert_xy_WARNING.dat";
    weld_filename_triplets=dire+"from_Carlo_weld_factor/newscreen/weld_factor_HLLHC_pos2p8_2mm_from_CZannini_invert.dat"; # updated value
    weld_filename_rest=dire+"from_Carlo_weld_factor/presentscreen/weld_factor_current_arcBS_from_CZannini.dat";
    
    namesRW=read_ncol_file_identify_header(param_filename_RW,'name');
    namesBS=select_LHC_names(namesRW,pattern='BS'); # all beam screens
    namesBS_triplets=select_LHC_names(namesBS,pattern='BS_1') # beam screens in triplets
    namesBS_rest=invert_selection(namesBS,namesBS_triplets) # other beam screens (those as in the LHC)
    namesrest=invert_selection(namesRW,namesBS) # the rest
    #print namesRW,namesBS,namesBS_triplets,namesBS_rest,namesrest
    
    imp_mod_RW_BS_triplets,wake_mod_RW_BS_triplets=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,Bfield=Bfield_triplet,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesBS_triplets,
	freq_dep_factor_file=weld_filename_triplets,
	lxplusbatch=lxplusbatch,comment=optionTCu_triplet+'_'+Estr+Bfieldstr,dire='BS_HLLHC_v2_'+Estr+'/');

    imp_mod_RW_BS_rest,wake_mod_RW_BS_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesBS_rest,
	freq_dep_factor_file=weld_filename_rest,
	lxplusbatch=lxplusbatch,comment='_'+Estr,dire='BS_HLLHC_v2_'+Estr+'/');

    # total RW beam screens impedance
    imp_mod_RW_BS=[];wake_mod_RW_BS=[];
    add_impedance_wake(imp_mod_RW_BS,imp_mod_RW_BS_triplets,1,1);
    add_impedance_wake(wake_mod_RW_BS,wake_mod_RW_BS_triplets,1,1);
    add_impedance_wake(imp_mod_RW_BS,imp_mod_RW_BS_rest,1,1);
    add_impedance_wake(wake_mod_RW_BS,wake_mod_RW_BS_rest,1,1);
    
    imp_mod_RW_rest,wake_mod_RW_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesrest,
	lxplusbatch=lxplusbatch,comment='_'+Estr,dire='Rest_HLLHC_v2_'+Estr+'/');
    
    # individual broad-band contributions
    param_filename_BB=dire+'triplets_HLLHC_BB_param_50GHz.dat';
    namesBB=read_ncol_file_identify_header(param_filename_BB,'name');
    namestaper=select_LHC_names(namesBB,pattern='taper');
    namesBPMs=select_LHC_names(namesBB,pattern='BPM');
    
    imp_mod_triplets_BB_taper,wake_mod_triplets_BB_taper=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB,
    	beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=namestaper,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    
    imp_mod_triplets_BB_BPMs,wake_mod_triplets_BB_BPMs=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB,
    	beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=namesBPMs,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    
    # HOMs contributions
    # WARNING: check RF, and check convention for shumt impedance (omega/c...)
    param_filename_HOMs=dire+'HLLHC_HOMs_param.dat';
    imp_mod_HOMs,wake_mod_HOMs=LHC_many_resonators(avbetax,avbetay,param_filename_HOMs,
    	beta_filename_rest,beta=1,wake_calc=wake_calc,namesref=None,
	fpar=freq_param(ftypescan=1,fmin=3e8,fmax=3e9,fsamplin=2e4,fadded=np.concatenate((10**np.arange(2,7,0.2),np.arange(1e7,5e10+1e7,1e7),[1e12]))),
	zpar=zpar);
    
    # holes contributions
    param_filename_holes=dire+'HLLHC_pumpingholes_param.dat';
    names_holes=read_ncol_file_identify_header(param_filename_holes,'[nM]ame');
    names_holes_triplets=select_LHC_names(names_holes,pattern='BS_1');
    names_holes_rest=invert_selection(names_holes,names_holes_triplets);
    
    imp_mod_holes_triplets,wake_mod_holes_triplets=LHC_many_holes(avbetax,avbetay,param_filename_holes,
    	beta_filename_rest,fcutoff=fcutoffBB,namesref=names_holes_triplets,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog));
    imp_mod_holes_rest,wake_mod_holes_rest=LHC_many_holes(avbetax,avbetay,param_filename_holes,
    	beta_filename_rest,fcutoff=fcutoffBB,namesref=names_holes_rest,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog));
    
    # compute broad-band model from design (the contributions not already taken into account elsewhere)
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband_less(wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar,fcutoffBB=fcutoffBB);

    # add up
    imp_mod=[];wake_mod=[];
    add_impedance_wake(imp_mod,imp_mod_coll_RW,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll_RW,1,1);

    add_impedance_wake(imp_mod,imp_mod_coll_geom,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll_geom,1,1);

    add_impedance_wake(imp_mod,imp_mod_RW_BS,1,1);
    add_impedance_wake(wake_mod,wake_mod_RW_BS,1,1);

    add_impedance_wake(imp_mod,imp_mod_RW_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_RW_rest,1,1);

    add_impedance_wake(imp_mod,imp_mod_triplets_BB_taper,1,1);
    add_impedance_wake(wake_mod,wake_mod_triplets_BB_taper,1,1);

    add_impedance_wake(imp_mod,imp_mod_triplets_BB_BPMs,1,1);
    add_impedance_wake(wake_mod,wake_mod_triplets_BB_BPMs,1,1);

    add_impedance_wake(imp_mod,imp_mod_HOMs,1,1);
    add_impedance_wake(wake_mod,wake_mod_HOMs,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod,wake_mod_BB,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_holes_triplets,1,1);
    add_impedance_wake(wake_mod,wake_mod_holes_triplets,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_holes_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_holes_rest,1,1);
    
    # additional optional stuff (Crab cavities, beam-beam wire compensator)
    # default value
    noption=0;
    imp_mod_BB_Crab=[];wake_mod_BB_Crab=[];
    imp_mod_BBC=[];wake_mod_BBC=[];
    
    if (optionCrab!=None):
	noption += 1;
        # OLD VERSION: crab cavities as one broad-band model
    	#param_filename_BB_Crab=dire+'HLLHC_BB_param.dat';
	#imp_mod_BB_Crab,wake_mod_BB_Crab=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB_Crab,
    	#    beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=['Crab'],
	#    fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
	
	# crab cavities as list of HOMs for specific kind of cavities (for now, '', '_SLAC' or '_BNL')
	param_filename_HOMs_Crab=dire+'HLLHC_crab_HOMs_param'+optionCrab+'.dat';
	imp_mod_BB_Crab,wake_mod_BB_Crab=LHC_many_resonators(avbetax,avbetay,param_filename_HOMs_Crab,
    	    beta_filename_rest,beta=1,wake_calc=wake_calc,namesref=None,
	    fpar=freq_param(ftypescan=1,fmin=3e8,fmax=3e9,fsamplin=2e4,fadded=np.concatenate((10**np.arange(2,7,0.2),np.arange(1e7,5e10+1e7,1e7),[1e12]))),
	    zpar=zpar);

	add_impedance_wake(imp_mod,imp_mod_BB_Crab,1,1);
	add_impedance_wake(wake_mod,wake_mod_BB_Crab,1,1);
    
	    
    if optionBBC>=1:
        # beam-beam wire compensators
	noption += 1;
	
	if optionBBC==1:
	    # Here modelled as stripline BPMs, assuming
	    # it is at 9.5 sigma with beta=4.34819e+03, 7 TeV, 2mm.mrad emittances,
	    # and wire radius=1mm)
	    param_filename_BBC=dire+'HLLHC_stripline_param.dat';
    	    imp_mod_BBC,wake_mod_BBC=LHC_many_striplineBPMs(avbetax,avbetay,param_filename_BBC,
		    beta_filename_rest,beta=1,wake_calc=wake_calc,namesref=['BBC'],
		    fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);

	elif optionBBC==2:
	    # here wire considered as embedded in a tungsten collimator at 9.5 sigmas,
	    # eps=2 mm.mrad, 7 TeV, taking minimum possible beta functions to compute half-gap)
	    # (see above)
	    param_filename_BBC=dire+'HLLHC_TCT_BBC_param.dat';
	    
	    imp_mod_BBC_RW,wake_mod_BBC_RW,imp_mod_BBC_geom,wake_mod_BBC_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,param_filename_BBC,param_filename_BBC,
		beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,zpar=zpar,namesref=None,
		BPM=False,fcutoffBB=fcutoffBB,lxplusbatch=lxplusbatch,comment=commentcoll,dire=direcoll);

	    add_impedance_wake(imp_mod_BBC,imp_mod_BBC_RW,1,1);
	    add_impedance_wake(wake_mod_BBC,wake_mod_BBC_RW,1,1);    
	    add_impedance_wake(imp_mod_BBC,imp_mod_BBC_geom,1,1);
	    add_impedance_wake(wake_mod_BBC,wake_mod_BBC_geom,1,1);    

	
	add_impedance_wake(imp_mod,imp_mod_BBC,1,1);
	add_impedance_wake(wake_mod,wake_mod_BBC,1,1);    

    # margin factor
    multiply_impedance_wake(imp_mod,margin_factor);
    multiply_impedance_wake(wake_mod,margin_factor);

    
    if flagplot:
    	imp_mod_list=[imp_mod,imp_mod_coll_RW,imp_mod_coll_geom,imp_mod_RW_BS,
		imp_mod_RW_rest,imp_mod_triplets_BB_taper,
		imp_mod_triplets_BB_BPMs,imp_mod_HOMs,imp_mod_holes_triplets,
		imp_mod_holes_rest,imp_mod_BB,imp_mod_BB_Crab,imp_mod_BBC];
	
	for iw in imp_mod_list: multiply_impedance_wake(iw,margin_factor);

	leg=['Total','RW from coll','Geom. from coll','RW from beam-screen',
		'RW from warm pipe','Tapers in triplets','BPMs in triplets',
		'RF, ATLAS, CMS, ALICE & LHCb','Pumping holes (triplets)',
		'Pumping holes (rest)','Other broad-band contributions',
		'Crab cavities','Wire compensators'];	
	
	plot_compare_imp_model(imp_mod_list[:len(imp_mod_list)-(2-noption)],
	    leg[:len(imp_mod_list)-(2-noption)],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad'],
	    saveimp=root_result+'/plot_imp_HLLHC_v2_'+commentsave+'details',
	    saveratio=root_result+'/plot_imp_ratio_HLLHC_v2_'+commentsave+'details',
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3,plotpercent=True,legpercentpos=(0.8,0.75));

	# plot of the total
	plot_compare_imp_model([imp_mod],[''],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad','Zxcst','Zycst'],
	    saveimp=root_result+'/plot_imp_HLLHC_v2_'+commentsave,
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3);

	if wake_calc:

    	    wake_mod_list=[wake_mod,wake_mod_coll_RW,wake_mod_coll_geom,wake_mod_RW_BS,
		    wake_mod_RW_rest,wake_mod_triplets_BB_taper,
		    wake_mod_triplets_BB_BPMs,wake_mod_HOMs,wake_mod_holes_triplets,
		    wake_mod_holes_rest,wake_mod_BB,wake_mod_BB_Crab,wake_mod_BBC];
	
	    for iw in wake_mod_list: multiply_impedance_wake(iw,margin_factor);

	    plot_compare_imp_model(wake_mod_list[:len(imp_mod_list)-(2-noption)],
	        leg[:len(imp_mod_list)-(2-noption)],listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad'],
		saveimp=root_result+'/plot_wake_HLLHC_v2_'+commentsave+'details',
		saveratio=root_result+'/plot_wake_ratio_HLLHC_v2_'+commentsave+'details',
		xlim=[1e-1,1e6],ylim=[1e12,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3,
		plotpercent=True,legpercentpos=(0.8,0.8));
	    
	    # plot of the total
	    plot_compare_imp_model([wake_mod],[''],listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad','Wxcst','Wycst'],
		saveimp=root_result+'/plot_wake_HLLHC_v2_'+commentsave,
		xlim=[1e-2,1e6],ylim=[1e10,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3);
	

    return imp_mod,wake_mod;
    
