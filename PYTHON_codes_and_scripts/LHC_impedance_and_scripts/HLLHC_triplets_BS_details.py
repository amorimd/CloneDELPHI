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
from HLLHC_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/HLLHC_triplets_BS';
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=False;
    
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # coll. files definition
    param_filename_coll='../Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat';
    settings_filename_coll='../Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat';
    squeeze='0p15m_round';
    ftypescan=0;nflog=100;fcutoffBB=50e9;zpar=z_param();
    dire="../LHC_elements/"    
    
    leg=['Total','RW from triplet & D1 beam-screens','Tapers in triplets','Pumping holes in triplets & D1'];	
	
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    # compute total imp. model
    imp_mod_tot,wake_mod_tot=HLLHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
	    settings_filename_coll,dire="../LHC_elements/",
	    commentcoll='_HLLHC',direcoll='Coll_HLLHC_v2/',
	    lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeeze,
	    wake_calc=wake_calc,optionCrab='',margin_factor=1,
	    optionBBC=2,flagplot=flagplot,root_result=root_result,commentsave='_HLLHC')

    imp_mod_list=[imp_mod_tot];
    wake_mod_list=[wake_mod_tot];
    
    # beta functions for all the rest
    beta_filename_rest=dire+"HLLHC_beta_length_B"+str(beam)+"_sq"+squeeze+".dat"

    # compute model for the rest of the RW, with weld factor for beam screens
    if E>=1e12: Estr=float_to_str(E/1e12)+'TeV';
    else: Estr=float_to_str(E/1e9)+'GeV';
    param_filename_RW=dire+"HLLHC_RW_param_aC_Cu_ss.dat"
    weld_filename_triplets=dire+"from_Carlo_weld_factor/newscreen/weld_factor_HLLHC_pos2p8_2mm_from_CZannini_invert.dat"; # updated value
    
    namesRW=read_ncol_file_identify_header(param_filename_RW,'name');
    namesBS=select_LHC_names(namesRW,pattern='BS'); # all beam screens
    namesBS_triplets=select_LHC_names(namesBS,pattern='BS_1') # beam screens in triplets
    
    imp_mod_RW_BS_triplets,wake_mod_RW_BS_triplets=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesBS_triplets,
	power_loss_param=[7.5e-2,3.5e11,2*1404],freq_dep_factor_file=weld_filename_triplets,
	lxplusbatch=lxplusbatchImp,comment='_'+Estr,dire='BS_HLLHC_v2_'+Estr+'/');

    imp_mod_list.append(imp_mod_RW_BS_triplets);
    wake_mod_list.append(wake_mod_RW_BS_triplets);
    
    # test with more Cu
    param_filename_RW2=dire+"HLLHC_RW_param_aC_Cu_ss_moreCu.dat"
    imp_mod_RW_BS_triplets2,wake_mod_RW_BS_triplets2=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW2,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesBS_triplets,
	power_loss_param=[7.5e-2,3.5e11,2*1404],freq_dep_factor_file=weld_filename_triplets,
	lxplusbatch=lxplusbatchImp,comment='_'+Estr+'_100mum_Cu',dire='BS_HLLHC_v2_'+Estr+'/');


    # individual broad-band contributions
    param_filename_BB=dire+'triplets_HLLHC_BB_param.dat';
    namesBB=read_ncol_file_identify_header(param_filename_BB,'name');
    namestaper=select_LHC_names(namesBB,pattern='taper');
    
    imp_mod_triplets_BB_taper,wake_mod_triplets_BB_taper=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB,
    	beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=namestaper,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
        
    imp_mod_list.append(imp_mod_triplets_BB_taper);
    wake_mod_list.append(wake_mod_triplets_BB_taper);
	
    
    # holes contributions
    param_filename_holes=dire+'HLLHC_pumpingholes_param.dat';
    names_holes=read_ncol_file_identify_header(param_filename_holes,'[nM]ame');
    names_holes_triplets=select_LHC_names(names_holes,pattern='BS_1');
    names_holes_rest=invert_selection(names_holes,names_holes_triplets);
    
    imp_mod_holes_triplets,wake_mod_holes_triplets=LHC_many_holes(avbetax,avbetay,param_filename_holes,
    	beta_filename_rest,fcutoff=fcutoffBB,namesref=names_holes_triplets,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog));

    imp_mod_list.append(imp_mod_holes_triplets);
    wake_mod_list.append(wake_mod_holes_triplets);
	
    param_filename_holes2=dire+'HLLHC_pumpingholes_param_feb2014_RKersevan.dat';

    imp_mod_holes_triplets2,wake_mod_holes_triplets2=LHC_many_holes(avbetax,avbetay,param_filename_holes2,
    	beta_filename_rest,fcutoff=fcutoffBB,namesref=names_holes_triplets,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog));
	

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot contributions of coll. families w.r.t the total impedance
	maxratio=plot_compare_imp_model(imp_mod_list,leg,listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_HLLHC_v2_triplets_details',
	    saveratio=root_result+'/plot_imp_ratio_HLLHC_v2_triplets_details',
	    xlim=[1e3,5e9],ylim=[1e1,1e10],bounds=[8e3,2e9],legpos=1,plotpercent=True,legpercentpos=(1,0.9),maxpercent=5);
    
	print "all triplets stuff: ratios ",maxratio;
	
	# same with only holes (pessimistic case from R. Kersevan)
	imp_mod_list2=[imp_mod_list[0],imp_mod_holes_triplets2];
	leg2=[leg[0],leg[-1]+", pessimistic case"];
	maxratio_holes=plot_compare_imp_model(imp_mod_list2,leg2,listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_HLLHC_v2_triplets_holes_pessimistic_details',
	    saveratio=root_result+'/plot_imp_ratio_HLLHC_v2_triplets_holes_pessimistic_details',
	    xlim=[1e3,5e9],ylim=[1e1,1e10],bounds=[8e3,2e9],legpos=1,plotpercent=True,legpercentpos=(1,0.9),maxpercent=1);
    
	print "holes triplets, pessimistic: ratio ",maxratio_holes;
	
	# same with only RW (more copper)
	imp_mod_list3=[imp_mod_list[0],imp_mod_RW_BS_triplets2];
	leg3=[leg[0],leg[1]+", 100 $ \mu $m Cu"];
	maxratio_100mumCu=plot_compare_imp_model(imp_mod_list3,leg3,listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_HLLHC_v2_triplets_RW_100mumCu_details',
	    saveratio=root_result+'/plot_imp_ratio_HLLHC_v2_triplets_RW_100mumCu_details',
	    xlim=[1e3,5e9],ylim=[1e1,1e10],bounds=[8e3,2e9],legpos=1,plotpercent=True,legpercentpos=(1,0.9),maxpercent=5);
    
	print "RW, with 100mum Cu ",maxratio_100mumCu;
	
