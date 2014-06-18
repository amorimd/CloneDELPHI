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
import time
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
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/LHC_coll_family_details';
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # coll. files definition
    scenario='_mm_kept';
    TCL6str='';TCL6gap=10.; # gap of TCL6 in sigmas
    #TCL6str='_TCL6_15sig';TCL6gap=15.; # gap of TCL6 in sigmas
    #TCL6str='_TCL6_20sig';TCL6gap=20.; # gap of TCL6 in sigmas
    param_filename_coll='../Coll_settings/collgaps_fromRoderik_modifNico_materialnames_feb2014'+scenario+TCL6str+'.dat';
    settings_filename_coll='../Coll_settings/collgaps_fromRoderik_modifNico_materialnames_feb2014'+scenario+TCL6str+'.dat';
    beta_filename_coll=param_filename_coll;
    squeeze='0p55m_10m_0p55m_10m';
    BPMflag=False;BPMstr='_noBPM';
    assym_fact_TCL6=1.;assym_fact_TCL6str=''; # symmetric TCL6 jaws
    #assym_fact_TCL6=round(40./TCL6gap*10)/10;assym_fact_TCL6str='_assymTCL6_40sig'; # other jaw of TCL6 at 40 sigmas
    
    
    # read coll. file
    namesref,material,thick,angle,length,halfgap,betax,betay=read_coll_files_several_mat(param_filename_coll,settings_filename_coll,beta_filename_coll,namesref=None);

    # create the coll. family scan
    # injection protection collimators
    inj_coll=select_LHC_names(namesref,pattern='TDI');
    inj_coll.extend(select_LHC_names(namesref,pattern='TCLI'));
    TCDQ=select_LHC_names(namesref,pattern='TCDQ');
    TCP=select_LHC_names(namesref,pattern='TCP');
    TCS=select_LHC_names(namesref,pattern='TCS');
    TCT=select_LHC_names(namesref,pattern='TCT');
    TCLA=select_LHC_names(namesref,pattern='TCLA');
    TCL=select_LHC_names(namesref,pattern='TCL.');
    TCL4=select_LHC_names(namesref,pattern='TCL.4');
    TCL6=select_LHC_names(namesref,pattern='TCL.6');
    names_scan=[inj_coll,TCDQ,TCP,TCS,TCT,TCLA,TCL,TCL4,TCL6];
    leg=['Total LHC imp. model','Inj. prot. coll.','TCDQ','TCP','TCS','TCT','TCLA','TCL','TCL4','TCL6'];
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    # compute total imp. model
    imp_mod_tot,wake_mod_tot=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
	    settings_filename_coll,TDIcoating='preLS1',dire="../LHC_elements/",
	    commentcoll=scenario,direcoll='Coll'+scenario+'/',
	    lxplusbatch=lxplusbatchImp,BPM=BPMflag,beam=beam,squeeze=squeeze,
	    wake_calc=wake_calc,assymetry_factor_TCL6=assym_fact_TCL6)

    imp_mod_list=[imp_mod_tot];
    wake_mod_list=[wake_mod_tot];
#    imp_mod_list=[];wake_mod_list=[];
    
    # loops to construct model for each collimator family
    for inames,names in enumerate(names_scan):
    #if True:
    #	names=TCL6;
	
	imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,
		avbetax,avbetay,param_filename_coll,settings_filename_coll,beta_filename_coll,
		wake_calc=wake_calc,ftypescan=0,nflog=100,zpar=z_param(),namesref=names,
		TDIcoating='preLS1',BPM=BPMflag,fcutoffBB=50e9,lxplusbatch=lxplusbatchImp,
		comment=scenario,dire='Coll'+scenario+'/',assymetry_factor_TCL6=assym_fact_TCL6);

	# add up
	imp_mod=[];wake_mod=[];
	add_impedance_wake(imp_mod,imp_mod_coll_RW,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll_RW,1,1);
	add_impedance_wake(imp_mod,imp_mod_coll_geom,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll_geom,1,1);
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot contributions of coll. families w.r.t the total impedance
	plot_compare_imp_model(imp_mod_list[:-2],leg[:-2],listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_LHC_v2_coll_family_details'+BPMstr+TCL6str+assym_fact_TCL6str,
	    saveratio=root_result+'/plot_imp_ratio_LHC_v2_coll_family_details'+BPMstr+TCL6str+assym_fact_TCL6str,
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3,plotpercent=True);
    
	# same with only TCL4 and TCL6
	imp_mod_list2=[imp_mod_list[0]];imp_mod_list2.extend(imp_mod_list[-2:]);
	leg2=[leg[0]];leg2.extend(leg[-2:]);
	plot_compare_imp_model(imp_mod_list2,leg2,listcomp=['Zlong','Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_LHC_v2_coll_TCL4_6_details'+BPMstr+TCL6str+assym_fact_TCL6str,
	    saveratio=root_result+'/plot_imp_ratio_LHC_v2_coll_TCL4_6_details'+BPMstr+TCL6str+assym_fact_TCL6str,
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3,plotpercent=True);
    
