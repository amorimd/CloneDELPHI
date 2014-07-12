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
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_optim_IR3';
    os.system("mkdir -p "+root_result);
    
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
    scenario='_mm_kept';
    param_filename_coll=path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames'+scenario+'.dat';
    settings_filename_coll_prefix=path_here+'Coll_settings/../Coll_settings/collgaps_settings_sigma'+scenario+'_offsetIR3_';
    beta_filename_coll=param_filename_coll;
    squeeze='0p6m_3m_0p6m_3m';
    fminrefine=1.e11;fmaxrefine=5.e12;nrefine=5000;ftypescan=2;nflog=100;
    BPMflag=True

    # sigma offset scan
    sigmascan=['0sig','-4sig','-3sig','-2sig','-1sig','1sig','2sig','3sig','4sig'];

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
    
    imp_mod_list=[]; # complete list of impedance scenarios
    
    for isig,sig in enumerate(sigmascan):
    
	settings_filename_coll=settings_filename_coll_prefix+sig+'.dat';

	# compute total imp. model
	imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
		settings_filename_coll,TDIcoating='preLS1',dire=path_here+"LHC_elements/",
		commentcoll=scenario,direcoll='Coll'+scenario+'/',
		lxplusbatch=lxplusbatchImp,BPM=BPMflag,beam=beam,squeeze=squeeze,
		wake_calc=wake_calc)
	
	imp_mod_list.append(imp_mod);
    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios
	maxratio=plot_compare_imp_model(imp_mod_list,sigmascan,listcomp=['Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_IR3_sig_scenarios",
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_IR3_sig_scenarios",xlim=[1e3,1e12],legpos=[0.8,0.08]);


	# check details of coll. impedance

	# compute impedance of all collimators
	settings_filename_coll=settings_filename_coll_prefix+'0sig.dat';

	imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,
    	    param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=None,
	    TDIcoating='preLS1',BPM=BPMflag,lxplusbatch=lxplusbatchImp,comment=scenario,dire='Coll'+scenario+'/');

	imp_mod_coll=[];
	add_impedance_wake(imp_mod_coll,imp_mod_coll_RW,1,1);
	add_impedance_wake(imp_mod_coll,imp_mod_coll_geom,1,1);

	# imp. of coll. in IR3
	namesref,material,angle,length,halfgap,betax,betay=read_coll_files(param_filename_coll,settings_filename_coll,beta_filename_coll,namesref=None);
	namesIR3=select_LHC_coll_IR(namesref,IRlist=[3]);print "namesIR3=",namesIR3,len(namesIR3)

	imp_mod_coll_IR3_RW,wake_mod_coll_IR3_RW,imp_mod_coll_IR3_geom,wake_mod_coll_IR3_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,
    	    param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesIR3,
	    TDIcoating='preLS1',BPM=BPMflag,lxplusbatch=lxplusbatchImp,comment=scenario,dire='Coll'+scenario+'/');

	imp_mod_coll_IR3=[];
	add_impedance_wake(imp_mod_coll_IR3,imp_mod_coll_IR3_RW,1,1);
	add_impedance_wake(imp_mod_coll_IR3,imp_mod_coll_IR3_geom,1,1);

	# other collimators
	namesanti=invert_selection(namesref,namesIR3);print "namesanti=",namesanti,len(namesanti)
	print "namesref=",namesref,len(namesref)

	imp_mod_coll_anti_RW,wake_mod_coll_anti_RW,imp_mod_coll_anti_geom,wake_mod_coll_anti_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,
    	    param_filename_coll,settings_filename_coll,
	    beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesanti,
	    TDIcoating='preLS1',BPM=BPMflag,lxplusbatch=lxplusbatchImp,comment=scenario,dire='Coll'+scenario+'/');

	imp_mod_coll_anti=[];
	add_impedance_wake(imp_mod_coll_anti,imp_mod_coll_anti_RW,1,1);
	add_impedance_wake(imp_mod_coll_anti,imp_mod_coll_anti_geom,1,1);

	imp_mod_list_0sig=[imp_mod_coll,imp_mod_coll_IR3,imp_mod_coll_anti]


	# plot and compare
	maxratio=plot_compare_imp_model(imp_mod_list_0sig,['Total coll.','IR3 only','Other coll.'],listcomp=['Zxdip','Zydip'],
	    saveimp=root_result+'/plot_imp_'+machine+"_coll_details",
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_coll_details",plotpercent=True);


    if not(flagsave): pylab.show();
