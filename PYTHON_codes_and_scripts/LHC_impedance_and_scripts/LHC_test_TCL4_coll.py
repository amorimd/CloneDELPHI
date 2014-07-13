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
import numpy as np
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param
from LHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    E=6.5e12;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    flagsave=1;
    
    # scan definition
    scenarioscan=['','_onlyTCL4'];
    model='_relaxed';
    root_result=path_here+'../../../DELPHI_results/'+machine+'/TCL4';
    os.system("mkdir -p "+root_result);
    suffix='_only_coll' # suffix for output files 
 
    # compute non changing part of the impedance
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
    if len(sys.argv)>1: lxplusbatch=str(sys.argv[1]);
    else: lxplusbatch=None;
    print lxplusbatch
        
    # rest (i.e. not coll.) of the machine wall impedance
    param_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_param.dat"
    #beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p6m_3m_0p6m_3m.dat"
    beta_filename_rest=path_here+"LHC_elements/beam_screens_warm_pipe_LHC_beta_length_B1_sq0p55m_10m_0p55m_10m.dat"
    
    imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
    	wake_calc=False,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,
	comment='_'+Estr,dire='Rest_'+Estr+'/');
    
    # broad-band model
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=False,
    	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());
    

    param_filename_coll=path_here+"Coll_settings/collgaps_fromRoderik_modifNico_materialnames"+model+".dat";
    beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;

    # compute model for all collimators
    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,
	comment=model,dire='Coll'+model+'/');

    # compute model for TCL4 alone
    imp_mod_TCL4,wake_mod_TCL4=LHC_singlecoll_iw_model('TCL4','CU',5e-3,0,gamma,1,
    	wake_calc=False,fpar=freq_param(),zpar=z_param(),lxplusbatch=lxplusbatch,comment='',dire='Coll_TCL4/')

    # beta function weighting (560m in both planes: pessimistic)
    imp_mod_TCL4weighted=[];wake_mod_TCL4weighted=[];
    add_impedance_wake(imp_mod_TCL4weighted,imp_mod_TCL4,560/avbetax,560/avbetay);
    add_impedance_wake(wake_mod_TCL4weighted,wake_mod_TCL4,560/avbetax,560/avbetay);

    # add up
    imp_mod=[];wake_mod=[];
    add_impedance_wake(imp_mod,imp_mod_coll,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll,1,1);
    add_impedance_wake(imp_mod,imp_mod_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_rest,1,1);
    add_impedance_wake(imp_mod,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod,wake_mod_BB,1,1);

    if (lxplusbatch==None)or(lxplusbatch.startswith('retrieve')):

	# dump into a file
	filemodel=open('impedances'+model+suffix+'.txt','w');
	pick.dump(imp_mod,filemodel);
	pick.dump(imp_mod_TCL4weighted,filemodel);
	pick.dump(imp_mod_coll,filemodel);
	filemodel.close();

	#compare_vs_zbase()

	plot_compare_imp_model([imp_mod_coll,imp_mod_TCL4],
    	    [' total collimators',r' one TCL4 with $ \beta=560 $ m'],
	    listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_contrib_'+machine+suffix+'_TCL4',
	    saveratio=root_result+'/plot_imp_ratio_contrib_'+machine+suffix+'_TCL4');

	plot_compare_imp_model([imp_mod,imp_mod_TCL4],
    	    [' total impedance',r' one TCL4 with $ \beta=560 $ m'],
	    listcomp=['Zlong','Zxdip','Zydip'],saveimp=root_result+'/plot_imp_contrib_'+machine+'_TCL4',
	    saveratio=root_result+'/plot_imp_ratio_contrib_'+machine+'_TCL4');

    if not(flagsave): pylab.show();
