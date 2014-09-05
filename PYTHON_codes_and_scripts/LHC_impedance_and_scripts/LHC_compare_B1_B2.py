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

    e,m0,c,E0=proton_param()
    E=4e12

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E)

    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_compare_B1_B2_BBcutoff50GHz'
    os.system("mkdir -p "+root_result);
    
    flagsave=1 # 1 to save figure instead of plotting on screen
    flagplot=True # to plot impedances

    wake_calc=False # True -> compute wake as well (otherwise only imp.)
       
    kmax=1 # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60 # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y'] # colors
    linetype=['-','--',':']

    # scan definition
    beamscan=np.array(['1','2'])
    strsubscan='_LHC_v2_B1_B2'
    margin_factor=1
    scenario='_2012_v2'
    legbeam=np.array(['LHC B1 2012','LHC B2 2012'])
    squeeze='0p6m_3m_0p6m_3m'
    
    param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'
    suffix_param='_2012.txt'
    
    settings_filename_coll=path_here+'Coll_settings/coll_settings_physics_fill_3265_B'
    suffix_settings='.txt'
    

    imp_mod_list=[] # complete list of impedance scenarios
    wake_mod_list=[]# complete list of wake scenarios
    
    for ibeam,beam in enumerate(beamscan):
    	
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E)
    	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	print "scenario: ",scenario,", beam: ",beam
	
	imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll+beam+suffix_param,
		settings_filename_coll+beam+suffix_settings,dire=path_here+"LHC_elements/",
		commentcoll=scenario,direcoll='Coll'+scenario+'/',
		lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeeze,
		wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave=scenario+'_B'+beam)
	
	imp_mod_list.append(imp_mod);
	wake_mod_list.append(wake_mod);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    #filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    #pick.dump(imp_mod,filemodel);
	    #filemodel.close();
	    
	    # write Ascii files with each component
	    write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
	    	dire=root_result+'/')
	    
	    if (wake_calc):
		# write Ascii files with each component
		write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
		    dire=root_result+'/')
	    
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# dip only
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)		
    

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	# plot and compare all scenarios with 2012 impedance
	maxratio_sb=plot_compare_imp_model(imp_mod_list,legbeam,listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
	    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
	    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
	    ylim=[1e5,1e10],ylim_ratio=[0,3],bounds=[40e6,2e9],legpos=3);
	
	# plot and compare all scenarios with 2012 wake
	if wake_calc:
	    maxratio_w=plot_compare_imp_model(wake_mod_list,legbeam,listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad'],
		saveimp=root_result+'/plot_wake_'+machine+"_scenarios"+strsubscan,
		saveratio=root_result+'/plot_wake_ratio_'+machine+"_scenarios"+strsubscan,
		xlim=[1e-5,1e6],ylim=[1e8,1e19],yliml=[1e6,1e19],bounds=[8e3,5e10],legpos=0);


    if not(flagsave): pylab.show();
