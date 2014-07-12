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

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    flagsave=1;
    
    # scan definition
    scenarioscan=['','_onlyTCRYO'];
    model='_nominal';
    root_result=path_here+'../../../DELPHI_results/'+machine+'/TCRYO';
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
    
    #imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,beta_filename_rest,
    #	wake_calc=False,ftypescan=0,nflog=100,namesref=None,lxplusbatch=lxplusbatchImp,comment='_'+Estr);
    
    # broad-band model
    #imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=True,wake_calc=False,
    #	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param());
    

    figratio,axratio=init_figure(); # figure for impedance ratio
    fig,ax=init_figure();
   
    for iscenario,scenario in enumerate(scenarioscan):
    
    	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    	param_filename_coll=path_here+"Coll_settings/collgaps_fromRoderik_modifNico_materialnames"+model+".dat";
        beta_filename_coll=param_filename_coll;settings_filename_coll=param_filename_coll;
	
	# select coll. names
	namestot=read_ncol_file_identify_header(param_filename_coll,'name');
	namesref=select_LHC_coll_IR(namestot,pattern='TCRYO',IRlist=[7]);

	# compute model for collimators
	if iscenario==0:
	    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
		beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=namestot,lxplusbatch=lxplusbatch,
		comment=model,dire='Coll'+model+'/');
	else:
	    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
		beta_filename_coll,wake_calc=wake_calc,ftypescan=0,nflog=100,namesref=namesref,
		lxplusbatch=lxplusbatch,comment=model,dire='Coll'+model+'/');	

	# add up
	imp_mod=[];wake_mod=[];
	add_impedance_wake(imp_mod,imp_mod_coll,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll,1,1);
	#add_impedance_wake(imp_mod,imp_mod_rest,1,1);
	#add_impedance_wake(wake_mod,wake_mod_rest,1,1);
	#add_impedance_wake(imp_mod,imp_mod_BB,1,1);
	#add_impedance_wake(wake_mod,wake_mod_BB,1,1);
	
	if (iscenario==0): imp_mod0=deepcopy(imp_mod);wake_mod0=deepcopy(wake_mod);
    
	if (lxplusbatch.startswith('retrieve'))or(lxplusbatch==None):

	    # dump into a file
	    filemodel=open('impedances'+model+scenario+'.txt','w');
	    pick.dump(imp_mod,filemodel);
	    filemodel.close();

	    #compare_vs_zbase()

	    # first plot all components
	    pat=['-','--'];units=["","/m","/m$^2$","/m$^3$","/m$^4$"];
	    listcomp=['Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad'];
	    # corresponding a,b,c,d and planes in imp_mod
	    lista=[1,0,0,0,0,0];listb=[0,1,0,0,1,0];listc=[0,0,1,0,0,0];listd=[0,0,0,1,0,1];
	    listplane=['x','y','x','y','x','x'];

            #fig,ax=init_figure();
	    for icomp,comp in enumerate(listcomp[:2]):
		# find corresponding term in imp_mod
		for kiw,iw in enumerate(imp_mod):
		    if test_impedance_wake_comp(iw,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw_comp=kiw;
		for ir,r in enumerate(['Re','Im']):
        	    plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,ir]),r+'('+comp+'), coll. imp. post LS1,'+scenario.replace('_',' '),mark[iscenario]+pat[ir]+col[icomp],"Z [$\Omega$/m]",ax,3,xlab='Frequency [Hz]');

	    # plot dipolar impedances ratio (reference = 4TeV case)
	    if (iscenario>0):
		for icomp,comp in enumerate(listcomp[:2]):
		    # find corresponding term in imp_mod & imp_mod0
		    for kiw,iw in enumerate(imp_mod):
			if test_impedance_wake_comp(iw,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw_comp=kiw;
		    for kiw0,iw0 in enumerate(imp_mod0):
			if test_impedance_wake_comp(iw0,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw0_comp=kiw0;
		    # plot
        	    plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,0]+1j*imp_mod[kiw_comp].func[:,1])/np.abs(imp_mod0[kiw0_comp].func[:,0]+1j*imp_mod0[kiw0_comp].func[:,1]),comp+', post LS1,'+scenario.replace('_',' '),pat[icomp]+col[iscenario],"Imp. ratio w.r.t total coll. imp.",axratio,3,xlab='Frequency [Hz]');


    axratio.set_xlim([1e3,1e11])
    if flagsave: end_figure(fig,ax,save=root_result+'/plot_imp_'+machine+model+suffix)
    else: end_figure(fig,ax);
    if flagsave: end_figure(figratio,axratio,save=root_result+'/plot_imp_ratio_'+machine+model+suffix)
    else: end_figure(figratio,axratio);

    if not(flagsave): pylab.show();
