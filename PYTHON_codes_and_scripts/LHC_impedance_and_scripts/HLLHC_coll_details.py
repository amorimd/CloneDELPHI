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
from HLLHC_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine+'/LHC_coll_details';
    
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
    scenario='_HLLHC';optics='baseline';
    param_filename_coll='../Coll_settings/collgaps_HLLHC_'+optics+'_from_Roderik_modifNico_material_names.dat';
    settings_filename_coll=param_filename_coll;
    beta_filename_coll=param_filename_coll;
    squeeze='0p15m_round';
    fminrefine=1.e11;fmaxrefine=5.e12;nrefine=5000;ftypescan=2;nflog=100;
    BPMflag=True
    
    # read coll. file
    namesref,material,thick,angle,length,halfgap,betax,betay=read_coll_files_several_mat(param_filename_coll,settings_filename_coll,beta_filename_coll,namesref=None);

    # scans
    partscan=['_RW','_geom','_tot',''];
    Qpscan=np.array([0]);
    #Qpaver=np.array([14,15,16]);
    #Qpscanplot=np.array([0,15]);
    #iQpaver=select_in_table(Qpaver,Qpscan);# print iQpaver

    dampscan=np.array([0]); # damper gain scan
    #Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.7e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        

    tuneshiftQp=np.zeros((len(partscan),len(namesref)+1,2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftQpm0=np.zeros((len(partscan),len(namesref)+1,2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    kickfactor=np.zeros((len(partscan),len(namesref)+1,2));
    
    # compute total imp. model
    imp_mod,wake_mod=HLLHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
	    settings_filename_coll,TDIcoating='preLS1',dire="../LHC_elements/",
	    commentcoll=scenario,direcoll='Coll'+scenario+'_v2/',
	    lxplusbatch=lxplusbatchImp,BPM=BPMflag,beam=beam,squeeze=squeeze,
	    wake_calc=wake_calc,flagplot=True,root_result=root_result)

    imp_mod_list=[imp_mod];
    
    # loop to construct model for each collimators
    imp_mod_RW=[];imp_mod_geom=[];imp_mod_tot=[];
    imp_mod_RW_list=[];imp_mod_geom_list=[];imp_mod_tot_list=[]; # list of models of each coll.
    for iname,name in enumerate(namesref):
    
	if name.startswith('TDI'):
	    # special case of the TDI
	    imp_RW,wake_RW,imp_geom,wake_geom=LHC_TDI_iw_model_with_geom(name,halfgap[iname],angle[iname],gamma,
		    TDIcoating='preLS1',wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,
		    lxplusbatch=lxplusbatchImp,comment=scenario,dire='Coll'+scenario+'/');

	else:
	    # reorder materials and thicknesses
	    materials=[material[n][iname] for n in range(len(material))];
	    thicks=[thick[n][iname] for n in range(0,len(material)-1)];thicks.append(np.inf);
	
	    imp_RW,wake_RW,imp_geom,wake_geom=LHC_singlecoll_iw_model_with_geom(name,materials,
	    	halfgap[iname],angle[iname],gamma,length[iname],wake_calc=wake_calc,
		fpar=freq_param(ftypescan=ftypescan,nflog=nflog,fminrefine=fminrefine,
		fmaxrefine=fmaxrefine,nrefine=nrefine),BPM=BPMflag,lxplusbatch=lxplusbatchImp,
		comment=scenario+'_'+materials[0]+'_'+float_to_str(round(halfgap[iname]*1e5)/1e2)+'mm',dire='Coll'+scenario+'/');

	imp_tot=[];
	add_impedance_wake(imp_tot,imp_RW,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(imp_tot,imp_geom,betax[iname]/avbetax,betay[iname]/avbetay);
	imp_mod_tot_list.append(imp_tot);
	imp_mod_RW_list.append(imp_RW);
	imp_mod_geom_list.append(imp_geom);
	
	add_impedance_wake(imp_mod_RW,imp_RW,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(imp_mod_geom,imp_geom,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(imp_mod_tot,imp_tot,betax[iname]/avbetax,betay[iname]/avbetay);

    # last element of lists is total model of collimators (RW, geom or sum of the 2)
    imp_mod_RW_list.append(imp_mod_RW);
    imp_mod_geom_list.append(imp_mod_geom);
    imp_mod_tot_list.append(imp_mod_tot);
    
    namesref.append('total');
	
    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

        # DELPHI loops now
	for ipart,part in enumerate(partscan):
	    
	    for imod,imp_mod in enumerate(eval('imp_mod'+part+'_list')):
	    
		if ipart<=2: namestr=namesref[imod];
		else: namestr='total_model';

		for iplane,plane in enumerate(['x','y']):
		    # select Zxdip or Zydip
		    for iw in imp_mod:
			if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);

		    kickfactor[ipart,imod,iplane]=transverse_kick_factor(imp_mod,sigmaz,powerspectrum='gaussian',compname='Z'+plane+'dip');

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

			# to put loop on Qpscan outside the lxplus job
			for iQp,Qp in enumerate(Qpscan):
			    tuneshiftnx=np.zeros((1,len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
			    tuneshiftQp[ipart,imod,iplane,iM,iQp,:,:,:,:,:],tuneshiftnx,tuneshiftQpm0[ipart,imod,iplane,iM,iQp,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
		    		    nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
		    		    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
		    		    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
		    		    kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,flagm0=True,
		    		    lxplusbatch=lxplusbatchDEL,comment=machine+scenario+part+'_'+namestr+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
		    		    queue='8nh',dire=root_result+'/');
			#print "Waiting 20 minutes...";
			#time.sleep(1200);


	# now the plots
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    # plot for each case histograms of the contribution of each collimator to
	    # total tuneshift
	    rotationangle=75; # angle for labels
	    #sizexlab='xx-small' # size for xticks
	    sizexlab=10;
	    
	    for ipart,part in enumerate(partscan[:-1]):

		for iplane,plane in enumerate(['x','y']):

		    # plot kickfactor vs collimators
					
		    # output files
		    suffix=machine+'_'+float_to_str(round(E/1e9))+'GeV'+scenario+part+'_'+plane;
		    fileoutplot=root_result+'/plot_kickfactor_vs_coll_'+suffix;
		    fileoutdata=root_result+'/data_kickfactor_vs_coll_'+suffix+'.dat';

		    # tuneshift and ratio vs total (of coll., and full model)
		    kick=-kickfactor[ipart,:,iplane];
		    kicktot=-kickfactor[ipart,-1,iplane];
		    kickfulltot=-kickfactor[-1,0,iplane];
		    ratio_coll=kick/kicktot;
		    ratio_full=kick/kickfulltot;

		    # plot tuneshift
		    fig,ax=init_figure();					
		    ax.bar(np.arange(len(namesref)-1),kick[:-1],facecolor='b',edgecolor='k',width=0.8);
		    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
		    ax.set_ylabel('Transverse kick factor');
		    end_figure(fig,ax,save=flagsave*(fileoutplot))

		    # plot ratios
		    fig,ax=init_figure();					
		    ax.bar(np.arange(len(namesref)-1),ratio_coll[:-1],facecolor='b',edgecolor='k',width=0.8);
		    ax.set_ylabel('Transverse kick factor / total coll. kick factor');
		    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
		    end_figure(fig,ax,save=flagsave*(fileoutplot+'_ratio_coll'))

		    fig,ax=init_figure();					
		    ax.bar(np.arange(len(namesref)-1),ratio_full[:-1],facecolor='b',edgecolor='k',width=0.8);
		    ax.set_ylabel('Transverse kick factor / total model kick factor');
		    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
		    end_figure(fig,ax,save=flagsave*(fileoutplot+'_ratio_full'))

		    # write data
		    fid=open(fileoutdata,'w');
		    print >>fid, "Coll_name\tKickfactor\tRatio_vs_total_coll\tRatio_vs_full_model";
		    for iname,name in enumerate(namesref): print >>fid, name,"\t",kick[iname],"\t",ratio_coll[iname],"\t",ratio_full[iname];
		    fid.close();

		    for iM,M in enumerate(Mscan):

			for iQp,Qp in enumerate(Qpscan):

			    for idamp,damp in enumerate(dampscan):

				for iNb,Nb in enumerate(Nbscan):
				    
				    partstr=['\Re','\Im'];
			    	    for ir,r in enumerate(['real','imag']):

					m0ind=[',0','']; # additional index needed if not tuneshift of mode 0
					for im0,m0 in enumerate(['','m0']):
					
					    # plot tuneshifts vs collimators

					    # output files
					    suffix=m0+'_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+scenario+part+'_'+str(M)+'b_Qp'+float_to_str(Qp)+'_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+'_'+r;
					    fileoutplot=root_result+'/plot_vs_coll_'+suffix;
					    fileoutdata=root_result+'/data_vs_coll_'+suffix+'.dat';

				            # tuneshift and ratio vs total (of coll., and full model)
				    	    ts=getattr(eval('tuneshiftQp'+m0+'[ipart,:,iplane,iM,iQp,idamp,iNb,0,0'+m0ind[im0]+']'),r);
					    tstot=getattr(eval('tuneshiftQp'+m0+'[ipart,-1,iplane,iM,iQp,idamp,iNb,0,0'+m0ind[im0]+']'),r);
					    tsfulltot=getattr(eval('tuneshiftQp'+m0+'[-1,0,iplane,iM,iQp,idamp,iNb,0,0'+m0ind[im0]+']'),r);
					    ratio_coll=ts/tstot;
					    ratio_full=ts/tsfulltot;

					    # plot tuneshift
					    fig,ax=init_figure();					
					    ax.bar(np.arange(len(namesref)-1),-ts[:-1],facecolor='b',edgecolor='k',width=0.8);
		    			    ax.set_ylabel(" $ -"+partstr[ir]+"(\Delta Q) $ ");
					    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
					    end_figure(fig,ax,save=flagsave*(fileoutplot))

					    # plot ratios
					    fig,ax=init_figure();					
					    ax.bar(np.arange(len(namesref)-1),ratio_coll[:-1],facecolor='b',edgecolor='k',width=0.8);
					    ax.set_ylabel('Tuneshift / total coll. tuneshift');
					    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
					    end_figure(fig,ax,save=flagsave*(fileoutplot+'_ratio_coll'))

					    fig,ax=init_figure();					
					    ax.bar(np.arange(len(namesref)-1),ratio_full[:-1],facecolor='b',edgecolor='k',width=0.8);
		    			    ax.set_ylabel('Tuneshift / total model tuneshift');
					    pylab.xticks(np.arange(len(namesref)-1), namesref[:-1], rotation=rotationangle, size=sizexlab);
					    end_figure(fig,ax,save=flagsave*(fileoutplot+'_ratio_full'))

					    # write data
					    fid=open(fileoutdata,'w');
					    print >>fid, "Coll_name\tTuneshift\tRatio_vs_total_coll\tRatio_vs_full_model";
					    for iname,name in enumerate(namesref): print >>fid, name,"\t",ts[iname],"\t",ratio_coll[iname],"\t",ratio_full[iname];
					    fid.close();



    if not(flagsave): pylab.show();
