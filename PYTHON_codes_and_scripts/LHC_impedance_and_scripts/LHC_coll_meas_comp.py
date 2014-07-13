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

    e,m0,c,E0=proton_param();
    E=4e12;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);
    avbetax=R/Qx;avbetay=R/Qy;

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_coll_meas_comp';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=True; # to plot impedances
    nevery=1; # downsampling of the impedance (take less points than in the full model)

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=60; # number of kept and plotted eigenvalues (in TMCI plot)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition (including measured parameters)
    scenarioscan=np.array(['_04h24m00_TCSGclosed','_04h28m00_TCSGopened','_05h12m00_TCPclosed',
    	'_05h15m00_TCPopened','_04h26m00_TCSGclosed','_04h31m00_TCSGopened','_05h12m30_TCPclosed',
	'_05h15m00_TCPopened']);
    tickscan=np.array(['TCSG closed','TCSG opened','TCP closed','TCP opened',
    	'TCSG closed','TCSG opened','TCP closed','TCP opened']);
    beamscan=np.array([1,1,1,1,2,2,2,2]);
    Qpxscan=np.array([5., 5., 2., 2., 1.8, 1.8, 1.8, 1.8]);
    Qpyscan=np.array([2., 2., 1.5, 1.5, 1., 1., 1., 1.]);
    Nbscan=np.array([88., 88., 82., 82., 85., 85., 76., 76.])*1e9; # intensity scan
    
    damp=0;M=1;
    squeeze='0p6m_3m_0p6m_3m';
    settings_filename_coll_root=path_here+'Coll_settings/coll_settings_';
   
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    zbasesuf=np.array(['_AllColl_4TeV_B1_20120624_TCSGclosed.dat','_AllColl_4TeV_B1_20120624_TCSGopened.dat',
    	'_AllColl_4TeV_B1_20120624_TCP_closed.dat','_AllColl_4TeV_B1_20120624_TCP_opened.dat',
	'_AllColl_4TeV_B2_20120624_TCSG_closed.dat','_AllColl_4TeV_B2_20120624_TCSG_opened.dat',
    	'_AllColl_4TeV_B2_20120624_TCP_closed.dat','_AllColl_4TeV_B2_20120624_TCP_opened.dat'])

    # for HEADTAIL files
    sufHEADTAIL="_BS_prt.dat1to8192.sussix_dQ";
    rootHEADTAIL="/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/HDTL/output/LHC_I";
    HEADTAILscan=['088_AllColl_4TeV_B1_20120624_TCSGclosed_xiH5V2_epsx2_epsy1p4',
    	'088_AllColl_4TeV_B1_20120624_TCSGopened_xiH5V2_epsx2_epsy1p4',
    	'082_AllColl_4TeV_B1_20120624_TCP_closed_xiH2V1p5',
	'082_AllColl_4TeV_B1_20120624_TCP_opened_xiH2V1p5',
	'085_AllColl_4TeV_B2_20120624_TCSG_closed_xiH1p8V1',
	'085_AllColl_4TeV_B2_20120624_TCSG_opened_xiH1p8V1',
    	'076_AllColl_4TeV_B2_20120624_TCP_closed_xiH1p8V1',
	'076_AllColl_4TeV_B2_20120624_TCP_opened_xiH1p8V1'];

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((len(scenarioscan),2,1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(scenarioscan),2,1,1),dtype=complex);
    tuneshiftQp_RW=np.zeros((len(scenarioscan),2,1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp_RW=np.zeros((len(scenarioscan),2,1,1),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    for iscenario,scenario in enumerate(scenarioscan):
    	
	print "scenario: ",scenario;
	beam=str(beamscan[iscenario]);
	
	param_filename_coll=path_here+'Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
	settings_filename_coll=settings_filename_coll_root+'B'+beam+'_4000GeV_20120624'+scenario+'.txt';
	beta_filename_coll=param_filename_coll;
	
	# compute imp. model for collimators (RW + geom)
    	imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,
		avbetax,avbetay,param_filename_coll,settings_filename_coll,beta_filename_coll,
		wake_calc=wake_calc,ftypescan=0,nflog=100,zpar=z_param(),namesref=None,
		TDIcoating='preLS1',BPM=False,fcutoffBB=50e9,lxplusbatch=lxplusbatchImp,
		comment='',dire='Coll_2012_v2/');
	
	# add up
	imp_mod=[];wake_mod=[];
	add_impedance_wake(imp_mod,imp_mod_coll_RW,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll_RW,1,1);

	add_impedance_wake(imp_mod,imp_mod_coll_geom,1,1);
	add_impedance_wake(wake_mod,wake_mod_coll_geom,1,1);
	
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    #filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    #pick.dump(imp_mod,filemodel);
	    #filemodel.close();
	    
	    # write Ascii files with each component
	    #write_imp_wake_mod(imp_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    #	listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zyxdip','Zxyquad','Zyxquad','Zxcst','Zycst'],
	    #	dire=root_result+'/')
	    
	    # plot and compare with zbase
	    compare_imp_vs_zbase(imp_mod_coll_RW,root_zbase=zbaseroot,suffix_zbase=zbasesuf[iscenario],
	    	save=root_result+'/plot_imp_RW_vs_zbase_'+machine+scenario);

	    if (wake_calc):
		# write Ascii files with each component
		#write_imp_wake_mod(wake_mod,"_Allthemachine_"+Estr+"_B"+str(beam)+scenario,
	    	#    listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wyxdip','Wxyquad','Wyxquad','Wxcst','Wycst'],
		#    dire=root_result+'/')
	    
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
		# dip only
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'_dip.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=2)		

    

	    # DELPHI scans now

	    # write HEADTAIL input file
	    # first do the longitudinal matching
	    delta_lin,eps0_lin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=True);
	    print "Longitudinal linear matching: delta_p/p=",delta_lin,", eps0=",eps0_lin,", Qs=",Qs;
	    delta_nlin,eps0_nlin=long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=False);
	    print "Longitudinal non-linear matching: delta_p/p=",delta_nlin,", eps0=",eps0_nlin,", Qs=",Qs;
	    # check that Qs is correct
	    if ((Qs-Qs_from_RF_param(V,h,gamma,eta))>1e-8): print " Pb in Qs !!";sys.exit();
	    
	    # write .cfg files (linear and non-linear matching)
	    cfgnamelin=root_result+'/'+machine+"_"+Estr+scenario+'_lin.cfg';
	    write_CFG_HEADTAIL(cfgnamelin,Nb=Nbscan[iscenario],betax=R/Qx,betay=R/Qy,
	    	sigmaz=sigmaz,delta=delta_lin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		gamma=gamma,nturns=20000,Qx=Qx,Qy=Qy,Qpx=Qpxscan[iscenario],Qpy=Qpyscan[iscenario],
		isyn=1,start_turn=19900,end_turn=19910,VRF=V,
		dratex=damp,dratey=damp);
	    
	    cfgnamenlin=root_result+'/'+machine+"_"+Estr+scenario+'_nlin.cfg';
	    write_CFG_HEADTAIL(cfgnamenlin,Nb=Nbscan[iscenario],betax=R/Qx,betay=R/Qy,
	    	sigmaz=sigmaz,delta=delta_nlin,Qs=Qs,alphap=eta+1./gamma**2,circ=2*np.pi*R,
		gamma=gamma,nturns=20000,Qx=Qx,Qy=Qy,Qpx=Qpxscan[iscenario],Qpy=Qpyscan[iscenario],
		isyn=4,start_turn=19900,end_turn=19910,VRF=V,
		dratex=damp,dratey=damp);
	    
	    # DELPHI computation
	    for iplane,plane in enumerate(['x','y']):
    		
		Qp=eval('Qp'+plane+'scan[iscenario]');print scenario, plane, Qp;
		Nb=Nbscan[iscenario];
	    
		# total impedance

		# select Zxdip+Zxquad or Zydip+Zyquad
		for iw in imp_mod:
		    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Zd=deepcopy(iw.func[::nevery,:]);freqd=deepcopy(iw.var[::nevery]);
		for iw in imp_mod:
		    if test_impedance_wake_comp(iw,0,0,1-iplane,iplane,plane): Zq=deepcopy(iw.func[::nevery,:]);freqq=deepcopy(iw.var[::nevery]);
		# sum the two
	    	freq=sort_and_delete_duplicates(np.concatenate((freqd,freqq)));
	    	Z=np.zeros((len(freq),2),dtype=float);
	    	for icol in range(2): Z[:,icol]=np.interp(freq,freqd,Zd[:,icol],right=0.)+np.interp(freq,freqq,Zq[:,icol],right=0.);
		
		tuneshiftnx=np.zeros((1,1,1,kmaxplot),dtype=complex);
		
	    	tuneshiftQp[iscenario,iplane,:,:,:],tuneshiftnx,tuneshiftm0Qp[iscenario,iplane,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    [0],[damp],[Nb],[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			    flag_trapz=1,flagdamperimp=0,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5e-2,abseps=1e-4,flagm0=True,
			    lxplusbatch=lxplusbatchDEL,comment=scenario+'_'+plane,
			    queue='8nh',dire=root_result+'/');
	    
		# RW impedance alone

		# select Zxdip+Zxquad or Zydip+Zyquad
		for iw in imp_mod_coll_RW:
		    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): ZdRW=deepcopy(iw.func[::nevery,:]);freqdRW=deepcopy(iw.var[::nevery]);
		for iw in imp_mod_coll_RW:
		    if test_impedance_wake_comp(iw,0,0,1-iplane,iplane,plane): ZqRW=deepcopy(iw.func[::nevery,:]);freqqRW=deepcopy(iw.var[::nevery]);
		# sum the two
	    	freqRW=sort_and_delete_duplicates(np.concatenate((freqdRW,freqqRW)));
	    	ZRW=np.zeros((len(freq),2),dtype=float);
	    	for icol in range(2): ZRW[:,icol]=np.interp(freqRW,freqdRW,ZdRW[:,icol],right=0.)+np.interp(freqRW,freqqRW,ZqRW[:,icol],right=0.);
		
		#fig,ax=init_figure();
		#ax.loglog(freq,Z[:,0],'b',freq,Z[:,1],'--b');
		#ax.loglog(freqRW,ZRW[:,0],'r',freqRW,ZRW[:,1],'--r');
		#end_figure(fig,ax,save=root_result+'/tmp_Z_plot_'+plane+scenario);
		
	    	tuneshiftQp_RW[iscenario,iplane,:,:,:],tuneshiftnx,tuneshiftm0Qp_RW[iscenario,iplane,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    [0],[damp],[Nb],[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,ZRW,freqRW,particle='proton',flagnorm=flagnorm,
			    flag_trapz=1,flagdamperimp=0,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5e-2,abseps=1e-4,flagm0=True,
			    lxplusbatch=lxplusbatchDEL,comment=scenario+'_RW_'+plane,
			    queue='8nh',dire=root_result+'/');
			    
		#print scenario,plane,tuneshiftm0Qp_RW[iscenario,iplane,:,:],tuneshiftm0Qp[iscenario,iplane,:,:]
	    
	    
    # now the plots (outside loop on scenarios)
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	rotationangle=75; # angle for labels
	sizexlab='x-small' # size for xticks
	#sizexlab=10;
	yoff=0.2; # offset for position of axis y=0 (w.r.t. window)

	# plot comparison between HEADTAIL and DELPHI
	fileoutplot_compHEAD_ts=root_result+'/plot_HEADTAIL_vs_DELPHI_tuneshifts_RW';
	fileoutplot_compHEAD_ratio=root_result+'/plot_HEADTAIL_vs_DELPHI_ratio_RW';
	tsHEAD=np.zeros((len(scenarioscan),2));
	tick_names=[];

	# read HEADTAIL data
    	for iscenario,scenario in enumerate(scenarioscan):

	    s=read_ncol_file(rootHEADTAIL+HEADTAILscan[iscenario]+sufHEADTAIL);
	    tsHEAD[iscenario,:]=s;
	    tick_names.append('B'+str(beamscan[iscenario])+'H, '+tickscan[iscenario]);
	    tick_names.append('B'+str(beamscan[iscenario])+'V, '+tickscan[iscenario]);

	# reshape tuneshifts in single-line arrays
	tsHEAD=np.squeeze(tsHEAD.reshape((1,-1)));
	tsDELPHI_RW=getattr(tuneshiftm0Qp_RW[:,:,0,0],'real');
	tsDELPHI_RW=np.squeeze(tsDELPHI_RW.reshape((1,-1)));
	
	# plot absolute tune shifts
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.plot(np.arange(len(tsHEAD)),tsHEAD,'xb',label='HEADTAIL',ms=10,lw=3);
	ax.plot(np.arange(len(tsDELPHI_RW)),tsDELPHI_RW,'or',label='DELPHI',ms=10,lw=3);
	pylab.xticks(np.arange(len(tsHEAD)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shift');
	end_figure(fig,ax,save=flagsave*(fileoutplot_compHEAD_ts),fontsize=25)

	# plot bars with discrepancy factor
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.bar(np.arange(len(tsHEAD)),tsHEAD/tsDELPHI_RW,facecolor='b',edgecolor='k',width=0.8);
	pylab.xticks(np.arange(len(tsHEAD)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shifts HEADTAIL / DELPHI');
	end_figure(fig,ax,save=flagsave*(fileoutplot_compHEAD_ratio),fontsize=25)

	
	sizexlab='small' # size for xticks
	
	# plot comparison between DELPHI and measurements
	fileoutplot_comp_ratio_RW=root_result+'/plot_DELPHI_vs_measurements_ratio_RW';
	fileoutplot_comp_ratio=root_result+'/plot_DELPHI_vs_measurements_ratio';
	fileoutplot_comp_ts=root_result+'/plot_DELPHI_vs_measurements_tuneshifts';
	fileoutplot_comp_ratio_total_RW=root_result+'/plot_DELPHI_ratio_total_RW';
	tick_names=['B1H, TCSG IR7','B1V, TCSG IR7','B2H, TCSG IR7','B2V, TCSG IR7',
		'B1H, TCP IR7','B1V, TCP IR7','B2V, TCP IR7'];
	#ind_conv=np.array([0,4,2,6,1,5,7]); # indices of conversion (from DELPHI arrays to measurement ones)
	ind_conv=np.array([0,1,4,5,2,3,7]); # indices of conversion (from DELPHI arrays to measurement ones)

	tsmeas=-np.array([0.3, 0.39, 0.4, 0.43, 0.09, 0.07, 0.1])*1e-3;
	errmeas=np.array([0.02, 0.04, 0.03, 0.07, 0.02, 0.01, 0.01])*1e-3;
	
	# DELPHI tuneshifts from resistive-wall
	tsRW=np.squeeze(getattr(tuneshiftm0Qp_RW[:,:,0,0],'real'));
	tsDELPHI_RW=tsRW[:-1:2,:]-tsRW[1::2,:];
	tsDELPHI_RW=np.squeeze(tsDELPHI_RW.reshape((1,-1)));
	tsDELPHI_RW=tsDELPHI_RW[ind_conv];
	# DELPHI total tuneshifts
	ts=np.squeeze(getattr(tuneshiftm0Qp[:,:,0,0],'real'));
	tsDELPHI=ts[:-1:2,:]-ts[1::2,:];
	tsDELPHI=np.squeeze(tsDELPHI.reshape((1,-1)));
	tsDELPHI=tsDELPHI[ind_conv];
	
	print ts;print tsRW;
	print (tsDELPHI-tsDELPHI_RW)/tsDELPHI;
	fig,ax=init_figure()
	ax.plot(np.arange(len(tsDELPHI)),tsDELPHI,'xb');
	ax.plot(np.arange(len(tsDELPHI_RW)),tsDELPHI_RW,'or');
	
	# test ind_conv
	namesB1=np.array(['B1'+el for el in scenarioscan[:4]]);
	namesB2=np.array(['B2'+el for el in scenarioscan[4:]]);
	names=np.hstack((namesB1,namesB2));
	names=np.hstack((names[:-1:2].reshape((-1,1)),names[1::2].reshape((-1,1))));print names;
	names[:,0]=np.array(['H'+el for el in names[:,0]]);
	names[:,1]=np.array(['V'+el for el in names[:,1]]);print names;
	names=np.squeeze(names.reshape((1,-1)));
	print names[ind_conv];
	print tick_names;
	
	# plot total absolute tuneshifts
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.errorbar(np.arange(len(tsmeas)),tsmeas,yerr=errmeas,fmt='xb',label='Measurements',ms=10,lw=3);
	ax.plot(np.arange(len(tsDELPHI)),tsDELPHI,'or',label='DELPHI',ms=10,lw=3);
	pylab.xticks(np.arange(len(tsmeas)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shift');
	end_figure(fig,ax,save=flagsave*(fileoutplot_comp_ts),fontsize=25)
	
	# plot bars with discrepancy factor (total tuneshifts)
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.bar(np.arange(len(tsmeas)),tsmeas/tsDELPHI,facecolor='b',edgecolor='k',width=0.8);
	pylab.xticks(np.arange(len(tsmeas)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shifts measured / simulated (DELPHI)');
	end_figure(fig,ax,save=flagsave*(fileoutplot_comp_ratio),fontsize=25)

	# plot bars with discrepancy factor (resistive-wall only)
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.bar(np.arange(len(tsmeas)),tsmeas/tsDELPHI_RW,facecolor='b',edgecolor='k',width=0.8);
	pylab.xticks(np.arange(len(tsmeas)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shifts measured / simulated (DELPHI, RW)');
	end_figure(fig,ax,save=flagsave*(fileoutplot_comp_ratio_RW),fontsize=25)

	# plot bars with ratio RW+geom/RW
	fig,ax=init_figure(axes=[0.15,yoff,0.8,0.95-yoff]);
	ax.bar(np.arange(len(tsDELPHI)),tsDELPHI/tsDELPHI_RW,facecolor='b',edgecolor='k',width=0.8);
	pylab.xticks(np.arange(len(tsDELPHI)), tick_names, rotation=rotationangle, size=sizexlab);
	ax.set_ylabel('Tune shifts Total/RW (simulated, DELPHI)');
	end_figure(fig,ax,save=flagsave*(fileoutplot_comp_ratio_total_RW),fontsize=25)


# tune slopes units are 10^-3/(10^11 p+/b)
#tuneslopemeas4000=[0.3/0.88 0.39/0.88 0.4/0.85 0.43/0.85 0.09/0.82 0.07/0.82 0.1/0.76];
#errmeas4000=[0.02/0.88 0.04/0.88 0.03/0.85 0.07/0.85 0.02/0.82 0.01/0.82 0.01/0.76]

#tuneslopeHEAD4000=[0.171/0.88 0.191/0.88 0.174/0.85 0.187/0.85 0.044/0.82 0.062/0.82 0.054/0.76];
#errHEAD4000=[0. 0. 0. 0. 0. 0. 0.]


    if not(flagsave): pylab.show();
