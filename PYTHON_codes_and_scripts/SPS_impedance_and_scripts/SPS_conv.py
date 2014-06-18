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
import numpy as np
import pickle as pick
from copy import deepcopy
import pylab,os,re
sys.path.append("../PYTHON/")
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *


def SPS_param(E0,E=26e9,optics='Q26'):

    # E is the energy in eV

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='SPS';
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=6911; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    
    if optics=='Q26':
	Qx=26.13;
	Qy=26.18;
	alphap=1.9181e-3; # momentum compaction factor
	if (E==26e9): Qs=0.00725;Estr='26GeV';taub=2.6e-9; # full length in s. This corresponds to 3MV
	elif (E==450e9): Qs=0.00467;Estr='450GeV';taub=1.5e-9; # full length in s
	else:
    	    print "SPS energy not recognized; Q26 injection parameters taken";
	    Qs=7.25e-3;Estr=float_to_str(E/1e9)+'GeV';taub=2.6e-9; # full length in s

    elif optics=='Q20':
	Qx=20.13;
	Qy=20.18;
	alphap=3.1e-3; # momentum compaction factor
	if (E==26e9): Qs=0.01513;Estr='26GeV';taub=2.9e-9; # full length in s. This corresponds to 3MV
	elif (E==450e9): Qs=0.00595;Estr='450GeV';taub=1.6e-9; # full length in s
	else:
    	    print "SPS energy not recognized; Q20 injection parameters taken";
	    Qs=0.01513;Estr=float_to_str(E/1e9)+'GeV';taub=2.9e-9; # full length in s

    else: print 'SPS_param: only Q20 & Q26 optics implemented';sys.exit();
    
    sigmaz=taub*beta*c/4.; # RMS bunch length (m)
    Qxfrac=Qx-np.floor(Qx);
    Qyfrac=Qy-np.floor(Qy);
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    dphase=0.; # additional damper phase

    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;

    
if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=26e9);

    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result='../DELPHI_results/'+machine;
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    flagplotQp=1; # 1 to do the plots vs Qp
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=50; # number of plotted eigenvalues (for TMCI)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];
    

    # scan definition
    scenarioscan=['_kickers_wall_BPMs_flanges_cavities_Q26_modifNico'];
    legscen=['full model'];
    csiscan=np.arange(-0.5,0.6,0.1);
    csiplotscan=np.array([0]);

    dampscan=np.array([0,0.05,0.1]); # damper gain scan
    Nbscan=np.arange(1.e10,5.1e11,5.e9); # intensity scan
    Nbscanplot=np.array([1.2e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    tuneshiftQp=np.zeros((len(scenarioscan),2,len(Mscan),len(csiscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    imp_mod_list=[]; # complete list of impedance scenarios
    wake_mod_list=[];# complete list of wake scenarios
    
    for iscenario,scenario in enumerate(scenarioscan):
    
	
	# read model
	imp_mod=imp_model_from_file('../Impedances/SPS/From_Carlo_impSPS/Zxdip'+scenario+'.txt','Zxdip');
	imp_mod1=imp_model_from_file('../Impedances/SPS/From_Carlo_impSPS/Zydip'+scenario+'.txt','Zydip');imp_mod.append(imp_mod1[0]);
	imp_mod1=imp_model_from_file('../Impedances/SPS/From_Carlo_impSPS/Zxquad'+scenario+'.txt','Zxquad');imp_mod.append(imp_mod1[0]);
	imp_mod1=imp_model_from_file('../Impedances/SPS/From_Carlo_impSPS/Zyquad'+scenario+'.txt','Zyquad');imp_mod.append(imp_mod1[0]);
	
	imp_mod_list.append(imp_mod);
	
	# dump into a file
	filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	pick.dump(imp_mod,filemodel);
	filemodel.close();

	# plot and compare all scenarios with 2012 impedance
	#maxratio=plot_compare_imp_model(imp_mod_list,legscen,listcomp=['Zxdip','Zydip'],
	#    saveimp=root_result+'/plot_imp_'+machine+"postLS1_scenarios",
	#    saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_scenarios");

	# plot and compare all scenarios with 2012 wake
	#maxratio_w=plot_compare_imp_model(wake_mod_list,legscen,listcomp=['Wxdip','Wydip'],
	#    savewake=root_result+'/plot_wake_'+machine+"postLS1_scenarios",
	#    saveratio=root_result+'/plot_wake_ratio_'+machine+"postLS1_scenarios",
	#    xlim=[1e-5,1e6],ylim=[1e8,1e19]);


    	# DELPHI loops now
	for iplane,plane in enumerate(['x','y']):
	    # select Zxdip or Zydip
	    for iw in imp_mod:
		if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);

	    Qpscan=csiscan*eval('Q'+plane);

	    for iM,M in enumerate(Mscan):

		# normalization factor for damper
		dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
    		    flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
		dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

		flag_trapz=0; # by default no trapz method

        	if (M==1): nxscan=np.array([0]);flag_trapz=1;
		else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,100),np.arange(M/2-10,M/2+11),
	    	    np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

		for iQp,Qp in enumerate(Qpscan):
		
		    tuneshiftQp[iscenario,iplane,iM,iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
			    lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_csi'+float_to_str(round(100.*csiscan[iQp])/100.)+strnorm[flagnorm]+'_'+plane,
			    queue='1nd',dire=root_result+'/');


    # now the plots (outside loop on scenarios)
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	for iplane,plane in enumerate(['x','y']):

	    Qpscan=csiscan*eval('Q'+plane);

	    for iM,M in enumerate(Mscan):

		for idamp,damp in enumerate(dampscan):

	            if flagplotQp:
			# plots vs Q'
			for Nb in Nbscanplot:

			    # initialize plots vs Qp
			    figQp=[];axQp=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

			    # output file name for plots vs Qp
			    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    strpart=['Re','Im'];
			    for ir,r in enumerate(['real','imag']):

				for iscenario,scenario in enumerate(scenarioscan):

				    # output file name for data vs Qp
				    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    #if flagcompute:
				    ts=getattr(tuneshiftQp[iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    data=np.hstack((csiscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="csi\t"+strpart[ir]+"_tuneshift")
				    #else:
				    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
				    #    Qpscan=s[:,0];ts=s[:,1];

				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(csiscan,np.squeeze(sgn*ts),'DELPHI, '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

				    if (M==1):
					# compare with HEADTAIL
					nsl=50;npr=100000;nlin=1; # HEADTAIL parameters for comparison
					fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip_quad';
					rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/SPS/SPS_1b_ntwake10_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_I"+float_to_str(Nb*1.08333333333/1e11)+"_pre1_drate"+float_to_str(damp)+"_flagdamp0";
					#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
					#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
					sufHEADTAIL="_Sussix_aver_most_tau.txt";
					s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
					fact=1;
                                	if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0
                                	plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
				if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
				else: end_figure(figQp[ir],axQp[ir]);



		    # TMCI plots
		    for csi in csiplotscan:

			icsi=pylab.mlab.find(abs(csiscan-csi)<=1e-8);icsi=icsi[0];print icsi;
			fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_csi'+float_to_str(round(1000.*csi)/1000.)+'_converged'+strnorm[flagnorm]+'_'+plane;
			patcol=['.b','b'];
			ylim=([-5,3],[-0.3,0.3]);
			
		    	for ir,r in enumerate(['real','imag']):

			    fig,ax=init_figure();

			    for iscenario,scenario in enumerate(scenarioscan):

				ts=tuneshiftQp[iscenario,iplane,iM,icsi,idamp,:,0,0,:];

				plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg=legscen[iscenario],patcol=patcol[ir],xlab='Nb [p+/b]',
				    title=machine+r", $ \xi = $ "+str(round(1000*csi)/1000.)+', '+str(1./damp)+' turns damping time',ms=1,ylim=ylim[ir]);

			    if flagsave: end_figure(fig,ax,save=fileoutplotTMCI+'_'+r,fontsize=25);
			    else: end_figure(fig,ax,fontsize=25);

		
		# imag. part of TMCI plot for several damping rates, for csi=0, all on same plot
		for csi in csiplotscan:
		
		    iscenario=0;icsi=pylab.mlab.find(abs(csiscan-csi)<=1e-8);icsi=icsi[0];
		
		    for ir,r in enumerate(['real','imag']):

			fig,ax=init_figure();

			fileoutplotTMCI2=root_result+'/plot_TMCI_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_dscan_csi'+float_to_str(round(1000.*csi)/1000.)+'_converged'+strnorm[flagnorm]+'_'+plane;
			pat=['.',''];
			ylim=([-5,3],[-0.3,0.01]);

			for idamp,damp in enumerate(dampscan):

			    ts=tuneshiftQp[iscenario,iplane,iM,icsi,idamp,:,0,0,:];

			    plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg=str(1./damp)+' turns damping time',patcol=pat[ir]+col[idamp],xlab='Nb [p+/b]',
				title=machine+r", $ \xi = $ "+str(round(1000*csi)/1000.),ms=1,ylim=ylim[ir]);

			if flagsave: end_figure(fig,ax,save=fileoutplotTMCI2+'_'+r,fontsize=25);
			else: end_figure(fig,ax,fontsize=25);


    if not(flagsave): pylab.show();
