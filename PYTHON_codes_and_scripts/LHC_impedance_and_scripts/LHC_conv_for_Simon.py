#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

if len(sys.argv)>1: lxplusbatchDEL=str(sys.argv[1]);
else: lxplusbatchDEL=None;
print lxplusbatchDEL;   

from string import *
import time
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
    E=4e12;
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=E);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/test_Simon';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
       
    kmax=2; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=20; # number of plotted eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition
    zbaseroot='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    model='';beam='1';
    Qpscan=np.arange(-20,31,1);

    dampscan=np.array([0,0.001,0.01]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.5e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    tuneshiftQp=np.zeros((2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    
    # DELPHI loops now

    for iplane,plane in enumerate(['x','y']):
        # select Zxdip or Zydip
        freq,Z=readZ(zbaseroot+'Z'+plane+plane+'dip_Allthemachine_'+Estr+'_B'+beam+model+'.dat')
	freqtot=freq[::10];Ztot=Z[::10,:];

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
	    #for iQp,Qp in enumerate(Qpscan):
	    #	tuneshiftnx=np.zeros((1,len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
	    #	tuneshiftQp[iscenario,iplane,iM,iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
	    #	    	nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
	    #	    	a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
	    #	    	flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
	    #	    	kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
	    #	    	lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
	    #	    	queue='8nh',dire=root_result+'/');
	    #print "Waiting 10 minutes...";
	    #time.sleep(600);

	    # case with Qpscan inside each lxplus job
	    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
	    
	    tuneshiftQp[iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
			kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
			lxplusbatch=lxplusbatchDEL,comment=machine+model+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
			queue='1nw',dire=root_result+'/');


    # now the plots (outside previous loop)
    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	# plots vs Q'
	for iplane,plane in enumerate(['x','y']):

	    for iM,M in enumerate(Mscan):

		for Nb in Nbscanplot:

		    # initialize plots vs Qp
		    figQp=[];axQp=[];
		    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

		    # output file name for plots vs Qp
		    fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

		    strpart=['Re','Im'];
		    for ir,r in enumerate(['real','imag']):

		    	for idamp,damp in enumerate(dampscan):

			    # output file name for data vs Qp
			    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+model+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    #if flagcompute:
			    ts=getattr(tuneshiftQp[iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
			    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
			    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
			    #else:
			    #    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
			    #    Qpscan=s[:,0];ts=s[:,1];

			    sgn=1;sgnstr='';
			    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
			    plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+model+' d='+str(damp),col[idamp],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			    #if (M==1)and(scenario=='_2012'):
			    #	# compare with HEADTAIL
			    #	nsl=500;npr=1000000;nlin=1; # HEADTAIL parameters for comparison
			    #	fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip';
			    #	rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_damper_dip_1b_ntwake20_nkick1_I"+float_to_str(Nb/1e11)+"_qsec0_oct0_drate"+float_to_str(damp)+"_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin);
			    #	#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
			    #	#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
			    #	sufHEADTAIL="_aver_Sussix_most_tau_finer.txt";
			    #	s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
			    #	fact=1;
                            #	if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0
                            #	plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

		    # finish plot vs Qp
		    for ir,r in enumerate(['real','imag']):
			if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
			else: end_figure(figQp[ir],axQp[ir]);


    if not(flagsave): pylab.show();
