#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from string import *
import numpy as np
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_conv_testTCTPmode import TCTP_modes
from LHC_param import LHC_param


def TCTP_modes_wake(wakefile,betaavx,betaavy,nmodes=1,scenario=0):

    # compute TCTP low freq. modes wake function and add to the wake in wakefile
    # betaavx and betaavy are the average betas taken for the beta function weighting (usually circ/(2*pi*tune) ),
    # nmodes is the number of resonances taken into account (up to 5), scenario is the
    # TCT settings & beta functions scenario (from Roderik Bruce):
    #	0 = first realistic (small beta*),
    #	1 = second realistic (large beta*),
    #	2 = first pessimistic (8 sigmas, beta*=0.55 m),
    #	3 = second pessimistic (12 sigmas, beta*=1.5 m).
    
    fmodes,Rmodes,Qmodes,halfgap,Rnorm_halfgap,namesTCT,settTCT,modestr=TCTP_modes(scenario=scenario);
    
    wakefilenew=wakefile.replace('.dat','_withTCTPmode_'+modestr+'.dat');
    
    for imode in range(nmodes):
    
        for iTCT,TCTname in enumerate(namesTCT):
	
	    hg=settTCT[iTCT,0];beta=settTCT[iTCT,1];
	    Rmode=np.interp(hg,halfgap,Rnorm_halfgap)*Rmodes[imode]*beta;
	    
	    if (TCTname.startswith('TCTH')): plane='x';
	    elif (TCTname.startswith('TCTV')): plane='y';
	    #print plane
	    
	    betaav=eval('betaav'+plane);
	    
	    if ((imode==0)and(iTCT==0)):
	    	add_resonator_wake(Rmode/betaav,fmodes[imode],Qmodes[imode],wakefile,
			plane=plane,save=wakefilenew);
	    	#print TCTname,", fr=",fmodes[imode],", Rmode=",Rmode/betaav;
	    
	    else:
	    	W=add_resonator_wake(Rmode/betaav,fmodes[imode],Qmodes[imode],wakefilenew,
			plane=plane,save=wakefilenew);
	    	#print TCTname,", fr=",fmodes[imode],", Rmode=",Rmode/betaav;


    return W,modestr;

if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0);
    beta=np.sqrt(1.-1./(gamma**2));
    
    flagsave=1; # 1 to save figure instead of plotting on screen
    os.system("mkdir -p ../../../DELPHI_results/"+machine);

    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    wakeroot='../Impedances/PostLS1/total_wakes/wakeforHDTL_Allthemachine_7TeV_B1_postLS1_'
    
    models=['baseline','relaxed','veryrelaxed'];
    #modelstr=['nominal','tight','relaxed']; # correct names

    for scenario in range(4):

	for imodel,model in enumerate(models):

	    for idip,dipstr in enumerate(['','_dip']):
	    
		wakefile=wakeroot+model+dipstr+'.dat';
	    	wake_old=read_ncol_file(wakefile);
		wake,modestr=TCTP_modes_wake(wakefile,R/Qx,R/Qy,nmodes=5,scenario=scenario);

	    modestrlong=modestr.replace('_',' ').replace('p5','.5').replace('sigmas'," $ \\sigma $").replace('betastar'," $ \\beta^* $ ");
    	    print model,modestrlong
	    modescan=['_old',''];
	    leg=['without TCTP mode','with TCTP mode, '+modestrlong];

    	    # figure, axis and output file name for wake plots
	    figW,axW=init_figure()
	    fileoutW=path_here+'../../../DELPHI_results/'+machine+'/plot_wake_'+machine+'_'+model+'_'+modestr;

	    for imode,mode in enumerate(modescan):

		# plot wake
		W=eval('wake'+mode);
		for iplane,plane in enumerate(['x','y']):
	    	    plot(W[:,0]*1e-9*beta*c,W[:,iplane+1]*1.e15,plane+'dip, '+leg[imode],linetype[iplane]+col[imode],"Wake [V/(C.m)] ",axW,3,xlab='Distance [m]');

	    # finish impedance plot (after the loop scanning with / without TCTP mode)
	    if flagsave: end_figure(figW,axW,save=fileoutW)
	    else: end_figure(figW,axW);

    if not(flagsave): pylab.show();
