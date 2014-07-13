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


def ALBA_imp(circ=268.,gamma=5871,avbetax=6.4917,avbetay=9.1260,dire='ALBA_lumped/',wake_calc=False,lxplusbatch=None):

    # all units are SI.
    # gamma is the relativistic mass factor and circ the total circumference
    # avbetax and avbetay are the average beta functions (where kick will be applied)
    # dire is the directory in ImpedanceWake2D where to put the results
    # wake_calc should be True to compute wake as well
    # lxplusbatch: if None, no use of lxplus batch system
    #            if 'launch' -> launch calculation on lxplus
    #           if 'retrieve' -> retrieve outputs

    imp_mod=[];wake_mod=[];

    beta=np.sqrt(1.-1./(gamma**2))

    zpar=z_param(ztypescan=2,zmin=0.01,zmax=2,nzlog=100,zminrefine=0.5e-5,zmaxrefine=0.03,zsamplin=0.5e-5)
    
    fparBB=freq_param(fmin=100,fmax=1e14,ftypescan=2,nflog=20,fminrefine=1e9,fmaxrefine=1e11,nrefine=10000)
    
    # broad band model parameters     
    Rt_list=np.array([0.203,0.014,0.00828,0.037])*1e6;
    Q_list=np.array([14.56,0.323,33.21,1.71]);
    fr_list=np.array([1.81,7.72,8.90,56.15])*1e9;
    betax=np.ones(4)*avbetax;
    betay=np.ones(4)*avbetay;

    # resonator BB model
    for iRt,Rt in enumerate(Rt_list):   

	imp_mod_BB,wake_mod_BB=imp_model_resonator(Rt,fr_list[iRt],Q_list[iRt],beta=beta,wake_calc=wake_calc,
        	fpar=fparBB,zpar=zpar,listcomp=['Zxdip','Zydip']);

        # add to the model
        add_impedance_wake(imp_mod,imp_mod_BB,betax[iRt]/avbetax,betay[iRt]/avbetay);
        add_impedance_wake(wake_mod,wake_mod_BB,betax[iRt]/avbetax,betay[iRt]/avbetay);



    # resistive-wall model
    if wake_calc: queue='1nd';
    else: queue='8nh'
    
    # parameters
    waketol=1.e11;
    
    fparRW=freq_param(fmin=1e2,fmax=1e14,ftypescan=2,nflog=40,fminrefine=1e11,fmaxrefine=5e12,nrefine=5000,fadded=[])
    layers_list=[[ss304L_layer()],
    	[Cu300K_layer(thickness=0.1e-3),vacuum_layer()],
	[layer(rhoDC=2.5e-7,tau=0.,thickness=1.5e-6),Al_layer(thickness=2e-3),vacuum_layer()],
	[Cu300K_layer(thickness=np.inf)],
	[ss304L_layer()],
	[ss304L_layer()]];
    b_list=np.array([14, 3, 4, 4.25, 12.4, 11.5])*1e-3;
    length_list=np.array([201.4, 4, 8.1, 2.5, 44.28, 4]);
    comment_list=['RRW1','RRW2','RRW3','RRW5','RRW6','RRW7'];
    betax=np.array([7.898,2.229,2.371,2.249,0.829,11.502]);
    betay=np.array([6.201, 1.6705, 1.876, 1.699, 21.53, 6.516]);

    for ilay,lay in enumerate(layers_list):   

	iw_input=impedance_wake_input(machine='ALBA',gamma=gamma,length=1,b=[b_list[ilay]],layers=lay,
		fpar=fparRW,zpar=zpar,waketol=waketol,freqlinbisect=1e11,geometry='flat',comment=comment_list[ilay]);

	imp_mod_RW,wake_mod_RW=imp_model_from_IW2D(iw_input,wake_calc=wake_calc,
	    	flagrm=True,lxplusbatch=lxplusbatch,queue=queue,dire=dire)

	multiply_impedance_wake(imp_mod_RW,length_list[ilay]);
	multiply_impedance_wake(wake_mod_RW,length_list[ilay]);

	# add to the model
	add_impedance_wake(imp_mod,imp_mod_RW,betax[ilay]/avbetax,betay[ilay]/avbetay);
	add_impedance_wake(wake_mod,wake_mod_RW,betax[ilay]/avbetax,betay[ilay]/avbetay);


    return imp_mod,wake_mod


if __name__ == "__main__":

    
    e,m0,c,E0=electron_param();
    
    Nbscan=np.hstack((np.array([0.001]),np.arange(0.05,15.05,0.05)))*1e10;
    #sigmaz_scan=np.array([0.0058,0.0067,0.0067,0.0072,0.0072,0.0075,0.0078,0.0081,0.0083,0.0085, 0.0085]);
    #Qs_scan=np.array([0.0078,0.0067,0.0067,0.0063,0.0063,0.006,0.0058,0.0056,0.0054,0.0053, 0.0053]);
    Qpscan=np.array([0]);Mscan=np.array([1]);
    sigmaz=4.6e-3;Qs=0.00845;

    gamma=5871.;beta=np.sqrt(1.-1./(gamma**2));E=E0*gamma/e;
    Qx=18.15;Qxfrac=Qx-np.floor(Qx);Qy=8.36;Qyfrac=Qy-np.floor(Qy);
    alphap=8.72e-4;circ=268.;R=circ/(2.*np.pi);h=448;eta=alphap-1./gamma**2
    f0=c*beta/circ; # rev. frequency
    omega0=2.*np.pi*f0;omegas=omega0*Qs;
    dphase=0;damp=0;
    machine='ALBA';Estr=float_to_str(round(E/1e8)/10.)+'GeV';
    taub=4.*sigmaz/(beta*c);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used -> we have to use these values in the normalization of the lumped impedance, to be consistent with DELPHI
    print avbetax,avbetay;
    
    kmax=5;kmaxplot=1000;
    root_result=path_here+'../../../DELPHI_results/'+machine;
    os.system("mkdir -p "+root_result);

    wake_calc=True;
    flagsave=1;
    
    imp_mod,wake_mod=ALBA_imp(circ=268.,gamma=5871,avbetax=avbetax,avbetay=avbetay,wake_calc=wake_calc,lxplusbatch=lxplusbatchImp);
    
    if (lxplusbatchImp=='retrieve'):
	plot_compare_imp_model([imp_mod],[''],listcomp=['Zxdip','Zydip','Zxquad','Zyquad'],
	    saveimp=root_result+'/plot_imp_tot',xlim=[1e8,1e11],ylim=[1,1e6]);
	#plot_compare_imp_model([imp_mod],[''],listcomp=['Zxdip','Zydip'],
	#	    saveimp='',xlim=[1e8,1e11],ylim=[1,1e6]);
	#pylab.show()
	
	if wake_calc: 
	    plot_compare_imp_model([wake_mod],[''],listcomp=['Wxdip','Wydip','Wxquad','Wyquad'],
		saveimp=root_result+'/plot_wake_tot',xlim=[0.5e-6,100],ylim=[1e10,1e18]);

	    write_wake_HEADTAIL(wake_mod,'total_lumped_wakeHD.dat',beta=1,ncomp=4)
	    

    	tuneshiftQp=np.zeros((1,len(Mscan),len(Qpscan),1,len(Nbscan),1,1,kmaxplot),dtype=complex);
    	tuneshiftQpm0=np.zeros((1,len(Mscan),len(Qpscan),1,len(Nbscan),1,1),dtype=complex);
	
	#for iNb,Nb in enumerate(Nbscan):

	    #Qs=Qs_scan[iNb];omegas=omega0*Qs;
	    #sigmaz=sigmaz_scan[iNb];
	    #taub=4.*sigmaz/(beta*c);
	    
   	    #print "Nb=",Nb/1e10,"e10, Qs=",Qs,", taub=",taub*1e9,"ns";
	
	g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

	# DELPHI computation
	tuneshiftQp,tuneshiftQpm0=DELPHI_wrapper(imp_mod,
	    Mscan,Qpscan,[damp],Nbscan,[omegas],[dphase],omega0,Qx,Qy,gamma,eta,a,b,
	    taub,g,planes=['y'],nevery=1,particle='electron',flagnorm=0,
	    flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
	    abseps=1.e-4,flagm0=True,lxplusbatch=lxplusbatchDEL,
	    comment=machine+'_'+Estr+'_'+float_to_str(taub*1e9)+'ns',
	    queue='2nd',dire=root_result+'/',flagQpscan_outside=False);
	    

	# now the plots
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['y']):

		for iM,M in enumerate(Mscan):
	    
		    # plots vs Nb, and TMCI plots
		    for Qp in Qpscan:

			# initialize plots vs Nb
			figNb=[];axNb=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figNb.append(fig);axNb.append(ax);
			figNbm0=[];axNbm0=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figNbm0.append(fig);axNbm0.append(ax);
			# initialize TMCI plots
			figTMCI=[];axTMCI=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figTMCI.append(fig);axTMCI.append(ax);

			pat=['.',''];

			# output file name for plots vs Q'
			fileoutplotNb=root_result+'/plot_vs_Nb_'+machine+'_'+Estr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged_'+plane;
			fileoutplotNbm0=root_result+'/plot_vs_Nb_m0_'+machine+'_'+Estr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged_'+plane;
			fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+Estr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged_'+plane;

			strpart=['Re','Im'];ylim=([-5,3],[-0.01,1]);
			for ir,r in enumerate(['real','imag']):

			    fileoutdataNb=root_result+'/data_vs_Nb_'+machine+'_'+Estr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged_'+plane;
			    fileoutdataNbm0=root_result+'/data_vs_Nb_m0_'+machine+'_'+Estr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged_'+plane;

			    ts_most=getattr(tuneshiftQp[iplane,iM,pylab.mlab.find(Qpscan==Qp),0,:,0,0,0],r);
			    data=np.hstack((Nbscan.reshape((-1,1)),ts_most.reshape((-1,1))));
			    write_ncol_file(fileoutdataNb+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_most_unstable")

			    tsm0=getattr(tuneshiftQpm0[iplane,iM,pylab.mlab.find(Qpscan==Qp),0,:,0,0],r);
			    data=np.hstack((Nbscan.reshape((-1,1)),tsm0.reshape((-1,1))));
			    write_ncol_file(fileoutdataNbm0+'_'+r+'.dat',data,header="Nb\t"+strpart[ir]+"_tuneshift_m0")

			    sgn=1;sgnstr='';
			    #if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
			    ylab="$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ";
			    if (ir==1): sgn=-omega0;ylab="Growth rate [s $ ^{-1} $ ] "; # invert sign of imaginary part & conversion to growth rate

			    plot(Nbscan[~np.isnan(np.squeeze(ts_most))],np.squeeze(sgn*ts_most)[~np.isnan(np.squeeze(ts_most))],'','b',ylab,axNb[ir],0,xlab="Intensity [p+/bunch]");
			    plot(Nbscan[~np.isnan(np.squeeze(tsm0))],np.squeeze(sgn*tsm0)[~np.isnan(np.squeeze(ts_most))],'','b',ylab+" (mode 0)",axNbm0[ir],0,xlab="Intensity [p+/bunch]");

			    # set plot axes
			    #axNb[ir].set_xlim([0,5e10]);
			    maxy=np.ceil(np.max(np.abs(sgn*ts_most))*1e5)/1e5;
			    if ir==0: axNb[ir].set_ylim([-maxy,0]);
			    else: ylim[ir][1]=np.ceil(maxy/omegas*5)/5.;axNb[ir].set_ylim([0,maxy]);

			    # TMCI plot
			    ts=tuneshiftQp[iplane,iM,pylab.mlab.find(Qpscan==Qp),0,:,0,0,:];
			    plot_TMCI(Nbscan,np.squeeze(ts/Qs),axTMCI[ir],part=r,leg='',patcol=pat[ir]+'b',xlab='Nb [p+/b]',
				title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.),ms=1,ylim=ylim[ir]);


			# finish plots vs Nb and TMCI plots
			for ir,r in enumerate(['real','imag']):
			    end_figure(figNb[ir],axNb[ir],save=flagsave*(fileoutplotNb+'_'+r))
			    end_figure(figNbm0[ir],axNbm0[ir],save=flagsave*(fileoutplotNbm0+'_'+r))
			    end_figure(figTMCI[ir],axTMCI[ir],save=flagsave*(fileoutplotTMCI+'_'+r));

