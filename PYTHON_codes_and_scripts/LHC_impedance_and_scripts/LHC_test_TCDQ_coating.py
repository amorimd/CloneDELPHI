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
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure,cmap
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

    if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
    elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
    else: lxplusbatchImp=None;lxplusbatchDEL=None;
    print lxplusbatchImp,lxplusbatchDEL;   

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=10; # number of converged eigenvalues (kmax most unstable ones are converged)
    root_result=path_here+'../../../DELPHI_results/'+machine+'/TCDQ';
    os.system("mkdir -p "+root_result);
    suffix='';#suffix='_only_TCSG_IR7' # suffix for output files 
    
    model='_TCDQ';
    name='TCDQ_CFC';material='CFC';halfgap=0.0041588; # 7TeV nominal settings (from R. Bruce)
    angle=0;length=10.4;

    # scan definitions
    coatingscan=1e-6*np.array([0,1,2,3,5,8,10,15,20,50,100]);
    Qpscan=np.arange(-10,21,1);
    dampscan=np.array([0,0.02]); # damper gain scan
    #Nbscan=np.arange(1.e10,5.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.5e11]); # intensity scan
    Mscan=np.array([1,1782]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
        
    figratio,axratio=init_figure(); # figure for impedance ratio
    fig=[];ax=[];
    for i in range(3): figu,axi=init_figure();fig.append(figu);ax.append(axi);

    tuneshiftQp=np.zeros((len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(coatingscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan)),dtype=complex);
    powloss=np.zeros((len(coatingscan),len(Mscan),len(Nbscan)));
   
    for icoat,coat in enumerate(coatingscan):
    
	# compute model for TCDQ
	strcoat='_coatCu_'+str(1e6*coat)+'mum';
	
	if (icoat==0):
	    imp_mod,wake_mod=LHC_singlecoll_iw_model(name,material,halfgap,angle,
	    	gamma,length,wake_calc=False,coatingmat=None,coatingthickness=0,
		lxplusbatch=lxplusbatchImp,comment='_coatCu_0mum',dire='Coll'+model+'/')
	else:
	    imp_mod,wake_mod=LHC_singlecoll_iw_model(name,material,halfgap,angle,
	    	gamma,length,wake_calc=False,coatingmat='Cu300K',coatingthickness=coat,
		lxplusbatch=lxplusbatchImp,comment='_coatCu_'+str(1e6*coat)+'mum',dire='Coll'+model+'/')
	
	
	if (icoat==0): imp_mod0=deepcopy(imp_mod);wake_mod0=deepcopy(wake_mod);
    
	if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):

	    # dump into a file
	    filemodel=open('impedances'+model+strcoat+'.txt','w');
	    pick.dump(imp_mod,filemodel);
	    filemodel.close();

	    #compare_vs_zbase()

	    # first plot all components
	    listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad'];
	    units=["","/m","/m","/m","/m","/m","/m"];
	    # corresponding a,b,c,d and planes in imp_mod
	    lista=[0,1,0,0,0,0,0];listb=[0,0,1,0,0,1,0];listc=[0,0,0,1,0,0,0];listd=[0,0,0,0,1,0,1];
	    listplane=['z','x','y','x','y','x','x'];

            #fig,ax=init_figure();
	    if (np.mod(icoat,3)==1):
		for icomp,comp in enumerate(listcomp[:3]):
		    # find corresponding term in imp_mod
		    for kiw,iw in enumerate(imp_mod):
			if test_impedance_wake_comp(iw,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw_comp=kiw;
		    for ir,r in enumerate(['Re','Im']):
			#print comp,iw_comp.plane,iw_comp.a,iw_comp.b;
        		plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,ir]),r+'('+comp+')'+strcoat.replace('_',' '),pat[ir],"Z [$ \Omega $"+units[icomp]+"]",ax[icomp],3,xlab='Frequency [Hz]');

    	    #if flagsave: end_figure(fig,ax,save=root_result+'/plot_imp_'+machine+model+suffix)
	    #else: end_figure(fig,ax);
	
	    # plot dipolar impedances ratio (reference = 4TeV case)
	    if (icoat>0)and(np.mod(icoat,3)==1):
		for icomp,comp in enumerate(listcomp[:3]):
		    # find corresponding term in imp_mod & imp_mod0
		    for kiw,iw in enumerate(imp_mod):
			if test_impedance_wake_comp(iw,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw_comp=kiw;
		    for kiw0,iw0 in enumerate(imp_mod0):
			if test_impedance_wake_comp(iw0,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw0_comp=kiw0;
		    # plot
        	    plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,0]+1j*imp_mod[kiw_comp].func[:,1])/np.abs(imp_mod0[kiw0_comp].func[:,0]+1j*imp_mod0[kiw0_comp].func[:,1]),comp+', post LS1,'+strcoat.replace('_',' '),mark[icomp],"Imp. ratio w.r.t. no coating",axratio,1,xlab='Frequency [Hz]');

	    # compute power loss
	    for iM,M in enumerate(Mscan):

		for iNb,Nb in enumerate(Nbscan):

		    powloss[icoat,iM,iNb]=power_loss(imp_mod,sigmaz,gamma,Nb,M,2*np.pi*R,powerspectrum='gaussian');


	    # DELPHI computation
	    for iplane,plane in enumerate(['x','y']):
	        # select Zxdip or Zydip
		for iw in imp_mod:
		    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);

		for iM,M in enumerate(Mscan):

		    flag_trapz=0; # by default no trapz method for DELPHI

        	    if (M==1): nxscan=np.array([0]);flag_trapz=1;
		    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,100),np.arange(M/2-10,M/2+11),
	    		np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

		    tuneshiftQp[icoat,iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			    nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
			    kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1.e-3,
			    lxplusbatch=lxplusbatchDEL,comment=machine+strcoat+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
			    dire=root_result+'/');

    		    if (lxplusbatchDEL.startswith('retrieve'))or(lxplusbatchDEL==None):
			# extract tune shift of m=0 mode
    			#for iQp,Qp in enumerate(Qpscan):
			for idamp,damp in enumerate(dampscan):
			    for iNb,Nb in enumerate(Nbscan):
				iQp=pylab.mlab.find(Qpscan==0);
				print coat,", M=",M,", d=",damp,", Nb=",Nb
				tuneshiftm0Qp[icoat,iplane,iM,iQp,idamp,iNb]=extract_between_bounds(np.squeeze(tuneshiftQp[icoat,iplane,iM,iQp,idamp,iNb,0,0,:]),-Qs/2,Qs/2);


    # finish impedance plots
    for icomp,comp in enumerate(listcomp[:3]):
	ax[icomp].set_xlim([1e3,1e11]);
	if flagsave: end_figure(fig[icomp],ax[icomp],save=root_result+'/plot_imp_'+comp+'_'+machine+model+suffix,fontsize=20)
	else: end_figure(fig[icomp],ax[icomp]);
    axratio.set_xlim([1e3,1e11]);axratio.set_ylim([0,1.5]);
    if flagsave: end_figure(figratio,axratio,save=root_result+'/plot_imp_ratio_'+machine+model+suffix,fontsize=20)
    else: end_figure(figratio,axratio);


    # now the plots (outside loop on coating)
    if (lxplusbatchDEL.startswith('retrieve'))or(lxplusbatchDEL==None):

	for iplane,plane in enumerate(['x','y']):

	    for iM,M in enumerate(Mscan):

		for idamp,damp in enumerate(dampscan):

		    for iNb,Nb in enumerate(Nbscan):

			# initialize plots vs Qp
			figQp=[];axQp=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

			# output file name for plots vs Qp
			fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			strpart=['Re','Im'];
			for ir,r in enumerate(['real','imag']):

			    for icoat,coat in enumerate(coatingscan):

				# output file name for data vs Qp
				strcoat='_coatCu_'+str(1e6*coat)+'mum';
				fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strcoat+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

			    #    if flagcompute:
				ts=getattr(tuneshiftQp[icoat,iplane,iM,:,idamp,iNb,0,0,0],r);
				data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
				write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				#else:
			    #	s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
			    #	Qpscan=s[:,0];ts=s[:,1];

				sgn=1;sgnstr='';
				if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				plot(Qpscan,sgn*ts,"CFC with Cu coating "+str(coat*1e6)+" $ \mu $ m",'',"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			# finish plot vs Qp
			for ir,r in enumerate(['real','imag']):
			    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r,fontsize=20)
			    else: end_figure(figQp[ir],axQp[ir]);

	
    if (lxplusbatchDEL.startswith('retrieve'))or(lxplusbatchDEL==None):

	for iM,M in enumerate(Mscan):

	    for iNb,Nb in enumerate(Nbscan):
	    
		# plot power loss vs coating thickness
		fileoutplotpow=root_result+'/plot_power_loss_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_Nb'+float_to_str(Nb/1.e11)+'e11';

		figpow,axpow=init_figure();
		plot(coatingscan*1e6,np.abs(powloss[:,iM,iNb]),"",'-b',"Power loss [W]",axpow,3,xlab="Cu thickness [ $ \mu $ m]");
		if flagsave: end_figure(figpow,axpow,save=fileoutplotpow)
		else: end_figure(figpow,axpow);

		for iplane,plane in enumerate(['x','y']):

		    for idamp,damp in enumerate(dampscan):

			# plot tuneshift of m=0 mode at Q'=0 vs coating thickness

			# output file name
			fileoutplottu=root_result+'/plot_tuneshift_vs_coat_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

		    	figtu,axtu=init_figure();
			#print np.squeeze(np.real(tuneshiftm0Qp[pylab.mlab.find(Qpscan==0),:,iplane,iM,idamp,iNb]))
			iQp=pylab.mlab.find(Qpscan==0);#print np.real(tuneshiftm0Qp[:,iplane,iM,iQp,idamp,iNb]);
			plot(coatingscan*1e6,np.real(tuneshiftm0Qp[:,iplane,iM,iQp,idamp,iNb]),"",'-b',"Real tuneshift of $ m=0 $ mode",axtu,1,xlab="Cu thickness [ $ \mu $ m]");
			if flagsave: end_figure(figtu,axtu,save=fileoutplottu)
			else: end_figure(figtu,axtu);


    if not(flagsave): pylab.show();
