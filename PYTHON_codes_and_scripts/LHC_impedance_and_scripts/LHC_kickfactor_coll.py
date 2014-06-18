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
from copy import deepcopy
import pylab,os,re
sys.path.append("../PYTHON/")
from plot_lib import plot,init_figure,end_figure,cmap
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_conv import LHC_param
from LHC_imp import *
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    
    mu0=4e-7*np.pi;Z0=mu0*c;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    alphap=eta+1./(gamma*gamma);

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
    flagSach=(lxplusbatchDEL=='launch'); # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    flagSach=0; # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    actionMOSES='retrieve'; # 'launch' or 'retrieve' (but acannot do 'launch' from lxplus...)
    
    kmax=5; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=20; # number of plotted eigenvalues
    root_result='../DELPHI_results/'+machine+'/coll_kick_factor';
    suffix='';#suffix='_only_TCSG_IR7' # suffix for output files 
    
    model='_test';
    angle=0;length=1.;
    
    # BB model for geometric impedance: parameters
    hgBBlist=np.array([1e-3,3e-3,11.5e-3]); # half-gaps where BB value given
    RtBBlist=np.array([4e5,6.7e4,4.7e3]); # shunt impedances at these half-gaps
    frBB=5e9;Q=1; # cutoff and quality factor (standard LHC BB model)

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # scan definitions
    materialscan=['CFC','W'];rhoscan=np.array([5e-6,5.4e-8]);
    scenarioscan=deepcopy(materialscan);scenarioscan.append('geometric');
    hgapscan=1e-3*np.array([1,2,3,4,5,6,7,8,10,12]);
    Nbscan=np.arange(1.e10,5.05e11,5.e9);
    Nbscanplot=np.array([1.e10,1.e11,2.e11]);
    dampscan=np.array([0]);
    Mscan=np.array([1]);
    Qpscan=np.arange(-10,21);
    Qpplothg=np.array([0,2,5,10]); # scan in Qp for plot vs half-gap (and for TMCI plot)

    # parameters for MOSES and Sacherer
    nNbMOS=min(120,len(Nbscan));lmaxMOS=3;nmaxMOS=3;lmaxSach=3;
    MOSESrootname='/home/nmounet/DFS/Documents/MOSES4W/Results_LHC/MZobov_geom_BB_coll_impedance/LHC'
    
    imp_mod_list={}; # complete list of impedance scenarios
    wake_mod_list={};# complete list of wake scenarios
    kickfactor=np.zeros((len(materialscan),len(hgapscan)));
    kickfactorgeom=np.zeros(len(hgapscan));
   
    for imat,material in enumerate(materialscan):
    
        name='coll_'+material;
	imp_mod_list[material]=[];

	for ihalfgap,halfgap in enumerate(hgapscan):

	    # compute model for collimator
	    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

	    imp_mod,wake_mod=LHC_singlecoll_iw_model(name,material,halfgap,angle,
		gamma,length,wake_calc=False,coatingmat=None,coatingthickness=0,
		lxplusbatch=lxplusbatchImp,comment=material+strhgap,dire='Coll'+model+'/')

	    if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):

	        imp_mod_list[material].append(imp_mod);
		
		# dump into a file
		filemodel=open(root_result+'/impedances'+model+name+strhgap+'.txt','w');
		pick.dump(imp_mod,filemodel);
		filemodel.close();
		
		# compute RW transverse kick factor
		kickfactor[imat,ihalfgap]=transverse_kick_factor(imp_mod,sigmaz,powerspectrum='gaussian');
		kickfactor[imat,ihalfgap]+=transverse_kick_factor(imp_mod,sigmaz,powerspectrum='gaussian',compname='Zxquad');

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):
	
	# compute BB model from M. Zobov value
	imp_mod_list['geometric']=[];
	for ihalfgap,halfgap in enumerate(hgapscan):

	    # compute BB model from M. Zobov values 
	    # first linear interpolation of his curve vs half-gap, in loglog
	    RtBB=np.exp(np.interp(np.log(halfgap),np.log(hgBBlist),np.log(RtBBlist)));
	    # compute model
	    imp_mod,wake_mod=imp_model_resonator(RtBB,frBB,Q);
	
	    imp_mod_list['geometric'].append(imp_mod);
	    # geometric kick factor
	    kickfactorgeom[ihalfgap]=transverse_kick_factor(imp_mod,sigmaz,powerspectrum='gaussian');
	
	# plot transverse kick factor vs halfgap
	fig,ax=init_figure();
	# first geometric one (rough estimate)
	plot(hgapscan*1e3,-kickfactorgeom,'M. Zobov geometric imp.',col[0],r" $ \kappa_\perp $ [V/(C.m)]",ax,2,xlab="half-gap [mm]");
	
	for imat,material in enumerate(materialscan):

            name='coll_'+material;

	    # plot and compare with M. Zobov and E. Metral formulas
	    gam5over4=0.90640247705547705
	    #formElias=-gam5over4*2.**(3/4.)*Z0*length*np.sqrt(c*rhoscan[imat]/(mu0*sigmaz))/(12.*hgapscan**3) # wrong by a factor 2^(1/4)
	    formZobov=gam5over4*c*length*np.sqrt(2*Z0*rhoscan[imat]/sigmaz)/(8.*hgapscan**3); # note: dip + quad here
	    
	    plot(hgapscan*1e3,-kickfactor[imat,:],material+', N. Mounet RW imp.',col[imat+1],r" $ \kappa_\perp $ [V/(C.m)]",ax,2,xlab="half-gap [mm]");
	    #plot(hgapscan*1e3,-formElias,material+', E. Metral formula','x'+col[imat+1],r" $ \kappa_\perp $ [V/(C.m)]",ax,2,xlab="half-gap [mm]");
	    plot(hgapscan*1e3,formZobov,material+', M. Zobov formula','o'+col[imat+1],r" $ \kappa_\perp $ [V/(C.m)]",ax,2,xlab="half-gap [mm]");

		
	ax.grid();
	end_figure(fig,ax,save=root_result+'/plot_kickfactor_'+machine+"coll_hgap_scan",fontsize=25)
	
    
	# DELPHI loops now
	tuneshiftQp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
	tuneshiftm0Qp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
	tuneshiftQpMOS=np.zeros((len(hgapscan),len(Qpscan),nNbMOS,1,(nmaxMOS+1)*(2*lmaxMOS+1)),dtype=complex);
	tuneshiftQpMOSm0=np.zeros((len(hgapscan),len(Qpscan),nNbMOS,1),dtype=complex);
	tuneshiftQpSach=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
	tuneshiftQpSachm0=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1),dtype=complex);

	# MOSES
	RtBBscan=np.exp(np.interp(np.log(hgapscan),np.log(hgBBlist),np.log(RtBBlist)));
	
	if (lxplusbatchDEL=='retrieve'):
	    tuneshiftQpMOS[:,:,:,:,:],tuneshiftQpMOSm0[:,:,:,:],Iscan=MOSES_wrapper(MOSESrootname,
		    RtBBscan,frBB,Q,Qpscan,Nbscan,[omegas],omega0,E,alphap,sigmaz,avbetax,lmaxMOS,nmaxMOS,
		    mmin=-3,mmax=3,taumin=-0.5,taumax=0.5,
		    firstline='MOSES input file # LHC collimator geometric impedance M. Zobov',
		    action=actionMOSES);

	for ihalfgap,halfgap in enumerate(hgapscan):
	
	    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

	    for iscenario,scenario in enumerate(scenarioscan):

        	imp_mod=imp_mod_list[scenario][ihalfgap];

		for iplane,plane in enumerate(['x']):
		    # select Zxdip or Zydip
		    for iw in imp_mod:
			if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);
		    
		    #if (scenario!='geometric'):
		    #	# total impedance from dip+quad
		    #	for iw in imp_mod:
		    #	    if test_impedance_wake_comp(iw,0,0,1-iplane,iplane,plane): Z+=iw.func;
		    

		    #if (halfgap==1e-3)and(scenario=='CFC'): pylab.loglog(freq,Z[:,0],freq,Z[:,1]);pylab.show();
		    
		    for iM,M in enumerate(Mscan):

			# normalization factor for damper
			dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
    			    flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
			dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

			flag_trapz=0; # by default no trapz method

        		if (M==1): nxscan=np.array([0]);flag_trapz=1;
			else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,100),np.arange(M/2-10,M/2+11),
	    		    np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

			tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
			tuneshiftnxSach=np.zeros((len(Qpscan),len(nxscan),len(Nbscan),1,2*lmaxSach+1),dtype=complex);
			ZeffSach=np.zeros((len(Qpscan),len(nxscan),1,2*lmaxSach+1),dtype=complex);

			# DELPHI
			tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
				nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
				a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
				flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
				kmax=kmax,kmaxplot=kmaxplot,crit=1.e-2,abseps=1e-5,flagm0=True,
				lxplusbatch=lxplusbatchDEL,comment=machine+strhgap+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
				queue='1nd',dire=root_result+'/');

			if flagSach:
			    # Sacherer (no damper)
			    tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,:,:,:],tuneshiftnxSach,tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,:,:],ZeffSach=sacherer(imp_mod,
				    Qpscan,nxscan,Nbscan,[omegas],M,omega0,eval('Q'+plane),gamma,eta,taub,lmaxSach,
				    particle='proton',modetype='sinusoidal',compname='Z'+plane+'dip');

	if flagSach:
	    # save Sacherer tuneshifts
	    fileSach=open(root_result+'/Sacherer.txt','w');
	    pick.dump(tuneshiftQpSach,fileSach);
	    pick.dump(tuneshiftQpSachm0,fileSach);
	    fileSach.close();
	else:
	    # load Sacherer tuneshifts
	    fileSach=open(root_result+'/Sacherer.txt','r');
	    tuneshiftQpSach=pick.load(fileSach);
	    tuneshiftQpSachm0=pick.load(fileSach);
	    fileSach.close();
	

	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x']):

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			for Nb in Nbscanplot:

	        	    # plots vs Q'
			    for ihalfgap,halfgap in enumerate(hgapscan):

				strhgapleg='halfgap '+str(1e3*halfgap)+'mm';
				strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

				# initialize plots vs Qp
				figQpm0,axQpm0=init_figure(axes=[0.15,0.1,0.8,0.85]);
				figQp=[];axQp=[];
				for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

				# output file name for plots vs Qp
				fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				strpart=['Re','Im'];
				for ir,r in enumerate(['real','imag']):

				    for iscenario,scenario in enumerate(scenarioscan):

					# output file name for data vs Qp
					fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
					fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

					ts=getattr(tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
					Sachstr='';MOSESstr='';
					if damp==0:
					    # compare with Sacherer most unstable mode
					    tsSach=getattr(tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
					    Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
					    
					    if (M==1)and(scenario=='geometric'):
					    	tsMOS=getattr(tuneshiftQpMOS[ihalfgap,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
					    	MOSESstr=" MOSES_"+strpart[ir]+"_tuneshift"
					    	data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1)),tsMOS.reshape((-1,1))));
					    else:
					    	data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
					
					else:
					    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					
					write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift"+Sachstr+MOSESstr)

					sgn=1;sgnstr='';
					if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
					if damp==0:
					    plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
					    if (M==1)and(scenario=='geometric'):
					    	plot(Qpscan,np.squeeze(tsMOS),'MOSES, '+strhgapleg+', '+scenario,':'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

					
					# real tune shift of mode 0
					if (ir==0):
					    ts=getattr(tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
					    Sachstr='';MOSESstr='';
					    
					    if damp==0:
						# compare with Sacherer most unstable mode
						tsSach=getattr(tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0],r);
						Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift_mode0"
						
						if (M==1)and(scenario=='geometric'):
					    	    tsMOS=getattr(tuneshiftQpMOSm0[ihalfgap,:,pylab.mlab.find(Nbscan==Nb),0],r);
					    	    MOSESstr=" MOSES_"+strpart[ir]+"_tuneshift_mode0"
					    	    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1)),tsMOS.reshape((-1,1))));
						else:
					    	    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
					    
					    else:
						data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					    
					    write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr+MOSESstr)

					    sgn=1;sgnstr='';
					    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					    plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");
					    if damp==0:
					    	plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");
					    	if (M==1)and(scenario=='geometric'):
					    	    plot(Qpscan,np.squeeze(tsMOS),'MOSES, '+strhgapleg+', '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");

					#if (M==1):
				    #	# compare with HEADTAIL
				    #	nsl=50;npr=100000;nlin=1; # HEADTAIL parameters for comparison
				    #	fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip_quad';
				    #	rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/SPS/SPS_1b_ntwake10_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_I"+float_to_str(Nb*1.08333333333/1e11)+"_pre1_drate"+float_to_str(damp)+"_flagdamp0";
				    #	#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
				    #	#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
				    #	sufHEADTAIL="_Sussix_aver_most_tau.txt";
				    #	s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
				    #	fact=1;
                                    #	if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0
                                    #	plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

				    # finish plots vs Qp
				    if (ir==0):
					if flagsave: end_figure(figQpm0,axQpm0,save=fileoutplotQpm0+'_'+r)
					else: end_figure(figQpm0,axQpm0);

				for ir,r in enumerate(['real','imag']):
				    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
				    else: end_figure(figQp[ir],axQp[ir]);

	        	    
			    # plot real tune shifts of mode 0 vs half-gap
			    for Qp in Qpplothg:

				iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];
				
				# initialize plots vs half-gap
				fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);

				# output file name for plots vs hg
				fileoutplothg=root_result+'/plot_vs_hg_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV_Qp'+float_to_str(Qp)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				strpart=['Re','Im'];
				r='real';ir=0;

				for iscenario,scenario in enumerate(scenarioscan):

				    # output file name for data vs halfgap
				    fileoutdatahg=root_result+'/data_vs_hg_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV_Qp'+float_to_str(Qp)+'_'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    # tune shift of mode 0
				    ts=getattr(tuneshiftm0Qp[:,iscenario,iplane,iM,iQp,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
				    Sachstr='';MOSESstr='';
				    if damp==0:
					# compare with Sacherer most unstable mode
					tsSach=getattr(tuneshiftQpSachm0[:,iscenario,iplane,iM,iQp,pylab.mlab.find(Nbscan==Nb),0],r);
					Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift_mode0"
					
					if (M==1)and(scenario=='geometric'):
					    tsMOS=getattr(tuneshiftQpMOSm0[:,iQp,pylab.mlab.find(Nbscan==Nb),0],r);
					    MOSESstr=" MOSES_"+strpart[ir]+"_tuneshift_mode0"
					    data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1)),tsMOS.reshape((-1,1))));
					else:
					    data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
				    
				    else:
				    	data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1))));
				    write_ncol_file(fileoutdatahg+'_'+r+'.dat',data,header="halfgap[m]\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr+MOSESstr)

				    sgn=1;sgnstr='-';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    
				    plot(hgapscan*1e3,np.abs(np.squeeze(sgn*ts)),'DELPHI, '+scenario,col[iscenario],"$ |"+strpart[ir]+"(Q-Q_0)| $ ",ax,2,xlab=" half-gap [mm] ");
				    if damp==0:
				    	plot(hgapscan*1e3,np.abs(np.squeeze(sgn*tsSach)),'Sacherer, '+scenario,'--'+col[iscenario],"$ |"+strpart[ir]+"(Q-Q_0)| $ ",ax,2,xlab=" half-gap [mm] ");
					if (M==1)and(scenario=='geometric'):
				    	    plot(hgapscan*1e3,np.abs(np.squeeze(tsMOS)),'MOSES, '+scenario,'x'+col[iscenario],"$ |"+strpart[ir]+"(Q-Q_0)| $ ",ax,2,xlab=" half-gap [mm] ");

				    #if (M==1):
				#	# compare with HEADTAIL
				#	nsl=50;npr=100000;nlin=1; # HEADTAIL parameters for comparison
				#	fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip_quad';
				#	rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/SPS/SPS_1b_ntwake10_nsl"+str(nsl)+"_npr"+str(npr)+"_nlin"+str(nlin)+"_I"+float_to_str(Nb*1.08333333333/1e11)+"_pre1_drate"+float_to_str(damp)+"_flagdamp0";
				#	#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
				#	#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
				#	sufHEADTAIL="_Sussix_aver_most_tau.txt";
				#	s=read_ncol_file(rootHEADTAIL+sufHEADTAIL,ignored_rows=1);
				#	fact=1;
                                #	if (ir==1): fact=1./omega0; # for imaginary part, divide by omega0
                                #	plot(s[:,0],fact*s[:,3*iplane+ir+1],'HEADTAIL, '+scenario,'x'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

				# finish plots vs half-gap
				if flagsave: end_figure(fig,ax,save=fileoutplothg+'_'+r)
				else: end_figure(fig,ax);

			    # TMCI plots
			    for Qp in Qpplothg:

				iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];
				
				for ihalfgap,halfgap in enumerate(hgapscan):

				    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

				    for iscenario,scenario in enumerate(scenarioscan):

					fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(round(10.*Qp)/10.)+'_converged'+strnorm[flagnorm]+'_'+plane;
					patcol=['.b','b'];
					if scenario=='CFC': ylim=([-0.01,0.01],[-0.001,0.001]);
					elif scenario=='W': ylim=([-0.001,0.001],[-0.0001,0.0001]);
					else: ylim=([-0.001,0.001],[-0.00001,0.00001]);

		    			for ir,r in enumerate(['real','imag']):

					    fig,ax=init_figure();

					    ts=tuneshiftQp[ihalfgap,iscenario,iplane,iM,iQp,idamp,:,0,0,:];

					    plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg='DELPHI',patcol=patcol[ir],xlab='Nb [p+/b]',
						title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.)+", half-gap="+str(halfgap*1e3)+"mm, "+scenario,ms=1,ylim=ylim[ir]);
					    
					    if (damp==0)and((M==1)and(scenario=='geometric')):
						tsMOS=getattr(tuneshiftQpMOS[ihalfgap,iQp,:,0,:],r);
						sgn=1;
						if (ir==1): sgn=-1; # invert sign of imaginary part for MOSES (different convention)
						plot_TMCI(Iscan,sgn*tsMOS/Qs,ax,part=r,leg='MOSES',patcol='xr',xlab='Nb [p+/b]',
							    title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.)+", half-gap="+str(halfgap*1e3)+"mm, "+scenario,ms=1,ylim=ylim[ir]);
	
					    if flagsave: end_figure(fig,ax,save=fileoutplotTMCI+'_'+r,fontsize=25);
					    else: end_figure(fig,ax,fontsize=25);


    if not(flagsave): pylab.show();
