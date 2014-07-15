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
import numpy as np
import pickle as pick
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *
from SPS_param import SPS_param
sys.path.append("../LHC_impedance_and_scripts/")
from LHC_coll_imp import *


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=450e9,optics='Q20');
    beta=np.sqrt(1.-1./(gamma**2))

    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'_coll';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    flagplotQp=1; # 1 to do the plots vs Qp
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    #flagSach=(lxplusbatchDEL=='launch'); # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    flagSach=0; # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    
    wake_calc=False
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=50; # number of plotted eigenvalues (for TMCI)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];
    lmaxSach=1; # max headtail mode for Sacherer

    # scan definition
    Escan=[26e9,450e9];
    scenarioscan=['_Q20','_Q20_coll','_Q26','_Q26_coll'];
    #modelscan=['_SPS_Q20_2012_full_modifNico','_SPS_Q20_2012_full_modifNico',
    #	'_kickers_wall_BPMs_flanges_cavities_Q26_modifNico','_kickers_wall_BPMs_flanges_cavities_Q26_modifNico'];
    modelscan=['_SPS_Q20_2012_full_refined_below1MHz_modifNico','_SPS_Q20_2012_full_refined_below1MHz_modifNico',
    	'_SPS_Q26_2012_full_refined_below1MHz_modifNico','_SPS_Q26_2012_full_refined_below1MHz_modifNico'];
    collscan=[False,True,False,True];
    opticscan=['Q20','Q20','Q26','Q26']
    legscen=['total Q20','total Q20 with collimators','total Q26','total Q26 with collimators'];
    csiscan=np.arange(-0.5,0.6,0.1);
    #csiscan=np.array([0]);
    csiplotscan=np.array([0]);

    dampscan=np.array([0,0.05,0.1]); # damper gain scan
    Nbscan=np.arange(1.e10,5.1e11,5.e9); # intensity scan
    Nbscanplot=np.array([1.2e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    tuneshiftQp=np.zeros((len(scenarioscan),2,len(Mscan),len(csiscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftQpSach=np.zeros((len(scenarioscan),2,len(Mscan),len(csiscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
    ZeffSach=np.zeros((len(scenarioscan),2,len(Mscan),len(csiscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
    tuneshiftQpSachm0=np.zeros((len(scenarioscan),2,len(Mscan),len(csiscan),len(Nbscan),1),dtype=complex);
    
    for iE,E in enumerate(Escan):
    
	imp_mod_list=[]; # complete list of impedance scenarios
	wake_mod_list=[];# complete list of wake scenarios

	for iscenario,scenario in enumerate(scenarioscan):

	    imp_mod_list_details=[];

	    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=E,optics=opticscan[iscenario]);
	    beta=np.sqrt(1.-1./(gamma**2))
    	    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	    print scenario,opticscan[iscenario],avbetax,avbetay

	    # read model from Carlo
	    imp_mod=imp_model_from_file(path_here+'imp_model_SPS_from_CZannini/Zxdip'+modelscan[iscenario]+'.txt','Zxdip');
	    imp_mod1=imp_model_from_file(path_here+'imp_model_SPS_from_CZannini/Zydip'+modelscan[iscenario]+'.txt','Zydip');imp_mod.append(imp_mod1[0]);
	    imp_mod1=imp_model_from_file(path_here+'imp_model_SPS_from_CZannini/Zxquad'+modelscan[iscenario]+'.txt','Zxquad');imp_mod.append(imp_mod1[0]);
	    imp_mod1=imp_model_from_file(path_here+'imp_model_SPS_from_CZannini/Zyquad'+modelscan[iscenario]+'.txt','Zyquad');imp_mod.append(imp_mod1[0]);

	    if collscan[iscenario]:
		# compute model for collimators
		param_filename_coll=path_here+'Coll_settings/collgaps_'+opticscan[iscenario]+'_modifNico_materials.dat';
		beta_filename_coll=param_filename_coll;
		settings_filename_coll=param_filename_coll;

		imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,
		    avbetax,avbetay,param_filename_coll,settings_filename_coll,
		    beta_filename_coll,wake_calc=wake_calc,ftypescan=2,nflog=100,zpar=z_param(),namesref=None,
		    BPM=False,lxplusbatch=lxplusbatchImp,comment='SPS_'+opticscan[iscenario]+'_'+Estr,dire='Coll_SPS_'+opticscan[iscenario]+'_'+Estr+'/');

		# total collimators impedance
		imp_mod_coll=[];
		add_impedance_wake(imp_mod_coll,imp_mod_coll_RW,1,1);
		add_impedance_wake(imp_mod_coll,imp_mod_coll_geom,1,1);
		# add this to total impedance
		add_impedance_wake(imp_mod,imp_mod_coll,1,1);

		imp_mod_list_details.append(imp_mod);
		imp_mod_list_details.append(imp_mod_coll_RW);
		imp_mod_list_details.append(imp_mod_coll_geom);

		legdetails=['total','coll. (RW)','coll. (geom.)'];

		# output coll. longitudinal impedance
		write_imp_wake_mod(imp_mod_coll,'_collimators_RW_geom_'+Estr+'_'+opticscan[iscenario],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],dire=root_result+'/')
		write_imp_wake_mod(imp_mod_coll_RW,'_collimators_RW_'+Estr+'_'+opticscan[iscenario],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],dire=root_result+'/')
		write_imp_wake_mod(imp_mod_coll_geom,'_collimators_geom_'+Estr+'_'+opticscan[iscenario],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],dire=root_result+'/')

	    else:
		imp_mod_list_details.append(imp_mod);
		legdetails=[''];


	    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):
		# plot details and percents
		plot_compare_imp_model(imp_mod_list_details,legdetails,listcomp=['Zxdip','Zydip'],
			saveimp=root_result+'/plot_imp_'+machine+'_'+Estr+scenario+"_details",
			saveratio=root_result+'/plot_imp_ratio_'+machine+'_'+Estr+scenario+"_details",
			xlim=[1e3,2e9],plotpercent=True);


	    # output impedance
	    write_imp_wake_mod(imp_mod,scenario+'_'+Estr,listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],dire=root_result+'/')

    	    imp_mod_list.append(imp_mod);


	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):
    	    # compare all scenarios
	    plot_compare_imp_model(imp_mod_list,legscen,listcomp=['Zxdip','Zydip'],
			saveimp=root_result+'/plot_imp_'+machine+'_'+Estr+"_scenarios",
			saveratio=root_result+'/plot_imp_ratio_'+machine+'_'+Estr+"_scenarios",
			xlim=[1e3,2e9]);

	    #filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    #pick.dump(imp_mod,filemodel);
	    #filemodel.close();

	    # plot and compare all scenarios with 2012 impedance
	    #maxratio=plot_compare_imp_model(imp_mod_list,legscen,listcomp=['Zxdip','Zydip'],
	    #    saveimp=root_result+'/plot_imp_'+machine+"postLS1_scenarios",
	    #    saveratio=root_result+'/plot_imp_ratio_'+machine+"postLS1_scenarios");

	    # plot and compare all scenarios with 2012 wake
	    #maxratio_w=plot_compare_imp_model(wake_mod_list,legscen,listcomp=['Wxdip','Wydip'],
	    #    savewake=root_result+'/plot_wake_'+machine+"postLS1_scenarios",
	    #    saveratio=root_result+'/plot_wake_ratio_'+machine+"postLS1_scenarios",
	    #    xlim=[1e-5,1e6],ylim=[1e8,1e19]);


	    #if False:
	    for iscenario,scenario in enumerate(scenarioscan):

    		# DELPHI loops now
		machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=E,optics=opticscan[iscenario]);
		beta=np.sqrt(1.-1./(gamma**2))
    		avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
		print scenario,opticscan[iscenario],avbetax,avbetay

		for iplane,plane in enumerate(['x','y']):
		    # select Zxdip or Zydip
		    for iw in imp_mod_list[iscenario]:
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

			tuneshiftQp[iscenario,iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
				nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
				a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
				flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
				kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
				lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
				queue='2nd',dire=root_result+'/');

			for iQp,Qp in enumerate(Qpscan):

			    #tuneshiftQp[iscenario,iplane,iM,iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			#	    nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
			#	    a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
			#	    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
			#	    kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,abseps=1e-3,
			#	    lxplusbatch=lxplusbatchDEL,comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV_'+str(M)+'b_csi'+float_to_str(round(100.*csiscan[iQp])/100.)+strnorm[flagnorm]+'_'+plane,
			#	    queue='1nd',dire=root_result+'/');

			    if flagSach:
				# Sacherer (no damper)
				tuneshiftQpSach[iscenario,iplane,iM,:,:,:,:],tuneshiftnxSach,tuneshiftQpSachm0[iscenario,iplane,iM,:,:,:],ZeffSach[iscenario,iplane,iM,:,:,:,:]=sacherer(imp_mod_list[iscenario],
					Qpscan,nxscan,Nbscan,[omegas],M,omega0,eval('Q'+plane),gamma,eta,taub,lmaxSach,
					particle='proton',modetype='sinusoidal',compname='Z'+plane+'dip');

	    if flagSach:
		# save Sacherer tuneshifts
		fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','w');
		pick.dump(tuneshiftQpSach,fileSach);
		pick.dump(tuneshiftQpSachm0,fileSach);
		pick.dump(ZeffSach,fileSach);
		fileSach.close();
	    else:
		# load Sacherer tuneshifts
		fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','r');
		tuneshiftQpSach=pick.load(fileSach);
		tuneshiftQpSachm0=pick.load(fileSach);
		ZeffSach=pick.load(fileSach);
		fileSach.close();

	# now the plots (outside loop on scenarios)
	#if (False)and((lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve'))):
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
					Sachstr='';
					if damp==0:
					    # compare with Sacherer most unstable mode
					    tsSach=getattr(tuneshiftQpSach[iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
					    ZSach=getattr(ZeffSach[iscenario,iplane,iM,:,0,0,lmaxSach],r);
					    Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift Sacherer_"+strpart[ir]+"_Zeff"
					    data=np.hstack((csiscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1)),ZSach.reshape((-1,1))));
					else:
				    	    data=np.hstack((csiscan.reshape((-1,1)),ts.reshape((-1,1))));

					write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="csi\t"+strpart[ir]+"_tuneshift"+Sachstr)
					#else:
					#    s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
					#    Qpscan=s[:,0];ts=s[:,1];

					sgn=1;sgnstr='';
					if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					plot(csiscan,np.squeeze(sgn*ts),'DELPHI, '+legscen[iscenario],col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
					if damp==0:
					    plot(csiscan,np.squeeze(sgn*tsSach),'Sacherer, '+legscen[iscenario],'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");

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

				# finish plot vs Qp
				for ir,r in enumerate(['real','imag']):
				    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
				    else: end_figure(figQp[ir],axQp[ir]);



			# TMCI plots
			for csi in csiplotscan:

			    icsi=pylab.mlab.find(abs(csiscan-csi)<=1e-8);icsi=icsi[0];print icsi;
			    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_d'+float_to_str(damp)+'_csi'+float_to_str(round(1000.*csi)/1000.)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    patcol=['.',''];
			    ylim=([-5,3],[-1.5,0.01]);

		    	    for ir,r in enumerate(['real','imag']):

				fig,ax=init_figure();

				for iscenario,scenario in enumerate(scenarioscan):

				    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=E,optics=opticscan[iscenario]);
				    ts=tuneshiftQp[iscenario,iplane,iM,icsi,idamp,:,0,0,:];

				    plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg=legscen[iscenario],patcol=patcol[ir]+col[iscenario],xlab='Nb [p+/b]',
					title=machine+r", $ \xi = $ "+str(round(1000*csi)/1000.)+', '+str(1./damp)+' turns damping time',ms=1,ylim=ylim[ir]);

				if flagsave: end_figure(fig,ax,save=fileoutplotTMCI+'_'+r,fontsize=25);
				else: end_figure(fig,ax,fontsize=25);


		    # imag. part of TMCI plot for several damping rates, for csi=0, all on same plot
		    for csi in csiplotscan:

			iscenario=0;icsi=pylab.mlab.find(abs(csiscan-csi)<=1e-8);icsi=icsi[0];
			machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=SPS_param(E0,E=E,optics=opticscan[iscenario]);

			for ir,r in enumerate(['real','imag']):

			    fig,ax=init_figure();

			    fileoutplotTMCI2=root_result+'/plot_TMCI_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+'_'+str(M)+'b_dscan_csi'+float_to_str(round(1000.*csi)/1000.)+'_converged'+strnorm[flagnorm]+'_'+plane;
			    pat=['.',''];
			    ylim=([-5,3],[-1,0.01]);

			    for idamp,damp in enumerate(dampscan):

				ts=tuneshiftQp[iscenario,iplane,iM,icsi,idamp,:,0,0,:];

				plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg='',patcol=pat[ir],xlab='Nb [p+/b]',
				    title=machine+r", $ \xi = $ "+str(round(1000*csi)/1000.),ms=1,ylim=ylim[ir]);

			    if flagsave: end_figure(fig,ax,save=fileoutplotTMCI2+'_'+r,fontsize=25);
			    else: end_figure(fig,ax,fontsize=25);


    if not(flagsave): pylab.show();
