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

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7e12);

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/HLLHC_HOM_crab';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=True; # to plot impedances

    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
        
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=20; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    patcol=['xb','or','+g','vm','^k','.c','dy']; # pattern + color
    linetype=['-','--',':'];
    nevery=20; # downsampling of the impedance (take less points than in the full model)

    # scan definition
    subscan=np.array([0]);strsubscan='_HLLHC_round_Crab';margin_factor=1;
    scenarioscan=np.array(['_HLLHC_round_Crab']);
    legscen=np.array(['HL-LHC round (crab cav., no wire)']);
    squeezescan=np.array(['0p15m_round']);
    Escan=np.array([7e12]);
    optionCrabscan=np.array(['']);
    optionBBCscan=np.array([0]);
    
    param_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat']);
    settings_filename_coll_scan=np.array([path_here+'Coll_settings/collgaps_HLLHC_baseline_from_Roderik_modifNico_material_names.dat']);
    
    # HOM parameters scan definitions
    Qscan=np.array([100.,1e3,1e4]);
    Rscan=np.concatenate(([1e4,1e5,5e5,1e6],np.arange(2e6,1e7+2e6,2e6)));
    Rscanplot=np.array([0,1e4,1e5,1e6,1e7]); # values taken for some plots
    fscan=np.array([1e6,2e6,5e6,1e7,2e7,5e7,1e8,2e8,5e8,1e9,1.5e9,2e9]);
    level=0.01; # level of tolerance for added growth rate or subtracted TMCI threshold, due to HOM.
    # Example: 0.01 means we tolerate 1% more growth rate or 1% less TMCI threshold.
    Rscan_withzero=np.concatenate(([0],Rscan)); # point without HOM added (0 shunt impedance))
    Rscanfine=np.arange(0,Rscan[-1]+1e2,1e2); # finer mesh for interpolation
    
    # chromaticity scan definitions
    Qpscan=np.array([0,14,15,16]);
    Qpaver=np.array([14,15,16]);
    Qpscanplot=np.array([0,15]);
    iQpaver=select_in_table(Qpaver,Qpscan);print iQpaver
    QpscanTMCI=np.array([0]);Qpscancoupledbunch=np.array([0]);

    dampscan=np.array([0,0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscanplot=np.array([1.7e11]); # intensity scan for plot vs Qp
    #Nbscan=np.array([]);dampscan=np.array([]);
    Mscan=np.array([1,3564]); # scan on number of bunches
    #Mscan=np.array([1]); # scan on number of bunches

    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    tuneshiftQp=np.zeros((len(scenarioscan),len(Qscan),len(fscan),len(Rscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp=np.zeros((len(scenarioscan),len(Qscan),len(fscan),len(Rscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    tuneshiftQp_noHOM=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    tuneshiftm0Qp_noHOM=np.zeros((len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
    
    for iscenario,scenario in enumerate(scenarioscan[subscan]):
    	
	# parameters
	machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=Escan[subscan[iscenario]]);
    	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
	beta=np.sqrt(1.-1./gamma**2);
	print "scenario: ",scenario
	
	# compute imp. model
	collscen='_HLLHC';
	imp_mod,wake_mod=HLLHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll_scan[subscan[iscenario]],
		settings_filename_coll_scan[subscan[iscenario]],dire=path_here+"LHC_elements/",
		commentcoll=collscen,direcoll='Coll'+collscen+'_v2/',
		lxplusbatch=lxplusbatchImp,beam=beam,squeeze=squeezescan[subscan[iscenario]],
		wake_calc=wake_calc,optionCrab=optionCrabscan[subscan[iscenario]],margin_factor=margin_factor,
		optionBBC=optionBBCscan[subscan[iscenario]],flagplot=flagplot,root_result=root_result,commentsave=scenario)

	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # dump into a file
	    filemodel=open(root_result+'/impedances'+scenario+'.txt','w');
	    pick.dump(imp_mod,filemodel);
	    filemodel.close();
	    
	    if (wake_calc):
	    	# wake for HEADTAIL
		write_wake_HEADTAIL(wake_mod,root_result+"/wakeforhdtl_PyZbase_Allthemachine_"+Estr+"_B"+str(beam)+scenario+'.dat',beta=np.sqrt(1.-1./(gamma**2)),ncomp=6)		
    

	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

            # DELPHI on model without HOM
	    tuneshiftQp_noHOM[iscenario,:,:,:,:,:,:,:,:],tuneshiftm0Qp_noHOM[iscenario,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod,
	    	Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],	omega0,Qx,Qy,gamma,eta,a,b,
		taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		abseps=1.e-3,flagm0=True,lxplusbatch=lxplusbatchDEL,
		comment=machine+scenario+'_'+float_to_str(round(E/1e9))+'GeV'+strnorm[flagnorm],
		queue='1nd',dire=root_result+'/',flagQpscan_outside=False);
	    
	    # scan in HOM parameters
	    for iQ,Q in enumerate(Qscan):
	    
	    	for ifr,fr in enumerate(fscan):
	    
		    for iRt,Rt in enumerate(Rscan):
		    
		    	imp_mod_tot=[];wake_mod_tot=[]; # total model (with HOM in)
			add_impedance_wake(imp_mod_tot,imp_mod,1.,1.);
			add_impedance_wake(wake_mod_tot,wake_mod,1.,1.);
			
			if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('launch')):
			
			    # Crab cavity parameters
			    namesref=['Crab'];length_ind=1.625; # individual length for 16 Crab cavities
			    beta_filename=path_here+"LHC_elements/HLLHC_beta_length_B"+str(beam)+"_sq"+squeezescan[subscan[iscenario]]+".dat";
			    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
			    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
			    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
			    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
			    # reorder such that the names match with namesref
			    ind=find_ind_names(namesref,names);
			    # calculate HOM
			    imp_mod_HOM,wake_mod_HOM=imp_model_resonator(Rt,fr,Q,beta=beta,wake_calc=wake_calc,
				    fpar=freq_param(ftypescan=1,fmin=fr*(1-2./np.sqrt(Q)),fmax=fr*(1+2./np.sqrt(Q)),fsamplin=fr/(Q*400),fadded=np.concatenate((10**np.arange(2,10.01,0.01),[10,1e11,1e12]))),
				    zpar=z_param(),listcomp=['Zlong','Zxdip','Zydip']);
			    # multiply HOM imp & wake by number of such elements and add (with beta function weight)
			    # to the total model
			    print "    nb of cavities=",length[ind[0]]/length_ind,", beta function ratio:",betax[ind[0]]/avbetax,betay[ind[0]]/avbetay;
			    multiply_impedance_wake(imp_mod_HOM,length[ind[0]]/length_ind);
			    multiply_impedance_wake(wake_mod_HOM,length[ind[0]]/length_ind);
			    add_impedance_wake(imp_mod_tot,imp_mod_HOM,betax[ind[0]]/avbetax,betay[ind[0]]/avbetay);
			    add_impedance_wake(wake_mod_tot,wake_mod_HOM,betax[ind[0]]/avbetax,betay[ind[0]]/avbetay);

			    # plot and compare with and without added HOM
			    maxratio_sb=plot_compare_imp_model([imp_mod,imp_mod_tot],['total model','total model + HOM, Q='+float_to_str(Q)+", fr="+float_to_str(fr)+", R="+float_to_str(Rt)],
				    listcomp=['Zlong','Zxdip','Zydip'],
				    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan+"_Q"+float_to_str(Q)+"_fr"+float_to_str(fr)+"_Rt"+float_to_str(Rt),
				    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan+"_Q"+float_to_str(Q)+"_fr"+float_to_str(fr)+"_Rt"+float_to_str(Rt),
				    ylim=[1e5,1e12],yliml=[1,1e9],bounds=[40e6,2e9]);

			#maxratio_mb=plot_compare_imp_model(imp_mod_list,legscen[subscan],listcomp=['Zlong','Zxdip','Zydip'],
			#    saveimp=root_result+'/plot_imp_'+machine+"_scenarios"+strsubscan,
			#    saveratio=root_result+'/plot_imp_ratio_'+machine+"_scenarios"+strsubscan,
			#    ylim=[1e5,1e10],bounds=[8e3,40e6],legpos=9);

			# plot and compare all scenarios with 2012 wake
			#if wake_calc:
			#    maxratio_w=plot_compare_imp_model(wake_mod_list,legscen[subscan],listcomp=['Wlong','Wxdip','Wydip'],
			#	saveimp=root_result+'/plot_wake_'+machine+"_scenarios"+strsubscan,
			#	saveratio=root_result+'/plot_wake_ratio_'+machine+"_scenarios"+strsubscan,
			#	xlim=[1e-5,1e6],ylim=[1e8,1e19],yliml=[1e6,1e19],legpos=0);

        		# DELPHI on total model with HOM
			tuneshiftQp[iscenario,iQ,ifr,iRt,:,:,:,:,:,:,:,:],tuneshiftm0Qp[iscenario,iQ,ifr,iRt,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod_tot,
	    		    Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],omega0,Qx,Qy,gamma,
			    eta,a,b,taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
			    flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
			    abseps=1.e-3,flagm0=True,lxplusbatch=lxplusbatchDEL,
			    comment=machine+scenario+'_HOM_Q'+float_to_str(Q)+"_fr"+float_to_str(fr)+"_Rt"+float_to_str(Rt)+"_"+float_to_str(round(E/1e9))+'GeV'+strnorm[flagnorm],
			    queue='2nd',dire=root_result+'/',flagQpscan_outside=False);
		

	    # now the plots (inside loop on scenarios)
	    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

		for iplane,plane in enumerate(['x','y']):

		    for iM,M in enumerate(Mscan):

			for idamp,damp in enumerate(dampscan):

			    if (M==1)and(damp==0):
			    
			    	# TMCI analysis
				for Qp in QpscanTMCI:
				    
				    iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];print "iQp=",iQp
				    
				    # initialize plot of threshold vs fr
				    figvsfr,axvsfr=init_figure();
				    # output file name for Rt threshold plot vs fr
				    fileoutplotfr=root_result+'/plot_TMCI_Rt_vs_fr_thres_level'+float_to_str(level)+machine+'_'+Estr+scenario+"_"+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
				
				    for iQ,Q in enumerate(Qscan):

	    				Rthres=np.zeros(len(fscan));

					for ifr,fr in enumerate(fscan):

					    Nbthres=np.zeros(len(Rscan)+1);

					    # initialize plot of threshold vs Rt
					    figvsRt,axvsRt=init_figure();
					    # output file name for threshold plot vs Rt
					    fileoutplotRt=root_result+'/plot_TMCI_Nb_vs_Rt_thres_'+machine+'_'+Estr+scenario+'_Q'+float_to_str(Q)+"_fr"+float_to_str(fr)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
					
					    for iRt,Rt in enumerate(Rscan_withzero):
			    	
						if Rt in Rscanplot:
						    # TMCI plot for a few values of Rt
						    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+Estr+scenario+'_Q'+float_to_str(Q)+"_fr"+float_to_str(fr)+"_Rt"+float_to_str(Rt)+"_"+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_converged'+strnorm[flagnorm]+'_'+plane;
						    patcol2=['.b','b'];
						    ylim=([-5,3],[-1,0.01]);

		    				    for ir,r in enumerate(['real','imag']):

							fig,ax=init_figure();

							if iRt==0: ts=tuneshiftQp_noHOM[iscenario,iplane,iM,iQp,idamp,:,0,0,:];
							else: ts=tuneshiftQp[iscenario,iQ,ifr,iRt-1,iplane,iM,iQp,idamp,:,0,0,:];

							plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg='DELPHI',patcol=patcol2[ir],xlab='Nb [p+/b]',
							    title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.)+', Q='+float_to_str(Q)+", fr="+float_to_str(fr)+", R="+float_to_str(Rt)+", "+scenario,ms=1,ylim=ylim[ir]);

							end_figure(fig,ax,save=flagsave*(fileoutplotTMCI+'_'+r),fontsize=25);

						# find intensity threshold
						if iRt==0: Nbthres[0]=find_intensity_threshold(Nbscan,tuneshiftQp_noHOM[iscenario,iplane,iM,iQp,idamp,:,0,0,0]*omega0,thresgrowth=0.01);
						else: Nbthres[iRt]=find_intensity_threshold(Nbscan,tuneshiftQp[iscenario,iQ,ifr,iRt-1,iplane,iM,iQp,idamp,:,0,0,0]*omega0,thresgrowth=0.01);

					    # plot threshold vs Rt
					    plot(Rscan_withzero/1e3,Nbthres,'DELPHI','-b',"TMCI threshold (nb p/b)",axvsRt,2,xlab="HOM transverse shunt impedance [k $ \Omega/ $ m]");

					    # finish plot vs Rt
					    end_figure(figvsRt,axvsRt,save=flagsave*(fileoutplotRt))
					    
					    # reinterpolate thresholds on a finer mesh
					    Nbthresfine=np.interp(Rscanfine,Rscan_withzero,Nbthres);
					    
					    # find where TMCI threshold becomes less than (1-level)*(threshold without HOM)
					    ratio=Nbthresfine/Nbthres[0];
					    ind=pylab.mlab.find(ratio<=(1-level));
					    if len(ind)>0: Rthres[ifr]=Rscanfine[ind[0]];
					    else:
					    	Rthres[ifr]=Rscan[-1];print "Warning (TMCI): Rt not high enough for Q=",Q," and fr=",fr," !"
					    
					# plot Rt threshold vs fr
					plot(fscan/1e6,Rthres/1e3,'Q='+float_to_str(Q),'-'+patcol[iQ],"Max. shunt impedance [k $\Omega/ $ m]",axvsfr,2,xlab="HOM resonance frequency [MHz]");

				    # finish plot Rt threshold vs fr
				    end_figure(figvsfr,axvsfr,save=flagsave*(fileoutplotfr))


			    if (M==3564)and(damp==0):
			    
			    	# Coupled-bunch modes analysis
				for Nb in Nbscanplot:
				    
				    iNb=pylab.mlab.find(Nbscan==Nb);iNb=iNb[0];print "iNb=",iNb
				
				    for Qp in Qpscancoupledbunch:

					iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];print iQp

					# initialize plot of threshold vs fr
					figvsfr,axvsfr=init_figure();
					# output file name for Rt threshold plot vs fr
					fileoutplotfr=root_result+'/plot_coupledbunch_Rt_vs_fr_thres_level'+float_to_str(level)+machine+'_'+Estr+scenario+"_"+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

					for iQ,Q in enumerate(Qscan):

	    				    Rthres=np.zeros(len(fscan));

					    for ifr,fr in enumerate(fscan):

						# initialize plot of threshold vs Rt
						figvsRt,axvsRt=init_figure();
						# output file name for threshold plot vs Rt
						fileoutplotRt=root_result+'/plot_coupledbunch_growthrate_vs_Rt_thres_'+machine+'_'+Estr+scenario+'_Q'+float_to_str(Q)+"_fr"+float_to_str(fr)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

						# plot maximum growth rate vs Rt
						gr=-omega0*np.imag(np.hstack((tuneshiftQp_noHOM[iscenario,iplane,iM,iQp,idamp,iNb,0,0,0],tuneshiftQp[iscenario,iQ,ifr,:,iplane,iM,iQp,idamp,iNb,0,0,0])));
						plot(Rscan_withzero/1e3,gr,'DELPHI','-b',"Growth rate [1/s]",axvsRt,2,xlab="HOM transverse shunt impedance [k $ \Omega/ $ m]");

						# finish plot vs Rt
						end_figure(figvsRt,axvsRt,save=flagsave*(fileoutplotRt))

						# reinterpolate thresholds on a finer mesh
						grfine=np.interp(Rscanfine,Rscan_withzero,gr);

						# find where growth rate becomes more than (1+level)*(growth rate without HOM)
						ratio=grfine/gr[0];
						ind=pylab.mlab.find(ratio>=(1+level));
						if len(ind)>0: Rthres[ifr]=Rscanfine[ind[0]];
						else:
					    	    Rthres[ifr]=Rscan[-1];print "Warning (coupled-bunch): Rt not high enough for Q=",Q," and fr=",fr," !"

					    # plot Rt threshold vs fr
					    plot(fscan/1e6,Rthres/1e3,'Q='+float_to_str(Q),'-'+patcol[iQ],"Max. shunt impedance [k $ \Omega/ $ m]",axvsfr,2,xlab="HOM resonance frequency [MHz]");

					# finish plot Rt threshold vs fr
					end_figure(figvsfr,axvsfr,save=flagsave*(fileoutplotfr))


			    if (damp>0):
			    
			    	# high chroma + damper modes analysis
				for Nb in Nbscanplot:
				    
				    iNb=pylab.mlab.find(Nbscan==Nb);iNb=iNb[0];print iNb
				    
				    # initialize plot of threshold vs fr
				    figvsfr,axvsfr=init_figure();
				    # output file name for Rt threshold plot vs fr
				    fileoutplotfr=root_result+'/plot_highchroma_Rt_vs_fr_thres_level'+float_to_str(level)+machine+'_'+Estr+scenario+"_Qpmin"+float_to_str(Qpaver[0])+"_Qpmax"+float_to_str(Qpaver[-1])+"_"+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				
				    for iQ,Q in enumerate(Qscan):

	    				Rthres=np.zeros(len(fscan));

					for ifr,fr in enumerate(fscan):

					    # initialize plot of threshold vs Rt
					    figvsRt,axvsRt=init_figure();
					    # output file name for threshold plot vs Rt
					    fileoutplotRt=root_result+'/plot_highchroma_growthrate_vs_Rt_thres_'+machine+'_'+Estr+scenario+"_Qpmin"+float_to_str(Qpaver[0])+"_Qpmax"+float_to_str(Qpaver[-1])+"_"+'_Q'+float_to_str(Q)+"_fr"+float_to_str(fr)+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
					
					    # find maximum growth rate for each shunt impedance (scan also Qpaver and look for maximum)
					    gr=np.zeros(len(Rscan_withzero));
					    for iRt,Rt in enumerate(Rscan_withzero):
					    
					    	if iRt==0: gr[0]=-omega0*np.min(np.imag(tuneshiftQp_noHOM[iscenario,iplane,iM,iQpaver,idamp,iNb,0,0,0]));
						else: gr[iRt]=-omega0*np.min(np.imag(tuneshiftQp[iscenario,iQ,ifr,iRt-1,iplane,iM,iQpaver,idamp,iNb,0,0,0]));
			    	    	    
					    # plot maximum growth rate vs Rt
					    plot(Rscan_withzero/1e3,gr,'DELPHI','-b',"Growth rate [1/s]",axvsRt,2,xlab="HOM transverse shunt impedance [k $ \Omega/ $ m]");

					    # finish plot vs Rt
					    end_figure(figvsRt,axvsRt,save=flagsave*(fileoutplotRt))
					    
					    # reinterpolate thresholds on a finer mesh
					    grfine=np.interp(Rscanfine,Rscan_withzero,gr);
					    
					    # find where growth rate becomes more than (1+level)*(growth rate without HOM)
					    ratio=grfine/gr[0];
					    ind=pylab.mlab.find(ratio>=(1+level));
					    if len(ind)>0: Rthres[ifr]=Rscanfine[ind[0]];
					    else:
					    	Rthres[ifr]=Rscan[-1];print "Warning (high chroma): Rt not high enough for Q=",Q," and fr=",fr," !"
					    
					# plot Rt threshold vs fr
					plot(fscan/1e6,Rthres/1e3,'Q='+float_to_str(Q),'-'+patcol[iQ],"Max. shunt impedance [k $ \Omega/ $ m]",axvsfr,2,xlab="HOM resonance frequency [MHz]");

				    # finish plot Rt threshold vs fr
				    end_figure(figvsfr,axvsfr,save=flagsave*(fileoutplotfr))


    if not(flagsave): pylab.show();
