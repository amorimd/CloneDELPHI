#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);
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
from LHC_coll_imp import *
from scipy import interpolate


if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=6.5e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    beam='1';
    # directory (inside DELPHI_results/[machine]) where to put the results
    root_result=path_here+'../../../DELPHI_results/'+machine+'/LHC_optim_IR3_IR7';
    os.system("mkdir -p "+root_result);
    
    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=0; # 0 to avoid computing (simply plot from existing data)
    flagplot=False;
    nevery=2;
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];

    # scan definition
    scenario='_mm_kept';
    param_filename_coll=path_here+'Coll_settings/collgaps_fromRoderik_modifNico_materialnames_feb2014'+scenario+'.dat';
    settings_filename_coll_prefix=path_here+'Coll_settings/collgaps_settings_sigma_feb2014'+scenario+'_offsetIR3_';
    beta_filename_coll=param_filename_coll;
    squeeze='0p55m_10m_0p55m_10m';
    BPMflag=True;

    # IR3 offset sigma scan
    sigmaIR3scan=[0,2,4];sigmaIR3scan_plot=[0,4];
    # IR7 retraction sigma scan
    sigmaIR7scan=[2,2.5,3,3.5,4,4.5,5,5.5,6];sigmaIR7scan_plot=[2,4,6];
    sigmaIR7scan_str=[];

    # relevant 2012 instability data - based on files ../Mesures_LHC/instability_table_B2V_pos_oct_flattop.csv &
    # ../Mesures_LHC/instability_table_B2H_pos_oct_flattop.csv (Q' > 10)
    # and file ../Mesures_LHC/instability_table_B2H_neg_oct_squeeze.csv (Q'>5)
    # There are also a 3 flat top instabilities taken from an LMC talk
    # by G. Arduini in August 2012 (slides in ../Docs/lmc_145c_talk_instabilities2012_LHC_Gianluigi.pdf,
    # slide 14&15, fills 2919, 2920 & 2932 - damper gain from Timber & trim editor)
    en2012=4;M_2012=1782;
    # negative octupole polarity
    dataoct_neg2012=np.array([58.,200.,200.]);
    meanoct_neg2012=np.average(dataoct_neg2012);erroct_neg2012=np.sqrt(np.var(dataoct_neg2012));
    dataemit_neg2012=np.array([2.3,2.5,2.5]);
    meanemit_neg2012=np.average(dataemit_neg2012);erremit_neg2012=0.5;
    dataNb_neg2012=np.array([1.4,1.47,1.46])*1e11;
    dataQp_neg2012=np.array([5.9,7.,7.]);
    dataplane_neg2012=np.array(['x','x','x']);
    datadamp_neg2012=1./np.array([50.,100.,100.]);
    # positive octupole polarity
    dataoct_pos2012=np.array([209.,487.,510.,35.,35.,510.]);
    meanoct_pos2012=np.average(dataoct_pos2012);erroct_pos2012=np.sqrt(np.var(dataoct_pos2012));
    dataemit_pos2012=np.array([2.3,2.3,2.3,2.2,2.2,2.5]);
    meanemit_pos2012=np.average(dataemit_pos2012);erremit_pos2012=0.5;
    dataNb_pos2012=np.array([1.64,1.64,1.63,1.64,1.64,1.44])*1e11;
    dataQp_pos2012=np.array([15.4,18.3,10.3,13.6,17.8,9.]);
    dataplane_pos2012=np.array(['x','x','x','x','y','x']);
    datadamp_pos2012=1./np.array([100.,100.,100.,100.,100.,100.]);
    
    # scan parameters for DELPHI computation at 4 TeV for those experimental cases
    oct_2012=np.hstack((dataoct_neg2012,dataoct_pos2012));
    emit_2012=np.hstack((dataemit_neg2012,dataemit_pos2012));
    Nb_2012=np.hstack((dataNb_neg2012,dataNb_pos2012));
    Qp_2012=np.hstack((dataQp_neg2012,dataQp_pos2012));
    damp_2012=np.hstack((datadamp_neg2012,datadamp_pos2012));
    plane_2012=np.hstack((dataplane_neg2012,dataplane_pos2012));
    ind_neg2012=np.arange(3);ind_pos2012=np.arange(3,9); # indices in the above tables, for resp. neg. and oct. polarities

    # scan parameters for DELPHI computation at 6.5TeV
    enLS1=6.5;octLS1=570.;
    dampscan=np.array([0.02]); # damper gain scan
    Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Mscan=np.array([1,1782,3564]); # scan on number of bunches
    Qpscan=np.array([14,15,16]);
    Qpaver=np.array([14,15,16]);
    iQpaver=select_in_table(Qpaver,Qpscan); print iQpaver
    
    wake_calc=False; # True -> compute wake as well (otherwise only imp.)
    
    # compute first with 2012 imp. model
    param_filename_coll_2012=path_here+'Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
    settings_filename_coll_2012=path_here+'Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';
    # fixed parameters
    machine,E_2012,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=en2012*1e12);
    avbetax_2012=R/Qx;avbetay_2012=R/Qy; # average beta functions used
    # model
    imp_mod_2012,wake_mod_2012=LHC_imp_model_v2(E_2012,avbetax_2012,avbetay_2012,param_filename_coll_2012,
	    settings_filename_coll_2012,TDIcoating='preLS1',dire=path_here+"LHC_elements/",commentcoll='_2012_v2',direcoll='Coll_2012_v2/',
	    lxplusbatch=lxplusbatchImp,BPM=False,beam=beam,squeeze='0p6m_3m_0p6m_3m',
	    wake_calc=wake_calc,flagplot=flagplot,root_result=root_result,commentsave='_2012_v2')
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    tuneshiftQp_2012=np.zeros((len(Qp_2012),1,1,1,1,kmaxplot),dtype=complex);
    factor_2012=np.zeros(len(Qp_2012));
    for iQp,Qp in enumerate(Qp_2012):

	# select parameters
	plane=plane_2012[iQp];iplane=int(plane=='y');print plane,iplane
	Nb=Nb_2012[iQp];damp=damp_2012[iQp];
	# select Zxdip or Zydip
	for iw in imp_mod_2012:
	    if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);
	
	flag_trapz=0; # by default no trapz method
	if (M_2012==1): nxscan=np.array([0]);flag_trapz=1;
	else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M_2012,M_2012/20),np.arange(M_2012/2-10,M_2012/2+11),
		np.arange(M_2012-10,M_2012),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);
    
    	tuneshiftnx=np.zeros((1,len(nxscan),1,1,1,1,kmaxplot),dtype=complex);
    
	tuneshiftQp_2012[iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
	    nxscan,[damp],[Nb],[omegas],[dphase],M_2012,omega0,eval('Q'+plane),
	    gamma,eta,a,b,taub,g,Z,freq,particle='proton',flagnorm=0,flag_trapz=flag_trapz,
	    flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
	    abseps=1.e-4,flagm0=False,lxplusbatch=lxplusbatchDEL,
	    comment=machine+'_2012_v2_'+float_to_str(round(E_2012/1e9))+'GeV_'+str(M_2012)+'b_Qp'+str(Qp)+'_'+plane,
	    queue='2nd',dire=root_result+'/');

	# "stability factor" for each case
	factor_2012[iQp]=-oct_2012[iQp]*emit_2012[iQp]/(en2012**2*np.imag(tuneshiftQp_2012[iQp,0,0,0,0,0]));
	print "all factors 2012:",factor_2012;
	factor_neg_oct_2012_mean=np.average(factor_2012[ind_neg2012]);
	factor_pos_oct_2012_mean=np.average(factor_2012[ind_pos2012]);
	factor_neg_oct_2012_sig=np.sqrt(np.var(factor_2012[ind_neg2012]));
	factor_pos_oct_2012_sig=np.sqrt(np.var(factor_2012[ind_pos2012]));
	print "averages & sigmas:",factor_neg_oct_2012_mean,factor_neg_oct_2012_sig,factor_pos_oct_2012_mean,factor_pos_oct_2012_sig;

    # end of computations with 2012 model
    
    
    # scans on post LS1 model
    tuneshiftQp=np.zeros((len(sigmaIR3scan),len(sigmaIR7scan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=enLS1*1e12);
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");
    
    # build models
    for isig3,sig3 in enumerate(sigmaIR3scan):
    
	imp_mod_list=[]; # complete list of impedance scenarios
	sigmaIR7scan_str=[];
	sig3str='_IR3_'+float_to_str(sig3)+'sig';

	for isig7,sig7 in enumerate(sigmaIR7scan):

	    sigmaIR7scan_str.append("IR7 TCS-TCP: "+str(sig7)+" $ \sigma $ ");
	    sig7str='_IR7_'+float_to_str(sig7)+'sig';
	    settings_filename_coll=settings_filename_coll_prefix+float_to_str(sig3)+'sig_retractionIR7_'+float_to_str(sig7)+'sig.dat';

	    # compute total imp. model
	    imp_mod,wake_mod=LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,
		    settings_filename_coll,TDIcoating='preLS1',dire=path_here+"LHC_elements/",
		    commentcoll=scenario,direcoll='Coll'+scenario+'/',
		    lxplusbatch=lxplusbatchImp,BPM=BPMflag,beam=beam,squeeze=squeeze,
		    wake_calc=wake_calc)

	    imp_mod_list.append(imp_mod);


	    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

		tuneshiftQp[isig3,isig7,:,:,:,:,:,:,:,:]=DELPHI_wrapper(imp_mod,
		    Mscan,Qpscan,dampscan,Nbscan,[omegas],[dphase],omega0,Qx,Qy,gamma,eta,a,b,
		    taub,g,planes=['x','y'],nevery=nevery,particle='proton',flagnorm=0,
		    flagdamperimp=0,d=None,freqd=None,kmax=kmax,kmaxplot=kmaxplot,crit=5.e-2,
		    abseps=1.e-4,flagm0=False,lxplusbatch=lxplusbatchDEL,
		    comment=machine+scenario+sig3str+sig7str+'_'+float_to_str(round(E/1e9))+'GeV',
		    queue='2nd',dire=root_result+'/',flagQpscan_outside=True);

	# now the plots (outside loop on IR7 sigma scan)
	if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):

	    # plot and compare all scenarios
	    maxratio=plot_compare_imp_model(imp_mod_list[::2],sigmaIR7scan_str[::2],listcomp=['Zxdip','Zydip'],
		saveimp=root_result+'/plot_imp_'+machine+sig3str+"_IR7_sig_scenarios",
		saveratio=root_result+'/plot_imp_ratio_'+machine+sig3str+"_IR7_sig_scenarios",xlim=[1e3,1e12],legpos=[0.8,0.08]);

    if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	# parameters
	oct2012signscan=['_neg_oct','_pos_oct'];	    

	# post-LS1 scenarios
	#legbeam=['Classical 25 ns','BCMS 25 ns', 'Classical 50 ns','BCMS 50 ns'];
	#emitbeam=[3.75,1.9,2.2,1.6];
	#intbeam=[1.35,1.15,1.65,1.6];
	#Mbeam=[3564,3564,1782,1782];
	#colbeam=['ok','xk','oc','xc'];
	# from Giovanni Rumolo LBOC talk, 8/4/2014, with 0.6 mm.mrad emittance blow-up in LHC
	# except emittance standard 25ns (3.75 -> nominal)
	legbeam=['25 ns, standard','25 ns, BCMS','25ns, standard, 8b+4e','50 ns, standard (2012)'];
	strbeam=['_classic_25ns','_BCMS_25ns','_classic_25ns_8b4e','_classic_50ns'];
	emitbeam=[3.75,1.9,2.9,2.2];
	intbeam=[1.3,1.3,1.8,1.7];
	Mbeam=[3564,3564,3564,1782];
	colbeam=['ok','xk','+k','dc'];
	# HL-LHC scenarios
	#legbeam=['PIC low emit.','PIC high emit.', 'US1',"US2 low int. & low emit.","US2 high int. & high emit."];
	#emitbeam=[1.8,2.22,2.62,2.26,2.5];
	#intbeam=[1.38,1.38,1.9,1.9,2.2];
	#colbeam=['or','+m','xg','dk','vb'];

	colscen=['b','r','g','m','k','c','y'];
	hatch=np.array(["/","\\","|","X","/","\\","|","X"]);

	# plots stabilizing emittance vs Nb for certain Qp
	for isign,octsign in enumerate(oct2012signscan):

	    # assumes oct*emit/en^2 =  fact *imag(-tuneshift(int));
	    #factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
	    fact=eval('factor'+octsign+'_2012_mean');
	    sig_fact=eval('factor'+octsign+'_2012_sig');
	    # emit=fact*en^2/oct *imag(-tuneshift)

	    for iM,M in enumerate(Mscan):

		for idamp,damp in enumerate(dampscan):

		    # initialize plot
		    figall,axall=init_figure(axes=[0.15,0.1,0.8,0.7]);

		    # output file name for plots vs emittance
		    fileoutplotemitall=root_result+'/plot_int_vs_emit_'+machine+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';

		    # plot postLS1 beam parameters
		    for ib,legb in enumerate(legbeam):
	  		if (M==1)or(M==Mbeam[ib]):
			    plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",axall,0,xlab="Norm. emittance ($ \mu $ m)");

		    # sigma retraction scans
		    i=0;
		    for isig3,sig3 in enumerate(sigmaIR3scan):

			if sig3 in sigmaIR3scan_plot:

			    sig3str='_IR3_'+float_to_str(sig3)+'sig';

			    for isig7,sig7 in enumerate(sigmaIR7scan):

			    	if sig7 in sigmaIR7scan_plot:

				    sig7str='_IR7_'+float_to_str(sig7)+'sig';

				    # output file name for data int. vs emittance
				    Estr=float_to_str(round(E/1e9))+'GeV';
				    fileoutplotemit=root_result+'/plot_int_vs_emit_'+machine+'_'+Estr+scenario+sig3str+sig7str+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';
				    fileoutdataemit=root_result+'/data_int_vs_emit_'+machine+'_'+Estr+scenario+sig3str+sig7str+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';

				    # initialize plot with only one scenario and its "error bar" (spread)
				    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.75])
				    # plot postLS1 beam parameters
				    for ib,legb in enumerate(legbeam):
	  				if (M==1)or(M==Mbeam[ib]):
					    plot([emitbeam[ib]],[intbeam[ib]],legb,colbeam[ib],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");

				    imtu=np.zeros(len(Nbscan));
				    for iNb,Nb in enumerate(Nbscan):
					ts=getattr(tuneshiftQp[isig3,isig7,:,iM,iQpaver,idamp,iNb,0,0,0],'imag');
					if len(np.isnan(ts))==0: imtu[iNb]=np.min(ts);
					else: imtu[iNb]=np.min(ts[~np.isnan(ts)]);
		    			#imtu[iNb]=np.max(np.abs(tuneshiftQp[iscenario+1,:,iM,iQpaver,idamp,iNb,0,0,0]));
				    print imtu

				    emitLS1=fact*enLS1**2/octLS1 * np.abs(imtu);
				    sig_emitLS1=sig_fact*enLS1**2/octLS1 * np.abs(imtu)
				    #plot(emitLS1,Nbscan/1e11,"Limit with "+legscen[subscan[iscenario+1]],colscen[iscenario],"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Norm. emittance ($ \mu $ m)");
				    #ax.errorbar(emitLS1,Nbscan/1e11,xerr=sig_emitLS1,colscen[i],label="Stab. limit "+scenario+", retraction TCS-TCP: IR3 "+str(sig3)+" $ \sigma $, IR7 "+str(sig7)+" $ \sigma $ ",lw=3);
			    	    plot(emitLS1,Nbscan/1e11,"Stab. limit "+scenario+", retraction IR3: "+str(sig3)+" $ \sigma $, IR7 "+str(sig7)+" $ \sigma $ ",colscen[i],"Intensity ($ 10^{11} $p+/b)",axall,0,lw=3.,xlab="Norm. emittance ($ \mu $ m)");
			    	    
				    # plot of each individual line with "error bars"
			    	    plot(emitLS1,Nbscan/1e11,"Stab. limit "+scenario+", retraction IR3: "+str(sig3)+" $ \sigma $, IR7 "+str(sig7)+" $ \sigma $ ",colscen[i],"Intensity ($ 10^{11} $p+/b)",ax,0,lw=3.,xlab="Norm. emittance ($ \mu $ m)");
			    	    ax.fill_betweenx(Nbscan/1e11,emitLS1-sig_emitLS1,emitLS1+sig_emitLS1,color=colscen[i],lw=1,alpha=0.3)
			    	    #ax.set_xlabel("Norm. emittance ($ \mu $ m)");
			    	    #ax.set_ylabel("Intensity ($ 10^{11} $p+/b)");
				    ax.set_xlim([0,6]);ax.set_ylim([0,6]);
			    	    axall.set_xlim([0,6]);axall.set_ylim([0,6]);
				    i+=1;
			    	    end_figure(fig,ax,save=flagsave*fileoutplotemit,legpos=(0.3,0.8))

				    data=np.hstack((emitLS1.reshape((-1,1)),Nbscan.reshape((-1,1)),sig_emitLS1.reshape((-1,1))));
				    write_ncol_file(fileoutdataemit+'.dat',data,header="emit\tNb\tsig_emit")

		    # finish plot
		    end_figure(figall,axall,save=flagsave*fileoutplotemitall,legpos=(0.3,0.5),legfontsize=20)


	# plot lowest unstable intensity vs sigmas of retraction TCS-TCP in IR7, for all IR3 retraction scenarios
	for isign,octsign in enumerate(oct2012signscan):

	    # assumes oct*emit/en^2 =  fact *imag(-tuneshift(int));
	    #factor2012=oct2012*emit2012/(en^2*imag(-tuneshift));
	    fact=eval('factor'+octsign+'_2012_mean');
	    sig_fact=eval('factor'+octsign+'_2012_sig');
	    # emit=fact*en^2/oct *imag(-tuneshift)

	    for ib,emit in enumerate(emitbeam):

		for idamp,damp in enumerate(dampscan):

		    # initialize plot
		    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.75]);

		    # output file name for plots int vs sig retraction in IR7
		    fileoutplotint=root_result+'/plot_int_vs_sigIR7_'+machine+scenario+strbeam[ib]+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';

		    # plot postLS1 beam parameters
		    plot(sigmaIR7scan,intbeam[ib]*np.ones(len(sigmaIR7scan)),legbeam[ib],'--k',"Intensity ($ 10^{11} $p+/b)",ax,0,xlab="Retraction TCS-TCP in IR7 ($ \sigma $ )");

		    # sigma retraction scans
		    for isig3,sig3 in enumerate(sigmaIR3scan):

			sig3str='_IR3_'+float_to_str(sig3)+'sig';

			# output file name for data int. vs emittance
			Estr=float_to_str(round(E/1e9))+'GeV';
			fileoutdataint=root_result+'/data_int_vs_sigIR7_'+machine+'_'+Estr+scenario+sig3str+strbeam[ib]+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qpaver[0])+'_'+float_to_str(Qpaver[-1])+octsign+'_converged';

			Nbthres=np.zeros(len(sigmaIR7scan));
			sig_Nbthres=np.zeros(len(sigmaIR7scan));

			for isig7,sig7 in enumerate(sigmaIR7scan):

			    sig7str='_IR7_'+float_to_str(sig7)+'sig';

			    imtu=np.zeros(len(Nbscan));
			    for iNb,Nb in enumerate(Nbscan):
				ts=getattr(tuneshiftQp[isig3,isig7,:,iM,iQpaver,idamp,iNb,0,0,0],'imag');
		    		if len(np.isnan(ts))==0: imtu[iNb]=np.min(ts);
				else: imtu[iNb]=np.min(ts[~np.isnan(ts)]);

			    #emit=fact*en**2/oct * np.abs(imtu);
			    # growth rate at stability limit
			    imtu_lim=-emit*octLS1/(fact*enLS1**2);
			    sig_imtu=sig_fact*emit*octLS1/(fact**2*enLS1**2);
			    
			    # small test plot
			    figtmp,axtmp=init_figure();
			    axtmp.plot(imtu,Nbscan,'xb',label='original');
			    end_figure(figtmp,axtmp,save='plottmp_orig'+sig3str+sig7str+strbeam[ib]+'_d'+float_to_str(damp));
			    # end small test plot
			    
			    # interpolation to get corresponding intensity at the limit
			    # cubic spline interpolation
			    spl=interpolate.splrep(imtu[::-1],Nbscan[::-1]);
			    Nbthres[isig7]=interpolate.splev(imtu_lim,spl,der=0);
			    sig_Nbthres[isig7]=np.abs(interpolate.splev(imtu_lim,spl,der=1)*sig_imtu);
			    #print strbeam[ib],", d=",damp,octsign,", sig IR3=",sig3," sig IR7=",sig7,": Nb & sig_Nb [1e11 p+/b]:", Nbthres[isig7]/1e11,sig_Nbthres[isig7]/1e11;
			    
			    # test plots
			    # interpolation control plot
			    figtmp,axtmp=init_figure();
			    axtmp.plot(imtu,Nbscan,'xb',label='original');
			    imtu_ref=np.arange(imtu[0],imtu[-1],(imtu[-1]-imtu[0])/(10.*len(Nbscan)));
			    Nbscan_ref=interpolate.splev(imtu_ref,spl,der=0);
			    axtmp.plot(imtu_ref,Nbscan_ref,'-r',label='spline');
			    end_figure(figtmp,axtmp,save='plottmp'+sig3str+sig7str+strbeam[ib]+'_d'+float_to_str(damp));
			    
			    # derivative plot
			    figtmp,axtmp=init_figure();
			    axtmp.plot(imtu[:-1],np.diff(Nbscan)/np.diff(imtu),'xb',label='rough finite differences');
			    dNbscan_ref=interpolate.splev(imtu_ref,spl,der=1);
			    axtmp.plot(imtu_ref,dNbscan_ref,'-r',label='spline');
			    end_figure(figtmp,axtmp,save='plottmp_der'+sig3str+sig7str+strbeam[ib]+'_d'+float_to_str(damp));
			    # end test plots

			ax.errorbar(sigmaIR7scan,Nbthres/1e11,yerr=sig_Nbthres/1e11,fmt=colscen[isig3],label="Stab. limit, retraction IR3: "+str(sig3)+" $ \sigma $ ",lw=3,ms=15);
			ax.set_xlabel("Retraction TCS-TCP in IR7 ($ \sigma $ )");
			ax.set_ylabel("Intensity ($ 10^{11} $p+/b)");
			ax.set_xlim([2,6]);ax.set_ylim([0,6]);

			data=np.hstack((np.array(sigmaIR7scan).reshape((-1,1)),np.array(Nbthres).reshape((-1,1)),np.array(sig_Nbthres).reshape((-1,1))));
			write_ncol_file(fileoutdataint+'.dat',data,header="sigma_ret_IR7\tNb\tsig_Nb")

		    # finish plot
		    end_figure(fig,ax,save=flagsave*fileoutplotint,legpos=(0.3,0.6),legfontsize=20)

    if not(flagsave): pylab.show();
