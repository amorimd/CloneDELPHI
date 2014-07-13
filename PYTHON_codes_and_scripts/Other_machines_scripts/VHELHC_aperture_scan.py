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
from copy import deepcopy
import pylab,os,re
path_here=os.getcwd()+"/";
from plot_lib import plot,init_figure,end_figure,cmap
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_imp import *
from LHC_coll_imp import *

def VHELHC_param(E=50e12,Qxfrac=0.9,Qyfrac=0.9,V=16e6):

    e,m0,c,E0=proton_param();
    # E is the energy in eV
    Estr=str(int(E/1e12))+'TeV';print Estr
    machine='VHELHC';
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=100e3; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    sigmaz=8e-2; # from VHE-LHC preliminary parameters
    Qxint=120;Qyint=120; # based on beta~R/Q with beta proportional to E^(1/3) and taking as a reference the LHC case
    Qx=Qxint+Qxfrac;
    Qy=Qyint+Qyfrac;
    alphap=6.9e-5; # based on alphap=1/gamma_t^2 with gamma_t~Q
    
    h=circ/c/2.5e-9 # approximate harmonic number
    taub=4.*sigmaz/(beta*c); # full length in s
    eta=alphap-1./(gamma*gamma); # slip factor
    Qs=Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton');
    omegas=Qs*omega0;
    
    dphase=0.; # additional damper phase
    
    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    
    mu0=4e-7*np.pi;Z0=mu0*c;
    E=3e12;V=12e6;
    #E=50e12;V=64e6;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=VHELHC_param(E=E,V=V);
    circ=2*np.pi*R;
    avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

    col=['b','r','g','m','k','c','y']; # colors
    mark=['x','o','+','d'];
    pat=['-','--'];

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    #flagcompute=1; # 0 to avoid computing (simply plot from existing data)
    flagSach=(lxplusbatchDEL=='launch'); # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    #flagSach=0; # 1 to compute Sacherer tuneshifts, 0 to retrieve them
    wake_calc=False;
    
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=20; # number of plotted eigenvalues
    root_result=path_here+'../../../DELPHI_results/'+machine+'/vacpipe_'+Estr;
    os.system("mkdir -p "+root_result);
    suffix=''; # suffix for output files 
    lmaxSach=1;
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # scan definitions
    materialscan=['Cu20K','Cu50K']; # material is cold Cu (5K should be ~same as 20K)
    scenarioscan=deepcopy(materialscan);
    hgapscan=1e-3*np.array([5,7,8,9,10,11,12,13,14,15,20,25,30,35,40]);
    Nbscan=np.arange(1.e10,30.05e11,5.e9);
    Nbscanplot=np.array([2.e10,1.e11,2.e11]);
    dampscan=np.array([0,0.02,0.05,0.1]);
    Mscan=np.array([1,6672,13344,66720]);
    #Qpscan=np.arange(-10,21);
    #Qpplothg=np.array([0,2,5,10]); # scan in Qp for plot vs half-gap (and for TMCI plot)
    Qpscan=np.array([0]);
    Qpplothg=np.array([0]); # scan in Qp for plot vs half-gap (and for TMCI plot)

    imp_mod_list={}; # complete list of impedance scenarios
    wake_mod_list={};# complete list of wake scenarios
   
    for imat,material in enumerate(materialscan):
    
        name='vacpipe_'+material;
	imp_mod_list[material]=[];leg=[];

	for ihalfgap,halfgap in enumerate(hgapscan):

    	    imp_mod=[];wake_mod=[];
    
	    # compute model for vacuum pipe
	    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';
	    
	    # B field found from energy and radius of curvature (deduced from VHE-LHC parameters: 
	    # arc filling factor=0.78, 1400*12 m of straigth sections)
	    radius=(circ-1440*12)*0.78/2/np.pi;
	    layers_iw=[eval(material+'_layer(thickness=np.inf,RRR=70,B=E/1e9/(radius*0.299792458))')];
	    print material,", resistivity=",layers_iw[0].rhoDC,"Ohm.m,tau=",layers_iw[0].tau,'s';

	    betax=2*avbetax;betay=2*avbetay;

	    fpar=freq_param(fmin=1,fmax=1e14,ftypescan=2,nflog=40,fminrefine=1e11,fmaxrefine=5e12,nrefine=2000)
	    zpar=z_param(ztypescan=2,zmin=0.1,zmax=1e7,nzlog=20,zminrefine=2e-6,zmaxrefine=0.02,zsamplin=2e-6)
	    waketol=1.e12;
	    
	    iw_input=impedance_wake_input(gamma=gamma,length=1,b=[halfgap],layers=layers_iw,
    		fpar=fpar,zpar=zpar,waketol=waketol,freqlinbisect=1e11,geometry='round',comment='_'+name+strhgap);

	    imp_mod_vac,wake_mod_vac=imp_model_elliptic(iw_input,halfgap,orientation='V',
    		wake_calc=wake_calc,flagrm=True,lxplusbatch=lxplusbatchImp,queue='1nh',dire=machine+'vacpipe_'+Estr+'/');

	    multiply_impedance_wake(imp_mod_vac,circ);
	    multiply_impedance_wake(wake_mod_vac,circ);

	    # add to the model
	    add_impedance_wake(imp_mod,imp_mod_vac,betax/avbetax,betay/avbetay);
	    add_impedance_wake(wake_mod,wake_mod_vac,betax/avbetax,betay/avbetay);

	    if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):

	        imp_mod_list[material].append(imp_mod);
		leg.append('r='+str(halfgap*1e3)+'mm');
		
		# dump into a file
		filemodel=open(root_result+'/impedances'+name+strhgap+'.txt','w');
		pick.dump(imp_mod,filemodel);
		filemodel.close();
		
		# write Ascii files with each component
		write_imp_wake_mod(imp_mod,"_"+machine+"_Allthemachine_"+Estr+strhgap,
	    	    listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
		    dire=root_result+'/')
		
		# impedance plot
		plot_compare_imp_model([imp_mod],['Vacuum pipe RW impedance'],listcomp=['Zydip'],
			saveimp=root_result+'/plot_imp_'+name+strhgap,
			xlim=[1e3,1e11],ylim=[1e4,1e10],legpos=3);

	if (lxplusbatchImp.startswith('retrieve'))or(lxplusbatchImp==None):

	    # impedance list plot
	    plot_compare_imp_model(imp_mod_list[material][::3],leg[::3],listcomp=['Zydip'],
		    saveimp=root_result+'/plot_imp_'+name,
		    xlim=[1e3,1e11],ylim=[1e4,1e10],legpos=3);

    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('retrieve')):
	
	# DELPHI loops now
	tuneshiftQp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
	tuneshiftm0Qp=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),1,1),dtype=complex);
	tuneshiftQpSach=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1,(2*lmaxSach+1)),dtype=complex);
	tuneshiftQpSachm0=np.zeros((len(hgapscan),len(scenarioscan),2,len(Mscan),len(Qpscan),len(Nbscan),1),dtype=complex);

	for ihalfgap,halfgap in enumerate(hgapscan):
	
	    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

	    for iscenario,scenario in enumerate(scenarioscan):

        	imp_mod=imp_mod_list[scenario][ihalfgap];

		for iplane,plane in enumerate(['x']):
		    # select Zxdip or Zydip
		    for iw in imp_mod:
			if test_impedance_wake_comp(iw,1-iplane,iplane,0,0,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);
		    
		    for iM,M in enumerate(Mscan):

			# normalization factor for damper
			dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
    			    flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
			dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

			flag_trapz=0; # by default no trapz method

        		if (M==1): nxscan=np.array([0]);flag_trapz=1;
			else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,M/20),np.arange(M/2-10,M/2+11),
	    		    np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

			tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),1,1,kmaxplot),dtype=complex);
			tuneshiftnxSach=np.zeros((len(Qpscan),len(nxscan),len(Nbscan),1,2*lmaxSach+1),dtype=complex);
			ZeffSach=np.zeros((len(Qpscan),len(nxscan),1,2*lmaxSach+1),dtype=complex);

			# DELPHI
			print "half-gap=",halfgap,", M=",M
			tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
				nxscan,dampscan,Nbscan,[omegas],[dphase],M,omega0,eval('Q'+plane),gamma,eta,
				a,b,taub,g,Z,freq,particle='proton',flagnorm=flagnorm,
				flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,freqd=None,
				kmax=kmax,kmaxplot=kmaxplot,crit=1.e-2,abseps=1e-5,flagm0=True,
				lxplusbatch=lxplusbatchDEL,comment=machine+strhgap+scenario+'_'+Estr+str(M)+'b'+strnorm[flagnorm]+'_'+plane,
				queue='1nd',dire=root_result+'/');

			if flagSach:
			    # Sacherer (no damper)
			    tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,:,:,:],tuneshiftnxSach,tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,:,:],ZeffSach=sacherer(imp_mod,
				    Qpscan,nxscan,Nbscan,[omegas],M,omega0,eval('Q'+plane),gamma,eta,taub,lmaxSach,
				    particle='proton',modetype='sinusoidal',compname='Z'+plane+'dip');

	if flagSach:
	    # save Sacherer tuneshifts
	    fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','w');
	    pick.dump(tuneshiftQpSach,fileSach);
	    pick.dump(tuneshiftQpSachm0,fileSach);
	    fileSach.close();
	else:
	    # load Sacherer tuneshifts
	    fileSach=open(root_result+'/Sacherer_'+Estr+'.txt','r');
	    tuneshiftQpSach=pick.load(fileSach);
	    tuneshiftQpSachm0=pick.load(fileSach);
	    fileSach.close();
	

	# now the plots (outside loop on scenarios)
	if (lxplusbatchDEL==None)or(lxplusbatchDEL.startswith('retrieve')):

	    for iplane,plane in enumerate(['x']):

		for iM,M in enumerate(Mscan):

		    for idamp,damp in enumerate(dampscan):

			for Nb in Nbscanplot:

	        	    if False:
				# plots vs Q'
				for ihalfgap,halfgap in enumerate(hgapscan):

				    strhgapleg='halfgap '+str(1e3*halfgap)+'mm';
				    strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

				    for iscenario,scenario in enumerate(scenarioscan):

					# initialize plots vs Qp
					figQpm0,axQpm0=init_figure(axes=[0.15,0.1,0.8,0.85]);
					figQp=[];axQp=[];
					for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

					# output file name for plots vs Qp
					fileoutplotQp=root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
					fileoutplotQpm0=root_result+'/plot_vs_Qp_m0_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

					strpart=['Re','Im'];
					for ir,r in enumerate(['real','imag']):

					    # output file name for data vs Qp
					    fileoutdataQp=root_result+'/data_vs_Qp_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
					    fileoutdataQpm0=root_result+'/data_vs_Qp_m0_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

					    ts=getattr(tuneshiftQp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
					    Sachstr='';
					    if damp==0:
						# compare with Sacherer most unstable mode
						tsSach=getattr(tuneshiftQpSach[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0,0],r);
						Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
						data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
					    else:
						data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));

					    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift"+Sachstr)

					    sgn=1;sgnstr='';
					    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
					    plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");
					    if damp==0:
						plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQp[ir],0,xlab=" $ Q^' $ ");


					    # real tune shift of mode 0
					    if (ir==0):
						ts=getattr(tuneshiftm0Qp[ihalfgap,iscenario,iplane,iM,:,idamp,pylab.mlab.find(Nbscan==Nb),0,0],r);
						Sachstr='';

						if damp==0:
						    # compare with Sacherer most unstable mode
						    tsSach=getattr(tuneshiftQpSachm0[ihalfgap,iscenario,iplane,iM,:,pylab.mlab.find(Nbscan==Nb),0],r);
						    Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift_mode0"
					    	    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));

						else:
						    data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));

						write_ncol_file(fileoutdataQpm0+'_'+r+'.dat',data,header="Qp\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr)

						sgn=1;sgnstr='';
						if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
						plot(Qpscan,np.squeeze(sgn*ts),'DELPHI, '+strhgapleg+', '+scenario,col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");
						if damp==0:
					    	    plot(Qpscan,np.squeeze(sgn*tsSach),'Sacherer, '+strhgapleg+', '+scenario,'--'+col[iscenario],"$ "+sgnstr+strpart[ir]+"(Q-Q_0) $ ",axQpm0,0,xlab=" $ Q^' $ ");

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
					    end_figure(figQpm0,axQpm0,save=flagsave*(fileoutplotQpm0+'_'+r))

				    for ir,r in enumerate(['real','imag']):
					end_figure(figQp[ir],axQp[ir],save=flagsave*(fileoutplotQp+'_'+r))

	        	    
			    # plot imag tune shifts (in terms of growth rates) of most unstable mode vs half-gap
			    for Qp in Qpplothg:

				iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];
				
				strpart=['Re','Im'];
				r='real';ir=0;fact=1;ylab=" $ |"+strpart[ir]+"(Q-Q_0)| $ ";
				r='imag';ir=1;fact=omega0;ylab="Growth rate [1/s]";

				for iscenario,scenario in enumerate(scenarioscan):

				    # initialize plots vs half-gap
				    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);

				    # output file name for plots vs hg
				    fileoutplothg=root_result+'/plot_vs_hg_m0_'+machine+'_'+Estr+'_Qp'+float_to_str(Qp)+'_'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    # output file name for data vs halfgap
				    fileoutdatahg=root_result+'/data_vs_hg_m0_'+machine+'_'+Estr+'_Qp'+float_to_str(Qp)+'_'+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				    # tune shift of most unstable mode
				    ts=getattr(tuneshiftQp[:,iscenario,iplane,iM,iQp,idamp,pylab.mlab.find(Nbscan==Nb),0,0,0],r);
				    Sachstr='';
				    if damp==0:
					# compare with Sacherer mode 0
					tsSach=getattr(tuneshiftQpSach[:,iscenario,iplane,iM,iQp,pylab.mlab.find(Nbscan==Nb),0,0],r);
					Sachstr=" Sacherer_"+strpart[ir]+"_tuneshift"
					data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1)),tsSach.reshape((-1,1))));
				    
				    else:
				    	data=np.hstack((hgapscan.reshape((-1,1)),ts.reshape((-1,1))));
				    
				    write_ncol_file(fileoutdatahg+'_'+r+'.dat',data,header="halfgap[m]\tDELPHI_"+strpart[ir]+"_tuneshift_mode0"+Sachstr)

				    sgn=1;sgnstr='-';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    
				    plot(hgapscan*1e3,np.abs(np.squeeze(sgn*ts))*fact,'DELPHI',col[iscenario],ylab,ax,2,xlab=" half-gap [mm] ");
				    if damp==0:
				    	plot(hgapscan*1e3,np.abs(np.squeeze(sgn*tsSach))*fact,'Sacherer','--'+col[iscenario],ylab,ax,2,xlab=" half-gap [mm] ");

				    # add on top of the plot some lines for the feedback gain (100, 50 and 10 turns)
				    for idamping,damping in enumerate([10.,20.,50.]):
					dampstr=str(int(damping))+' turns damping';
					plot(hgapscan*1e3,f0/damping*np.ones(len(hgapscan)),dampstr,'-'+col[iscenario+idamping+1],ylab,ax,2,xlab=" half-gap [mm] ");
				    
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
				    end_figure(fig,ax,save=flagsave*(fileoutplothg+'_'+r))

			# TMCI plots
			Nbthres=np.zeros((len(hgapscan),len(scenarioscan),len(Qpplothg)))
			for Qp in Qpplothg:

			    iQp=pylab.mlab.find(Qpscan==Qp);iQp=iQp[0];

			    for ihalfgap,halfgap in enumerate(hgapscan):

				strhgap='_halfgap_'+float_to_str(round(1e4*halfgap)/10.)+'mm';

				for iscenario,scenario in enumerate(scenarioscan):

				    fileoutplotTMCI=root_result+'/plot_TMCI_'+machine+'_'+Estr+strhgap+scenario+'_'+str(M)+'b_d'+float_to_str(damp)+'_Qp'+float_to_str(round(10.*Qp)/10.)+'_converged'+strnorm[flagnorm]+'_'+plane;
				    patcol=['.b','b'];
				    ylim=([-5,3],[-0.001,0.02]);

		    		    for ir,r in enumerate(['real','imag']):

					fig,ax=init_figure();

					ts=tuneshiftQp[ihalfgap,iscenario,iplane,iM,iQp,idamp,:,0,0,:];

					plot_TMCI(Nbscan,ts/Qs,ax,part=r,leg='DELPHI',patcol=patcol[ir],xlab='Nb [p+/b]',
					    title=machine+r", $ Q^' = $ "+str(round(100*Qp)/100.)+", half-gap="+str(halfgap*1e3)+"mm, "+scenario,ms=1,ylim=ylim[ir]);

					end_figure(fig,ax,save=flagsave*(fileoutplotTMCI+'_'+r),fontsize=25);

				    # initialize plot of threshold vs half-gap
				    fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);
				    # output file name for threshold plot vs hg
				    fileoutplothg=root_result+'/plot_vs_hg_thres_'+machine+'_'+Estr+scenario+'_Qp'+float_to_str(Qp)+'_'+str(M)+'b_d'+float_to_str(damp)+'_converged'+strnorm[flagnorm]+'_'+plane;

				    # find intensity threshold
				    Nbthres[ihalfgap,iscenario,iQp]=find_intensity_threshold(Nbscan,tuneshiftQp[ihalfgap,iscenario,iplane,iM,iQp,idamp,:,0,0,0]*omega0,thresgrowth=10.);
				    plot(hgapscan*1e3,np.squeeze(Nbthres[:,iscenario,iQp]),'DELPHI, '+scenario,col[iscenario],"TMCI threshold (nb p/b)",ax,2,xlab=" half-gap [mm] ");
				    ax.set_xlim([5,20]);ax.set_ylim([1e10,2e12]);
				
				    # finish plots vs half-gap
				    end_figure(fig,ax,save=flagsave*(fileoutplothg+'_'+r))

    if not(flagsave): pylab.show();
