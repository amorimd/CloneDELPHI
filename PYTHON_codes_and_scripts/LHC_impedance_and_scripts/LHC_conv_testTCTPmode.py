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
from io_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param


def TCTP_modes(scenario=0):

    # retrieve the TCTP low freq. modes parameters (resonance frequencies, shunt impedances
    # and quality factor), and the collimator settings (half-gaps and beta functions)
    # scenario is the TCT settings & beta functions scenario (from Roderik Bruce):
    #	0 = first realistic (small beta*),
    #	1 = second realistic (large beta*),
    #	2 = first pessimistic (8 sigmas, beta*=0.55 m),
    #	3 = second pessimistic (12 sigmas, beta*=1.5 m).
    
    # fr (Hz), R (Ohm/m) and Q from CST for 5 mm half-gap (correct materials) (from Benoit Salvant)
    fmodes=np.array([0.099,0.201,0.307,0.514,0.607])*1.e9;
    Rmodes=np.array([160.,51.,11.,11.,11.])/3.e-4;
    Qmodes=np.array([27.5,21.847826087,15.35,25.7,30.35]);

    # R vs. half-gap (from Benoit, CST simulation with PEC, normalized 
    # to unity at a halg-gap of 5mm)
    halfgap=np.arange(2,11); # in mm
    Rnorm_halfgap=np.array([637.,350.,192.,112.,75.,50.,40.,30.,25.])/112.;
    
    # read settings
    root='../Impedances/PostLS1/TCTs/TCT_postLS1_';
    suffix='_RBruce.dat';
    names=['realistic_scenario_small_betastar','realistic_scenario_large_betastar',
    	'pessimistic_scenario_8sigmas_betastar0p55','pessimistic_scenario_12sigmas_betastar1p5'];
    filename=root+names[scenario]+suffix;
    modestr=names[scenario];
    settTCT=read_ncol_file(filename,ignored_rows=1,ignored_cols=1);
    # extract TCT names
    namesTCT=[];
    for il,line in enumerate(open(filename)):
    	if (il>0): namesTCT.append(split(line)[0]);
	
    
    return fmodes,Rmodes,Qmodes,halfgap,Rnorm_halfgap,namesTCT,settTCT,modestr;
	

def TCTP_modes_impedance(betaavx,betaavy,nmodes=1,scenario=0,flagplot=True):

    # compute TCTP low freq. modes impedance
    # betaavx and betaavy are the average betas taken for the beta function weighting (usually circ/(2*pi*tune) ),
    # nmodes is the number of resonances taken into account (up to 5), scenario is the
    # TCT settings & beta functions scenario (from Roderik Bruce):
    #	0 = first realistic (small beta*),
    #	1 = second realistic (large beta*),
    #	2 = first pessimistic (8 sigmas, beta*=0.55 m),
    #	3 = second pessimistic (12 sigmas, beta*=1.5 m).
    # plot total modes impedance only when flagplot is True.
    
    fmodes,Rmodes,Qmodes,halfgap,Rnorm_halfgap,namesTCT,settTCT,modestr=TCTP_modes(scenario=scenario);
    
    # frequency span
    f=np.concatenate((10.**np.arange(-1,5),np.arange(1.e5,0.7001e9,1.e5),10.**np.arange(9.1,13.1,0.1),
    	10.**np.arange(14,15)));#print len(f)
    Zx=np.zeros((len(f),2));Zy=np.zeros((len(f),2));
    
    for imode in range(nmodes):
    
        for iTCT,TCTname in enumerate(namesTCT):
	
	    hg=settTCT[iTCT,0];beta=settTCT[iTCT,1];
	    Rmode=np.interp(hg,halfgap,Rnorm_halfgap)*Rmodes[imode]*beta;
	    
	    if (TCTname.startswith('TCTH')):
	    	Zx=Zx+resonator_impedance(Rmode/betaavx,fmodes[imode],Qmodes[imode],f);
	    	#print TCTname,", fr=",fmodes[imode],", Rmode=",Rmode/betaavx;
	    elif (TCTname.startswith('TCTV')):
	    	Zy=Zy+resonator_impedance(Rmode/betaavy,fmodes[imode],Qmodes[imode],f);
	    	#print TCTname,", fr=",fmodes[imode],", Rmode=",Rmode/betaavy;

    # plot
    if flagplot:
	for plane in ['x','y']:
	    Z=eval('Z'+plane);
	    figr,axr=init_figure()
	    plot(f,Z[:,0],'real part','b',"Z [ $ \Omega/ $m]",axr,0,xlab='Frequency [Hz]');
	    plot(f,Z[:,1],'imag. part','r',"Z [ $ \Omega/ $m]",axr,0,xlab='Frequency [Hz]');
	    axr.set_xlim([1e7,7e8]);
    	    end_figure(figr,axr,save=path_here+'../../../DELPHI_results/'+machine+'/plot_resonator_impedance_'+str(nmodes)+'modes_'+modestr)

    return f,Zx,Zy,modestr;

    
if __name__ == "__main__":

    e,m0,c,E0=proton_param();

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=7000e9);

    strnorm=['','_norm_current_chroma'];
    os.system("mkdir -p ../../../DELPHI_results/"+machine);

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
    flagcompute=0; # 0 to avoid computing (simply plot from existing data)
       
    kmax=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    kmaxplot=1; # number of converged eigenvalues (kmax most unstable ones are converged)
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

    # TCTP mode(s)
    #modestr=float_to_str(Rs/1e6)+'MOhmm_fr'+float_to_str(fr/1e6)+'MHz_Q'+float_to_str(Q);
    #modestrlong="$ R="+str(Rs/1e6)+"$ M$\Omega/ $m, $ f_r=$ "+str(fr/1e6)+" MHz, $ Q= $ "+str(Q);

    # MOSES results
    #fileMOSES='../MOSES_results/Results_'+machine+'/Rs'+float_to_str(Rs/1e6)+'e6_fr'+float_to_str(fr/1e6)+'MHz_Q'+float_to_str(Q)+'/'+machine+'_Qpscan';
    #print fileMOSES;

    # scan definition
    Qpscan=np.arange(-10,21,2);
    dampscan=np.array([0,0.02]); # damper gain scan
    #Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.5e11]); # intensity scan
    #Nbscan=np.array([]);dampscan=np.array([]);
    models=['baseline','relaxed','veryrelaxed'];
    #modelstr=['nominal','tight','relaxed']; # correct names
    
    for scenario in range(4):

	f,Zxmode,Zymode,modestr=TCTP_modes_impedance(R/Qx,R/Qy,nmodes=5,scenario=scenario);
	modestrlong=modestr.replace('_',' ').replace('p5','.5').replace('sigmas'," $ \\sigma $").replace('betastar'," $ \\beta^* $ ");
    	print modestr,modestrlong

	for iplane,plane in enumerate(['x','y']):

	    for imodel,model in enumerate(models):

		#if (imodel>0):
        	freq,Z=readZ('../Impedances/PostLS1/total_impedances/Z'+plane+plane+'dip_Allthemachine_7TeV_B1_postLS1_'+model+'.dat')
		modescan=['','_withTCTPmode_'+modestr];
		leg=['without TCTP mode','with TCTP mode, '+modestrlong];
    		Mscan=np.array([1,1782,3564]); # scan on number of bunches
		#else:
        	#    freq=f;Z=Zmode;
		#    modescan=['_withTCTPmode_'+modestr];leg=['with TCTP mode only, '+modestrlong];
    		#    Mscan=np.array([1]); # scan on number of bunches

    		# figure, axis and output file name for impedance plots
		figZ,axZ=init_figure()
		fileoutZ=path_here+'../../../DELPHI_results/'+machine+'/plot_imp_'+machine+'_'+model+'_'+modestr+'_'+plane;

		for iM,M in enumerate(Mscan):

		    # normalization factor for damper
		    dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
    			flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
		    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

		    flag_trapz=0; # by default no trapz method

        	    if (M==1): nxscan=np.array([0]);flag_trapz=1;
		    #elif (M==1782): nxscan=np.array([0, 1, 300, 600, 880, 890, 891, 892, 900, 910, 950, 1000, 1200, 1500, 1780, 1781])
		    #elif (M==3564): nxscan=np.array([0, 1, 300, 600, 900, 1200, 1500, 1770, 1780, 1781, 1782, 1785, 1790, 1800, 1900, 2000, 2300, 2600, 2900, 3200, 3500, 3560, 3561, 3562, 3563])
		    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,10),np.arange(M/2-20,M/2+21),
	    		np.arange(M-20,M),np.arange(0,20))));print "number of coupled-bunch modes=",len(nxscan);

		    tuneshiftQp=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(modescan)),dtype=complex);

		    for imode,mode in enumerate(modescan):

        		Zmode=eval('Z'+plane+'mode');
			# total impedance = imode * (resonator mode) + impedance model
    			freqtot=sort_and_delete_duplicates(np.concatenate((f,freq[::10])),tolerance=0.1);
			Ztot=np.zeros((len(freqtot),2),dtype=float);
			for icol in range(2): Ztot[:,icol]=imode*np.interp(freqtot,f,Zmode[:,icol])+np.interp(freqtot,freq,Z[:,icol]);
			print len(freqtot)

			if (iM==0):
			    # plot impedance
			    strpart=['Re','Im'];colpart=['b','r'];
			    for ir,r in enumerate(['real','imag']):
	    			plot(freqtot,Ztot[:,ir],r+" part, "+leg[imode],linetype[ir]+col[imode],"Z [ $ \Omega/ $m] ",axZ,3,xlab='Frequency [Hz]');
				axZ.set_xlim([1e3,1e11]);
				axZ.set_ylim([1e5,1e9]);

			for iQp,Qp in enumerate(Qpscan):

			    tuneshiftnx=np.zeros((len(nxscan),len(dampscan),len(Nbscan),kmaxplot),dtype=complex);
			    lmaxi=-1;nmaxi=-1; # to keep track of maximum lmax and nmax for all modes

			    for inx,nx in enumerate(nxscan):

    				omegaksi=Qp*omega0/eta;
				#print "Qp=",Qp,", chromatic freq. (not angular)=",omegaksi/(2.*np.pi);

				if flagnorm:
				    # normalization factor for damper
				    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    				flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
				    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

   				lmax=-1;nmax=-1;matZ=None;matdamper=None;

    				if flagcompute:
				
				    # computation
				    for idamp,damp in enumerate(dampscan):

					#lambdax=np.zeros((len(Nbscan),kmaxplot),dtype=complex);

					for iNb,Nb in enumerate(Nbscan):

					    #print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;

					    coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,eval('Q'+plane),particle='proton');

					    freqshift,v,lmax,nmax,matdamper,matZ=eigenmodesDELPHI_converged(nx,
						    M,omegaksi,omega0,eval('Q'+plane+'frac'),a,b,taub,g,Ztot,freqtot,coefdamper,
						    coefZ,omegas,flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=None,
						    freqd=None,kmax=kmax,crit=5.e-2,abseps=1.e-3,lmaxold=lmax,
						    nmaxold=nmax,matdamperold=matdamper,matZold=matZ);

					    tuneshiftnx[inx,idamp,iNb,:]=freqshift[:kmaxplot]/omega0;
					    #lambdax[iNb,:]=freqshift[:kmaxplot]/omegas;
					    lmaxi=max(lmax,lmaxi);nmaxi=max(nmax,nmaxi);	


    			    # finding the most unstable coupled-bunch mode
			    for idamp,damp in enumerate(dampscan):

				for iNb,Nb in enumerate(Nbscan):

				    inx=np.argmin(np.imag(tuneshiftnx[:,idamp,iNb,0]))
				    print plane,model,mode,", M=",M,", Qp=",Qp,", d=",damp,", Nb=",Nb;
				    print "   lmaxi=",lmaxi,", nmaxi=",nmaxi,", Most unstable coupled-bunch mode: ",nxscan[inx];
				    tuneshiftQp[iQp,idamp,iNb,imode]=tuneshiftnx[inx,idamp,iNb,0];


		    # plots vs Q'
		    for idamp,damp in enumerate(dampscan):

			for iNb,Nb in enumerate(Nbscan):

			    # initialize plots vs Qp
			    figQp=[];axQp=[];
			    for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

	    		    for imode,mode in enumerate(modescan):

				# output file name for plots vs Qp
				fileoutplotQp=path_here+'../../../DELPHI_results/'+machine+'/plot_vs_Qp_'+machine+'_'+str(round(E/1e9))+'GeV_'+model+'_'+modestr+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;
				# output file name for data vs Qp
				fileoutdataQp=path_here+'../../../DELPHI_results/'+machine+'/data_vs_Qp_'+machine+'_'+str(round(E/1e9))+'GeV_'+model+mode+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane;

				strpart=['Re','Im'];
				for ir,r in enumerate(['real','imag']):

				    if flagcompute:
					ts=getattr(tuneshiftQp[:,idamp,iNb,imode],r);
					data=np.hstack((Qpscan.reshape((-1,1)),ts.reshape((-1,1))));
					write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header="Qp\t"+strpart[ir]+"_tuneshift")
				    else:
				    	s=read_ncol_file(fileoutdataQp+'_'+r+'.dat',ignored_rows=1);
					Qpscan=s[:,0];ts=s[:,1];
					
				    sgn=1;sgnstr='';
				    if (ir==1): sgn=-1;sgnstr='-'; # invert sign of imaginary part
				    plot(Qpscan,sgn*ts/Qs,'DELPHI, '+leg[imode],col[imode],"$ "+sgnstr+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");


				    #if ((damp==0)and(M==1))and(imodel==0):
				#	# plot MOSES results as well in this case
				#	s=read_ncol_file(fileMOSES+'_'+float_to_str(Nb/1e11)+'e11.dat',ignored_rows=1);
				#	#sgn=1;
				#	#if (ir==1): sgn=-1; # invert sign of imaginary part
				#	plot(s[:,0],sgn*s[:,ir+1],'MOSES, '+leg[imode],'xg',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

				    if (scenario==2)and((M==1)and(imodel==0)):
					# plot HEADTAIL results as well in this case
					#fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_lin_dip'
					#rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin1_drate"+float_to_str(damp)+"_dip";
					fileoutplotQp=fileoutplotQp+'_vs_HEADTAIL_nonlin_all'
					rootHEADTAIL="/afs/cern.ch/work/n/nmounet/private/DATA_HEADTAIL/LHC_with_damper/LHC_testTCTPmodes/LHC_damper_1b_ntwake20_nkick1_nsl500_npr1000000_I1p5_qsec0_oct0_baseline_nlin4_drate"+float_to_str(damp);
					sufHEADTAIL="_aver_Sussix_most_tau.txt";
					s=read_ncol_file(rootHEADTAIL+mode+sufHEADTAIL,ignored_rows=1);
					fact=1;
					if (ir==1): fact=1/omega0; # for imaginary part, divide by omega0
					plot(s[:,0],fact*s[:,3*iplane+ir+1]/Qs,'HEADTAIL, '+leg[imode],'x'+col[imode+2],"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			    # finish plot vs Qp
			    for ir,r in enumerate(['real','imag']):
				if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r,fontsize=20)
				else: end_figure(figQp[ir],axQp[ir]);


		    # finish impedance plot (after the loop scanning with / without TCTP mode)
		    if flagsave: end_figure(figZ,axZ,save=fileoutZ)
		    else: end_figure(figZ,axZ);


    if not(flagsave): pylab.show();
