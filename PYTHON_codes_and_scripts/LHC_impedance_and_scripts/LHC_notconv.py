#!/usr/bin/python

import sys
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

    
if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    E=4e12;

    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr=LHC_param(E0,E=E);
    #eta=4.513e-4 # Xavier's case

    beam='1';
    # directory (insde DELPHI_results/[machine]) where to put the results
    root_result='/test_Xavier';os.system("mkdir -p "+root_result);
    # directory where are the impedances
    dirZ='/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/';
    
    lmax=10;nmax=10; # max azimuthal and radial mode numbers
    suffix='_lmax'+str(lmax)+'_nmax'+str(nmax); # suffix for output files

    strnorm=['','_norm_current_chroma'];

    flagsave=1; # 1 to save figure instead of plotting on screen
    flagdamperimp=0; # 1 to use frequency dependent damper gain (provided in Zd,fd)
    flagnorm=0; # 1 if damper matrix normalized at current chromaticity (instead of at zero chroma)
       
    col=['b','r','g','m','k','c','y']; # colors
    linetype=['-','--',':'];
    kmaxplot=1; # number of plotted eigenvalues
    
    # longitudinal distribution initialization
    g,a,b=longdistribution_decomp(taub,typelong="Gaussian");


    # MOSES results
    #fileMOSES='../MOSES_results/Results_'+machine+'/Rs'+float_to_str(Rs/1e6)+'e6_fr'+float_to_str(fr/1e6)+'MHz_Q'+float_to_str(Q)+'/'+machine+'_Qpscan';
    #print fileMOSES;

    # scan definition
    Qpscan=np.arange(-30,31,1);
    dampscan=np.array([0]); # damper gain scan
    #Nbscan=np.arange(1.e10,8.1e11,1.e10); # intensity scan
    Nbscan=np.array([1.5e11]); # intensity scan
    #Nbscan=np.array([]);dampscan=np.array([]);
    #models=['baseline','relaxed','veryrelaxed'];
    models=['_physics_fill_3265'];
    Mscan=np.array([1]); # scan on number of bunches
    
    for iplane,plane in enumerate(['x']):

	for imodel,model in enumerate(models):

           #freq,Z=readZ('../Impedances/PostLS1/total_impedances/Z'+plane+'dip_Allthemachine_7TeV_B1_postLS1_'+model+'.dat')
            freq,Z=readZ(dirZ+'Z'+plane+plane+'dip_Allthemachine_'+Estr+'_B'+beam+model+'.dat')
	    freqtot=freq[::10];Ztot=Z[::10,:];

    	    # figure, axis and output file name for impedance plots
	    figZ,axZ=init_figure()
	    fileoutZ=path_here+'../../../DELPHI_results/'+machine+root_result+'/plot_imp_'+machine+model+'_'+plane;

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

		tuneshiftQp=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),kmaxplot),dtype=complex);

		if (iM==0):
		    # plot impedance
		    strpart=['Re','Im'];colpart=['b','r'];
		    for ir,r in enumerate(['real','imag']):
	    		plot(freqtot,Ztot[:,ir],r+" part",linetype[ir],"Z [ $ \Omega/ $m] ",axZ,3,xlab='Frequency [Hz]');
			axZ.set_xlim([1e3,1e11]);
			axZ.set_ylim([1e5,1e9]);

		for iQp,Qp in enumerate(Qpscan):

		    tuneshiftnx=np.zeros((len(nxscan),len(dampscan),len(Nbscan),kmaxplot),dtype=complex);

		    for inx,nx in enumerate(nxscan):

    			omegaksi=Qp*omega0/eta;
			print "Qp=",Qp,", chromatic freq. (not angular)=",omegaksi/(2.*np.pi);

			if flagnorm:
			    # normalization factor for damper
			    dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    			flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-4);
			    dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

   			matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
				flagdamperimp=flagdamperimp,d=None,freqd=None,abseps=1e-3);

			matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,eval('Q'+plane+'frac'),a,b,taub,g,
				Ztot,freqtot,flag_trapz=flag_trapz,abseps=1e-3);

    			# computation
			for idamp,damp in enumerate(dampscan):

			    #lambdax=np.zeros((len(Nbscan),kmaxplot),dtype=complex);

			    for iNb,Nb in enumerate(Nbscan):

				print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;

				coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,eval('Q'+plane),particle='proton');

				freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas)
				
				# kmaxplot most unstable mode
				freqshift=freqshift.reshape(-1);
				ind=np.argsort(np.imag(freqshift));
				
				tuneshiftnx[inx,idamp,iNb,:]=freqshift[ind[:kmaxplot]]/omega0;
				#lambdax[iNb,:]=freqshift[:kmaxplot]/omegas;


    		    # finding the most unstable coupled-bunch mode
		    for idamp,damp in enumerate(dampscan):

			for iNb,Nb in enumerate(Nbscan):

			    for kmode in range(kmaxplot):
			    
				inx=np.argmin(np.imag(tuneshiftnx[:,idamp,iNb,kmode]))
				print plane,model,", M=",M,", Qp=",Qp,", d=",damp,", Nb=",Nb;
				print "kmode=",kmode,", Most unstable coupled-bunch mode: ",nxscan[inx];
				tuneshiftQp[iQp,idamp,iNb,kmode]=tuneshiftnx[inx,idamp,iNb,kmode];


		# plots vs Q'
		for idamp,damp in enumerate(dampscan):

		    for iNb,Nb in enumerate(Nbscan):

			# initialize plots vs Qp
			figQp=[];axQp=[];
			for ir in range(2): fig,ax=init_figure(axes=[0.15,0.1,0.8,0.85]);figQp.append(fig);axQp.append(ax);

			# output file name for plots vs Qp
			fileoutplotQp=path_here+'../../../DELPHI_results/'+machine+root_result+'/plot_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+model+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+suffix;
			# output file name for data vs Qp
			fileoutdataQp=path_here+'../../../DELPHI_results/'+machine+root_result+'/data_vs_Qp_'+machine+'_'+float_to_str(round(E/1e9))+'GeV'+model+'_'+str(M)+'b_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1.e11)+'e11_converged'+strnorm[flagnorm]+'_'+plane+suffix;

			strpart=['Re','Im'];
			for ir,r in enumerate(['real','imag']):

			    ts=getattr(tuneshiftQp[:,idamp,iNb,:],r);

			    plot(Qpscan,ts/Qs,'','',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			    data=np.hstack((Qpscan.reshape((-1,1)),ts));
			    # column headers
			    head="Qp";
			    for kmode in range(kmaxplot): head+="\t"+strpart[ir]+"_tuneshift_mode"+str(kmode);
			    
			    write_ncol_file(fileoutdataQp+'_'+r+'.dat',data,header=head)

			    #if ((damp==0)and(M==1))and(imodel==0):
			#	# plot MOSES results as well in this case
			#	s=read_ncol_file(fileMOSES+'_'+float_to_str(Nb/1e11)+'e11.dat',ignored_rows=1);
			#	sgn=1;
			#	if (ir==1): sgn=-1; # invert sign of imaginary part
			#	plot(s[:,0],sgn*s[:,ir+1],'MOSES, '+leg[imode],'xg',"$ "+strpart[ir]+"(Q-Q_0)/Q_s $ ",axQp[ir],0,xlab=" $ Q^' $ ");

			# finish plot vs Qp
			for ir,r in enumerate(['real','imag']):
			    if flagsave: end_figure(figQp[ir],axQp[ir],save=fileoutplotQp+'_'+r)
			    else: end_figure(figQp[ir],axQp[ir]);


		# finish impedance plot
		if flagsave: end_figure(figZ,axZ,save=fileoutZ)
		else: end_figure(figZ,axZ);


    if not(flagsave): pylab.show();
