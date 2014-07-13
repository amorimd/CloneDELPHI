#!/usr/bin/python

# library for computation of instabilities from the DELPHI code

import sys
import commands
# define user (for bsub command) (environment variable USER should be previously defined)
user=commands.getoutput("echo $USER");
if (len(user)>0): user_option=" -u "+user;
else: user_option='';
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);
import numpy as np
from string import *
import os,re
from io_lib import read_ncol_file,write_ncol_file
from string_lib import fortran_str,float_to_str
from tables_lib import *
from C_complex import *
from particle_param import *
from ctypes import *
from numpy import linalg as li
import time as ti  


libDELPHI=CDLL("libDELPHI.so");

prototypefact = CFUNCTYPE(c_long, c_int, c_int) # c_long is the type of the result
factpart = prototypefact(('factorialpart', libDELPHI));

prototypelag = CFUNCTYPE(c_double,c_int, c_int, c_double) # 1st c_double is the type of the result
laguerre = prototypelag(('Laguerre', libDELPHI));


def compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,Z,freqZ,flag_trapz=0,abseps=1e-4,lmaxold=-1,nmaxold=-1,couplold=None):

    ''' function to compute DELPHI impedance matrix '''
    
    aovertaub2=a/(taub*taub);
    bovertaub2=b/(taub*taub);
    Zsum=Complex(0.,0.);
    
    ng=len(g);nZ=len(freqZ);nZ2=2*nZ;
    glist=g.tolist();flist=freqZ.tolist();
    Z2=Z.reshape((1,-1))[0];Zlist=Z2.tolist();
    DArrayg=c_double * ng;
    DArrayf=c_double * nZ;
    DArrayZ=c_double * nZ2;
    
    g_arr=DArrayg(*glist);
    Z_arr=DArrayZ(*Zlist);
    freq_arr=DArrayf(*flist);

    # C function prototype
    prototypeisum = CFUNCTYPE(Complex, c_int, c_int, c_double, c_double, c_double, c_int, c_int, 
    	c_int, c_int, c_double, c_double, c_double, DArrayg, c_long, DArrayZ, DArrayf, c_long,
	c_int, c_double) # Complex is the type of the result
    isum = prototypeisum(('impedancesum', libDELPHI));
    
    coupl=np.zeros((2*lmax+1,nmax+1,2*lmax+1,nmax+1),dtype=complex);
    
    # compute matrix
    for l in range(lmax,-lmax-1,-1):
    
    	for lprime in range(lmax,-lmax-1,-1):
	
            coefj=1j**(lprime-l);
	    
    	    for n in range(nmax+1):
	    	
		fact=1./(factpart(n+abs(l),n)*2.**abs(l)); # =factorial(n)/(factorial(n+|l|)*2^|l|)
		
    	    	for nprime in range(nmax+1):
		
        	    if ( (l>=0)and(lprime>=0) ):

        		if (couplold==None) or ( ( (l>lmaxold)or(lprime>lmaxold) ) or ( (n>nmaxold)or(nprime>nmaxold) ) ):

                	    # matrix coefficient was never computed for this l, lprime, n, nprime
                	    Zsum=isum(nx, M, c_double(omegaksi), c_double(omega0),
			    	c_double(tunefrac), l, lprime, n, nprime, c_double(aovertaub2),
				c_double(bovertaub2), c_double(taub), g_arr, c_long(ng), Z_arr, freq_arr, c_long(nZ),
				flag_trapz, c_double(abseps));
			    #print Zsum.real,Zsum.imag

                	    coupl[l+lmax,n,lprime+lmax,nprime]=coefj*fact*(Zsum.real+1j*Zsum.imag);

        	        else:
		    	
                	    # matrix coefficient was already computed
                	    coupl[l+lmax,n,lprime+lmax,nprime]=couplold[l+lmaxold,n,lprime+lmaxold,nprime];

		    # other cases: uses the matrix symmetries
             	    elif ( (l<0)and(lprime>=0) ):

               		coupl[l+lmax,n,lprime+lmax,nprime]=coupl[-l+lmax,n,lprime+lmax,nprime];

             	    elif ( (l>=0)and(lprime<0) ):

               		coupl[l+lmax,n,lprime+lmax,nprime]=coupl[l+lmax,n,-lprime+lmax,nprime];

             	    elif ( (l<0)and(lprime<0) ):

               		coupl[l+lmax,n,lprime+lmax,nprime]=coupl[-l+lmax,n,-lprime+lmax,nprime];
			
    return coupl;
    

def compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,flagdamperimp=0,d=None,freqd=None,abseps=1e-4,lmaxold=-1,nmaxold=-1,damperold=None):

    ''' function to compute DELPHI damper matrix '''
    
    aovertaub2=a/(taub*taub);
    bovertaub2=b/(taub*taub);
    damp=Complex(0.,0.);
    
    ng=len(g);
    glist=g.tolist();
    DArrayg=c_double * ng;
    g_arr=DArrayg(*glist);

    if ((flagdamperimp==1)and(d!=None))and(freqd!=None):

	nd=len(freqd);nd2=2*nd;
	flist=freqd.tolist();
	d2=d.reshape((1,-1))[0];dlist=d2.tolist();

	DArrayf=c_double * nd;
	DArrayd=c_double * nd2;

	d_arr=DArrayd(*dlist);
	freq_arr=DArrayf(*flist);

	# C function prototype
	prototypedsum = CFUNCTYPE(Complex, c_int, c_int, c_double, c_double, c_double, c_int, c_int, 
    	    c_int, c_int, c_double, c_double, c_double, DArrayg, c_long, DArrayd, DArrayf, c_long,
	    c_double) # Complex is the type of the result
	dsum = prototypedsum(('dampersum', libDELPHI));

    else:
    
	# C function prototypes
	prototypeGln = CFUNCTYPE(c_double, c_int, c_int, c_double, c_double, c_double, c_long, DArrayg)
	Gln = prototypeGln(('Gln', libDELPHI))

	prototypeIln = CFUNCTYPE(c_double, c_int, c_int, c_double, c_double, c_double, c_double)
	Iln = prototypeIln(('Iln', libDELPHI))

    
    damper=np.zeros((2*lmax+1,nmax+1,2*lmax+1,nmax+1),dtype=complex);
    
    # compute matrix
    for l in range(lmax,-lmax-1,-1):
    
    	for lprime in range(lmax,-lmax-1,-1):
	
            coefj=1j**(lprime-l);
	    
    	    for n in range(nmax+1):
	    	
		fact=1./(factpart(n+abs(l),n)*2.**abs(l)); # =factorial(n)/(factorial(n+|l|)*2^|l|)
		
    	    	for nprime in range(nmax+1):
		
        	    if ( (l>=0)and(lprime>=0) ):

        		if (damperold==None) or ( ( (l>lmaxold)or(lprime>lmaxold) ) or ( (n>nmaxold)or(nprime>nmaxold) ) ):

                	    # matrix coefficient was never computed for this l, lprime, n, nprime

   		 	    if (flagdamperimp==1):
		   		# damper from an impedance-like function
             	    		damp=dsum(nx, M, c_double(omegaksi), c_double(omega0),
			    		c_double(tunefrac), l, lprime, n, nprime, c_double(aovertaub2),
					c_double(bovertaub2), c_double(taub), g_arr, c_long(ng), d_arr, freq_arr, c_long(nd),
					c_double(abseps));
			    	
				#print damp.real,damp.imag

                	    	damper[l+lmax,n,lprime+lmax,nprime]=coefj*fact*(damp.real+1j*damp.imag);


			    else:
				# bunch-by-bunch damper
				I_lprimenprime=Iln(lprime,nprime,c_double(-omegaksi),c_double(aovertaub2),c_double(bovertaub2),c_double(taub));
				G_ln=Gln(l,n,c_double(-omegaksi),c_double(aovertaub2),c_double(taub),c_long(ng),g_arr);
				damper[l+lmax,n,lprime+lmax,nprime]=coefj*fact*I_lprimenprime*G_ln;
		 
        	        else:
		    	
                	    # matrix coefficient was already computed
                	    damper[l+lmax,n,lprime+lmax,nprime]=damperold[l+lmaxold,n,lprime+lmaxold,nprime];

		    # other cases: uses the matrix symmetries
             	    elif ( (l<0)and(lprime>=0) ):

               		damper[l+lmax,n,lprime+lmax,nprime]=damper[-l+lmax,n,lprime+lmax,nprime];

             	    elif ( (l>=0)and(lprime<0) ):

               		damper[l+lmax,n,lprime+lmax,nprime]=damper[l+lmax,n,-lprime+lmax,nprime];

             	    elif ( (l<0)and(lprime<0) ):

               		damper[l+lmax,n,lprime+lmax,nprime]=damper[-l+lmax,n,-lprime+lmax,nprime];
			
    return damper;
    

    
def damper_imp_Karliner_Popov(L0,L1,L2,tauf,R,Q,f0,f):

    ''' evaluate (to a multiplicative factor) "impedance" of a damper
    (Karliner-Popov, Nucl. Inst. Meth. Phys. Res. A 2005)

    L0: distance from pickup to kicker (m),
    L1: pickup-length (m),
    L2: kicker length (m),
    tauf: 1/tauf ~ frequency of low-pass filter,
    R: machine radius,
    Q: tune,
    f0: revolution frequency,
    f: array of frequencies at which impedance is evaluated. '''

    c=299792458;

    omega0=2.*np.pi*f0;
    omega=2.*np.pi*f;

    tau=(L0-2.*L2)/c;

    # tmp=gamma_m-j*m/R
    tmp=1j*(-2.*omega/omega0-Q)/R;

    # filter
    K=1./(1.-1j*omega*tauf);

    Zimp=K*np.exp(-1j*omega0*Q*tau)*(1.-np.exp(-tmp*L1))*(1.-np.exp(-tmp*L2))/tmp;
    
    return Zimp;
    
def computes_coef(f0,dmax,b,g0,dnormfactor,taub,dphase,M,Nb,gamma,Q,particle='proton'):

    ''' compute coefficients in front of damper and imepdance matrices in the final eigenvalue system

    - f0: revolution frequency,
    - dmax: damper gain (inverse of number of damping turns) - depends also on normalization factor dnormfactor,
    - b: b parameter in DELPHI (for Laguerre poly. decomposition),
    - g0: first term in distribution decomposition over Laguerre polynomial,
    - dnormfactor: normalization factor for damper matrix. Usually such that damping rate of mode nx=0, l=0
        and n=0 (i.e. all coupled-bunch, azimuthal and radial mode numbers =0) is exactly dmax. It is precomputed
        beforehand (either at the current chromaticity, or at Q'=0),
    - taub: bunch length (total, or 4*RMS for Gaussian bunches) in seconds,
    - dphase: additional phase of the damper in radian (w.r.t. to purely resistive damper when dphase=0) 
    - M: number of bunches,
    - Nb: number of particles per bunch,
    - gamma: relativistic velocity factor,
    - Q: tune (with integer part included),
    - particle: 'proton' or 'electron'.
    '''
    
    # particle mass in kg
    if particle.startswith('proton'): m0=1.6726e-27;
    elif particle.startswith('electron'): m0=9.1094e-31;
    else: print "Pb with particle type";sys.exit();
    
    e=1.60218e-19; # elementary charge in C
    clight=299792458; # speed of light in m/s
    Ib=M*Nb*f0*e; # beam current
    
    bovertaub2=b/(taub*taub);
    
    coefdamper=(1j*f0*dmax*2.*bovertaub2/(g0*dnormfactor))*np.exp(1j*dphase);
    coefZ=1j*Ib*e/(2.*gamma*m0*clight*Q);
    
    return coefdamper,coefZ;
    

def eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=False):

    ''' compute and diagonalize the final matrix (including impedance, damper and Qs)

    - lmax: number of azimuthal modes in matrices (we go from -lmax to +lmax)
    - nmax: number of radial modes in matrices, minus one (we go from 0 to +nmax)
    - matdamper: precomputed damper matrix,
    - matZ: precomputed imepdance matrix,
    - coefdamper: coefficient in front of damper matrix
    - coefZ: coefficient in front of impedance matrix
    - omegas: synchrotron angular frequency
    - flageigenvect: True to compute eigenvect (faster without)
    
    '''
    
    Mat=np.zeros((2*lmax+1,nmax+1,2*lmax+1,nmax+1),dtype=complex);
    
    for l in range(lmax,-lmax-1,-1):
    
    	for lprime in range(lmax,-lmax-1,-1):
	
    	    for n in range(nmax+1):
	    	
    	    	for nprime in range(nmax+1):

		    Mat[l+lmax,n,lprime+lmax,nprime] = coefdamper*matdamper[l+lmax,n,lprime+lmax,nprime];
		    Mat[l+lmax,n,lprime+lmax,nprime]+= coefZ*matZ[l+lmax,n,lprime+lmax,nprime];
		    
		    if ( (l==lprime)and(n==nprime) ): Mat[l+lmax,n,lprime+lmax,nprime] += omegas*l;
    
  
    #t1=ti.clock()
    
    if (flageigenvect):
        # compute eigenvectors
    	eigenval,eigenvect=li.eig(Mat.reshape(((2*lmax+1)*(nmax+1),(2*lmax+1)*(nmax+1))));
    else:
    	# do not compute eigenvectors (~ 2 times faster)
	eigenval=li.eigvals(Mat.reshape(((2*lmax+1)*(nmax+1),(2*lmax+1)*(nmax+1))));
    	eigenvect=np.zeros(((2*lmax+1)*(nmax+1),(2*lmax+1)*(nmax+1)));
    
    #t2=ti.clock();
    #print "lmax=",lmax,", nmax=",nmax,", time for eigenvalue pb resolution [seconds]: ",t2-t1;
    
    return eigenval,eigenvect;


def eigenmodesDELPHI_converged(nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,Z,freqZ,coefdamper,coefZ,
	omegas,flag_trapz=0,flagdamperimp=0,d=None,freqd=None,kmax=5,crit=5.e-2,abseps=1.e-3,
	lmaxold=-1,nmaxold=-1,matdamperold=None,matZold=None,flageigenvect=False):

    ''' computes eigenmodes for increasing matrix size, until we get convergence on 
    the imaginary part of eigenvalues
    
    - nx: coupled-bunch mode considered (from 0 to M-1),
    - M: number of equidistant bunches,
    - omegaksi: chromatic angular frequency,
    - omega0: revolution angular frequency,
    - tunefrac: fractional part of the tune,
    - a,b: parameters for Laguerre polynomial decomposition in DELPHI,
    - taub: total bunch length (seconds) or 4*RMS for Gaussian bunches,
    - g: array with coef. of decomposition of initial distribution on Laguerre polynomials,
    - Z, freqZ: Z=damper impedance at the frequencies (NOT ANGULAR) freqZ
    - coefdamper, coefZ: coefficients in fromt resp. of damper and impedance matrix (precomputed beforehand),
    - omegas: synchrotron angular frequency.
    - flag_trapz: flag to use trapz method for computation of impedance matrix (if 1),
    - flagdamperimp: flag to use frequency dependent damper gain (if 1)(with d and freqd arrays),
    - d, freqd: d=damper impedance at the frequencies (NOT ANGULAR) freqd - used only if freqdamperimp==1,
    - kmax: number of eigenvalues to make converge,
    - crit: relative error tolerated on imaginary part of eigenvalues,
    - abseps: ~absolute error tolerated on imaginary part of eigenvalues (also used for absolute error determination
          on impedance and damper sums),
    - lmaxold, nmaxold: max azimuthal mode number of max radial mode number (radial mode begins at 0)
          already computed (in matdamperold and matZold),
    - matdamperold ,matZold: damper and impedance matrices already computed.
    '''

    eigenvalold=np.ones(kmax)*(1.e50+1j*1.e50);
    err=1e50;lmax=1;nmax=1;
    
    # absolute errors tolerated in damper and impedance sums computations 
    absepsd=abseps/max(abs(coefdamper),10);#print "absepsd=",absepsd
    absepsZ=abseps/max(abs(coefZ.real),10);#print "absepsZ=",absepsZ

    flagcrit=True;
    
    while flagcrit:
    
        #t1=ti.clock()
    	
	matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,
		flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=absepsd,lmaxold=lmaxold,
		nmaxold=nmaxold,damperold=matdamperold);

    	matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,Z,freqZ,
		flag_trapz=flag_trapz,abseps=absepsZ,lmaxold=lmaxold,nmaxold=nmaxold,couplold=matZold);

	#t2=ti.clock();
	#print "lmax=",lmax,", nmax=",nmax,", time for computing matrices [seconds]: ",t2-t1;
	
	if ((len(np.where(np.isinf(matdamper))[0])==0)and(len(np.where(np.isnan(matdamper))[0])==0))and((len(np.where(np.isinf(matZ))[0])==0)and(len(np.where(np.isnan(matZ))[0])==0)):

	    if (lmax<16)and(nmax<16):
	    	eigenval,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=flageigenvect);
	    else:
	    	print "	warning: eigenmodesDELPHI_converged: too big matrices, convergence stopped, error=",err,", lmax=",lmax,", nmax=",nmax;
		
	else:
	    print "	warning: eigenmodesDELPHI_converged: some nan or inf in matrix, convergence stopped, error=",err,", lmax=",lmax,", nmax=",nmax;
	    break;
	
	ind=np.argsort(np.imag(eigenval));
	
	kmaxtmp=min(kmax,len(eigenval));
	
	err=np.max(np.abs(np.imag(eigenval[ind[:kmaxtmp]])
		-np.imag(eigenvalold[:kmaxtmp]))
		/(np.abs(np.imag(eigenvalold[:kmaxtmp]))+abseps));
	
	flagcrit=(err>crit)or(kmaxtmp<kmax);
		
	for k in range(min(kmax,len(eigenval))): eigenvalold[k]=eigenval[ind[k]];
	
	if (lmax>lmaxold)and(nmax>nmaxold):
	    lmaxold=lmax;nmaxold=nmax;
	    matZold=matZ;matdamperold=matdamper;

	lmax+=1;nmax+=1;


    return eigenval[ind],v[:,ind],lmaxold,nmaxold,matdamperold,matZold;
	

def eigenmodesDELPHI_converged_scan(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,
	omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle='proton',flagnorm=0,
	flag_trapz=0,flagdamperimp=0,d=None,freqd=None,
	kmax=1,kmaxplot=10,crit=5.e-2,abseps=1.e-3,flagm0=False):
	
    ''' encapsulate DELPHI calculations, with scans on coupled-bunch modes, damper gain,
    nb of particles, synchrotron tune and damper phase.
    return tuneshifts (complex) for all these parameters scanned
    and also the tuneshifts of the most unstable coupled-bunch modes '''

    f0=omega0/(2.*np.pi);
    Qfrac=Q-np.floor(Q);

    if flagnorm==0:
	# normalization factor for damper
	dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,Qfrac,a,b,taub,g,
    	    flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=abseps);
	dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
    tuneshift_most=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
    if flagm0:
        tuneshiftm0nx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan)),dtype=complex);
	tuneshiftm0=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan)),dtype=complex);

    for iQp,Qp in enumerate(Qpscan):
    
   	omegaksi=Qp*omega0/eta;

	for inx,nx in enumerate(nxscan):

	    if flagnorm:
		# normalization factor for damper
		dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,
	    	    flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=abseps);
		dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

	    lmaxi=-1;nmaxi=-1; # to keep track of maximum lmax and nmax for all modes
   	    lmax=-1;nmax=-1;matZ=None;matdamper=None;

	    for idamp,damp in enumerate(dampscan):

		for iNb,Nb in enumerate(Nbscan):

		    for iomegas,omegas in enumerate(omegasscan):

			for idphase,dphase in enumerate(dphasescan):

			    #print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;

			    coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,
				    dphase,M,Nb,gamma,Q,particle=particle);

			    freqshift,v,lmax,nmax,matdamper,matZ=eigenmodesDELPHI_converged(nx,
				    M,omegaksi,omega0,Qfrac,a,b,taub,g,Z,freq,coefdamper,
				    coefZ,omegas,flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,
				    freqd=freqd,kmax=kmax,crit=crit,abseps=abseps,lmaxold=lmax,
				    nmaxold=nmax,matdamperold=matdamper,matZold=matZ);
			    
			    if flagm0:
				# extract mode m=0 tuneshift
				 tuneshiftm0nx[iQp,inx,idamp,iNb,iomegas,idphase]=extract_between_bounds(freqshift/omega0,-omegas/(2.*omega0),omegas/(2.*omega0));

			    if (len(freqshift)>=kmaxplot):
			        tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,:]=freqshift[:kmaxplot]/omega0;
			    else:
			        tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,:len(freqshift)]=freqshift[:]/omega0;
			        
			    #lambdax[iNb,:]=freqshift[:kmaxplot]/omegas;
			    lmaxi=max(lmax,lmaxi);nmaxi=max(nmax,nmaxi);	

	# find the most unstable coupled-bunch mode
	for idamp,damp in enumerate(dampscan):

	    for iNb,Nb in enumerate(Nbscan):

		for iomegas,omegas in enumerate(omegasscan):

		    for idphase,dphase in enumerate(dphasescan):

			for kmode in range(kmaxplot):

			    inx=np.argmin(np.imag(tuneshiftnx[iQp,:,idamp,iNb,iomegas,idphase,kmode]))
			    if kmode<kmax:
				print "Qp=",Qp,", M=",M,", d=",damp,", Nb=",Nb,", omegas=",omegas,", dphase=",dphase;
				print "   lmaxi=",lmaxi,", nmaxi=",nmaxi,", kmode=",kmode,", Most unstable coupled-bunch mode: ",nxscan[inx];
			    tuneshift_most[iQp,idamp,iNb,iomegas,idphase,kmode]=tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,kmode];

			if flagm0:

			    inx=np.argmax(np.abs(np.real(tuneshiftm0nx[iQp,:,idamp,iNb,iomegas,idphase])))
			    if kmode<kmax:
				print "Qp=",Qp,", M=",M,", d=",damp,", Nb=",Nb,", omegas=",omegas,", dphase=",dphase;
				print "   lmaxi=",lmaxi,", nmaxi=",nmaxi,", kmode=",kmode,", Most unstable coupled-bunch mode: ",nxscan[inx];
			    tuneshiftm0[iQp,idamp,iNb,iomegas,idphase]=tuneshiftm0nx[iQp,inx,idamp,iNb,iomegas,idphase];

    if flagm0: return tuneshift_most,tuneshiftnx,tuneshiftm0;
    else: return tuneshift_most,tuneshiftnx;
    

def eigenmodesDELPHI_converged_scan_lxplus(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,
	omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle='proton',flagnorm=0,
	flag_trapz=0,flagdamperimp=0,d=None,freqd=None,
	kmax=1,kmaxplot=10,crit=5.e-2,abseps=1.e-3,flagm0=False,
	lxplusbatch=None,comment='',queue='1nh',dire=''):
	
    ''' same as eigenmodesDELPHI_converged_scan with possibility to launch on lxplus
    lxplusbatch: if None, no use of any queuing system
                 if 'launch' -> launch calculation on lxplus (or any LSF batch system) on queue 'queue'
                 if 'retrieve' -> retrieve outputs
    comment is used to identify the batch job name and the pickle filename
    dire is the directory where to put/find the result; it should be a path 
    relative to the current directory (e.g. '../') '''
    
    import pickle as pick;
    import commands;
    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
    tuneshift_most=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
    if flagm0: tuneshiftm0=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan)),dtype=complex);
    
    if (lxplusbatch==None):
        if flagm0: tuneshift_most,tuneshiftnx,tuneshiftm0=eigenmodesDELPHI_converged_scan(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,
		omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=freqd,
		kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0);
	
	else: tuneshift_most,tuneshiftnx=eigenmodesDELPHI_converged_scan(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,
		omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=freqd,
		kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0);

    else:
    
	filename="DELPHI_pickle_"+comment+".dat";
	status,lsres=commands.getstatusoutput('ls '+dire+'out'+filename);
	#print status,lsres
        
	if lxplusbatch.startswith('launch')and(status!=0): # test if output file not already there
    
    	    # pickle file
    	    fileall=open(filename,'w');
    	    listnames=['Qpscan','nxscan','dampscan','Nbscan','omegasscan','dphasescan','M',
		    'omega0','Q','gamma','eta','a','b','taub','g','Z','freq',
		    'particle','flagnorm','flag_trapz','flagdamperimp','d','freqd',
		    'kmax','kmaxplot','crit','abseps','flagm0']
    	    for name in listnames: eval('pick.dump('+name+',fileall)');
	    fileall.close();

            # launch calculation on lxplus batch system
            filejob=open('batch'+comment+'.job','w');
            here=os.getcwd()+'/';#print here;
            print >> filejob, "mv "+here+"/"+filename+" .";
            #print >> filejob, "cp "+here+"/DELPHI.py .";
            #print >> filejob, "cp "+here+"/DELPHI_script.py .";
	    #print >> filejob, "cp "+here+"/plot_lib.py .";
	    #print >> filejob, "cp "+here+"/io_lib.py .";
	    #print >> filejob, "cp "+here+"/string_lib.py .";
	    #print >> filejob, "cp "+here+"/tables_lib.py .";
	    #print >> filejob, "cp "+here+"/C_complex.py .";
 	    #print >> filejob, "cp "+here+"/particle_param.py .";
            #print >> filejob, "chmod +x DELPHI_script.py";
            print >> filejob, "DELPHI_script.py "+filename+" > out_"+comment;
            print >> filejob, "cp out"+filename+" "+here+dire;
            print >> filejob, "cp out_"+comment+" "+here+dire;
            filejob.close();
            os.system("chmod 744 batch"+comment+".job");
            os.system("bsub"+user_option+" -e error1.out -q "+queue+" batch"+comment+".job");

    	else: # 'retrieve' case
    	    
	    status,lsres=commands.getstatusoutput('ls '+dire+'out'+filename);
	    if status==0:
		# pickle output file
    		fileall=open(dire+'out'+filename,'r');
    		tuneshift_most=pick.load(fileall);
		tuneshiftnx=pick.load(fileall);
		if flagm0: tuneshiftm0=pick.load(fileall);
		fileall.close();
	    
	    else:
	    	# fill tables with nans
		tuneshift_most.fill(np.nan+1j*np.nan);
		tuneshiftnx.fill(np.nan+1j*np.nan);
		if flagm0: tuneshiftm0.fill(np.nan+1j*np.nan);
		print " WARNING !!! no file",filename;
	    
	    os.system("rm -rf LSFJOB_* error1.out batch"+comment+".job");
	
    if flagm0: return tuneshift_most,tuneshiftnx,tuneshiftm0;
    else: return tuneshift_most,tuneshiftnx;
		
	
def DELPHI_wrapper(imp_mod,Mscan,Qpscan,dampscan,Nbscan,omegasscan,dphasescan,
	omega0,Qx,Qy,gamma,eta,a,b,taub,g,planes=['x','y'],nevery=1,particle='proton',flagnorm=0,
	flagdamperimp=0,d=None,freqd=None,
	kmax=1,kmaxplot=10,crit=5.e-2,abseps=1.e-3,flagm0=False,
	lxplusbatch=None,comment='',queue='1nh',dire='',flagQpscan_outside=True):
	
    ''' wrapper of eigenmodesDELPHI_converged_scan_lxplus, with more scans
    imp_mod is an impedance or wake model (see Impedance.py)
    from which one extract the planes given in 'planes'
    Mscan is the scan in number of bunches
    lxplusbatch: if None, no use of lxplus batch system
                 if 'launch' -> launch calculation on lxplus on queue 'queue'
                 if 'retrieve' -> retrieve outputs
    comment is used to identify the batch job name and the pickle filename
    dire is the directory where to put/find the result; it should be a path 
    relative to the current directory (e.g. '../')
    nevery indicates the downsampling (we take 1 frequency every "nevery" 
    frequencies), for the impedance.
    flagQpscan_outside: True to put the Qpscan outside the lxplus batch job,
    False so that it is inside the lxplus batch job.
    '''
    
    from Impedance import test_impedance_wake_comp,impedance_wake
    from copy import deepcopy

    strnorm=['','_norm_current_chroma'];

    tuneshiftQp=np.zeros((len(planes),len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
    if flagm0: tuneshiftm0Qp=np.zeros((len(planes),len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan)),dtype=complex);

    for iplane,plane in enumerate(planes):
    	
	# select Zxdip or Zydip
	flagx=int(plane=='x');
	for iw in imp_mod:
	    if test_impedance_wake_comp(iw,flagx,1-flagx,0,0,plane): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);

	for iM,M in enumerate(Mscan):

	    flag_trapz=0; # by default no trapz method

	    if (M==1): nxscan=np.array([0]);flag_trapz=1;
	    #elif (M==1782): nxscan=np.array([0, 1, 300, 600, 880, 890, 891, 892, 900, 910, 950, 1000, 1200, 1500, 1780, 1781])
	    #elif (M==3564): nxscan=np.array([0, 1, 300, 600, 900, 1200, 1500, 1770, 1780, 1781, 1782, 1785, 1790, 1800, 1900, 2000, 2300, 2600, 2900, 3200, 3500, 3560, 3561, 3562, 3563])
	    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,M/20),np.arange(M/2-10,M/2+11),
		np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

	    if flagQpscan_outside:
		# put loop on Qpscan outside the lxplus job
		for iQp,Qp in enumerate(Qpscan):
		    tuneshiftnx=np.zeros((1,len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
		    if flagm0:
			tuneshiftQp[iplane,iM,iQp,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[iplane,iM,iQp,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			    kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			    lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_Qp'+str(Qp)+'_'+plane,
			    queue=queue,dire=dire);

		    else:
			tuneshiftQp[iplane,iM,iQp,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			    kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			    lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
			    queue=queue,dire=dire);
	    
	    else:
		# put loop on Qpscan inside the lxplus job
		tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),kmaxplot),dtype=complex);
		if flagm0:
		    tuneshiftQp[iplane,iM,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[iplane,iM,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,omega0,eval('Q'+plane),gamma,eta,
			a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_'+plane,
			queue=queue,dire=dire);

		else:
		    tuneshiftQp[iplane,iM,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			nxscan,dampscan,Nbscan,omegasscan,dphasescan,M,omega0,eval('Q'+plane),gamma,eta,
			a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_'+strnorm[flagnorm]+'_'+plane,
			queue=queue,dire=dire);
	    

    if flagm0: return tuneshiftQp,tuneshiftm0Qp;
    else: return tuneshiftQp;
    

def dispersion_integral_oct_2D(Q,bx,bxy,distribution='gaussian'):

    ''' computes the dispersion integral in 2D from an octupole with detuning 
    per sigma given by bx (in the plane of the coherent motion) and bxy (in the
    other plane). Compute the integral for a certain coherent complex tune shift Q.
    The transverse distribution can be 'gaussian' or 'parabolic'.
    This is the integral vs Jx and Jy of Jx*dphi/dJx/(Q-bx*Jx-bxy*Jy-i0) (with phi the distribution function)
    
    NOTE: for stability diagrams, use -1/dispersion_integral, and usually the convention is to plot
    -Im[Q] vs Re[Q] (beware of all these minus signs).
    '''
    
    from scipy.special import exp1;
    
    if (np.imag(Q)==0): Q=Q-1e-15*1j;
    c=bxy/bx;
    q=-Q/bx;
    
    
    if (distribution=='gaussian'):
    
    	I1=(1-c-(q+c-c*q)*np.exp(q)* exp1(q) + c*np.exp(q/c)* exp1(q/c)) / ((1-c)**2);
	
	if np.isnan(I1): I1=1./q-(c+2)/q**2; # asymptotic form for large q (assuming c is of order 1)

    elif (distribution=='parabolic'):
    	
	xi=q/5.;
	
        #I1=((c+xi)**3 *np.log(1+xi)-(c+xi)**3.*np.log(c+xi)+(-1+c)*(c*(c+2*c*xi+(-1+2*c)*xi**2) + 
	#	(-1+c)*xi**2 *(3*c+xi+2*c*xi)*(np.log(xi)-np.log(1+xi))))/((-1+c)**2*c**2);
        I1=((c+xi)**3 *np.log((1+xi)/(c+xi))+(-1+c)*(c*(c+2*c*xi+(-1+2*c)*xi**2) + 
		(-1+c)*xi**2 *(3*c+xi+2*c*xi)*np.log(xi/(1+xi))))/((-1+c)**2*c**2);
	I1=-I1*4./5.;
	# I checked this is the same as in Scott Berg-Ruggiero CERN SL-AP-96-71 (AP)
	
	if (np.abs(xi)>100.): I1=1./q-(c+2)/q**2; # asymptotic form for large q (assuming c is of order 1) (obtained thanks to Mathematica - actually the same as for Gaussian)

    I=-I1/bx;
    I=-I; # additional minus sign because for DELPHI we want the integral with 
    # dphi/dJx (derivative of distribution) on the numerator, so -[the one of Berg-Ruggiero]
    
    return I;


def detuning_coef_oct_LHC(plane,gamma,epsnormx,epsnormy,current_foc_oct,current_defoc_oct):

    ''' compute detuning coefficients (per transverse sigma) bx (in the plane of the coherent motion)
    and bxy (in the other plane) for LHC octupoles with standard (i.e. non ATS) optics 
    (squeeze does not matter, but ATS does).
    Input parameters:
    - plane: plane of coherent motion: 'x' if horizontal, 'y' if vertical
     (sligthly different coefficients).
    - gamma: relativistic mass factor of the beam,
    - epsnormx: normalized emittance in the plane (x or y) studied (e.g. 3.75e-6 m at 7TeV),
    - epsnormy: normalized emittance in the plane (x or y) perpendicular to the plane studied (e.g. 3.75e-6 m at 7TeV),
    - current_foc_oct: current in the focusing octupoles (max is supposed to be 550A),
    - current_defoc_oct: current in the defocusing octupoles (max is supposed to be 550A),
    '''
    
    current_max=550.;

    beta=np.sqrt(1.-1./gamma**2); # relativistic velocity factor
    # reduction factor for the octupole current and the energy
    F=(current_foc_oct/current_max)*(7460.52/gamma);
    D=(current_defoc_oct/current_max)*(7460.52/gamma);

    eps1sigmax = epsnormx/(beta*gamma);
    eps1sigmay = epsnormy/(beta*gamma);

    # use madx value of O_3 (63100 T/m3) and madx beta functions (checked with S. Fartoukh).
    # We put the same values as in HEADTAIL.
    if plane=='y':
	# vertical
	ax=9789*F-277203*D;
    elif plane=='x':
	# horizontal
	ax=267065*F-7856*D;

    axy=-102261*F+93331*D;

    bx=ax*eps1sigmax;
    bxy=axy*eps1sigmay;
    
    return bx,bxy;
    
    	
def determinantDELPHI_tunespread(dQc,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'):

    ''' compute the determinant of the final matrix (including impedance, damper and transverse Landau damping) 
    as a function of coherent tune shift dQc. This function has to be solved (vs. dQc) to find 
    the eigenmodes.

    - dQc: coherent tune shift (from unperturbed tune Q0).
    - lmax: number of azimuthal modes in matrices (we go from -lmax to +lmax)
    - nmax: number of radial modes in matrices, minus one (we go from 0 to +nmax)
    - matdamper: precomputed damper matrix,
    - matZ: precomputed impedance matrix,
    - coefdamper: coefficient in front of damper matrix
    - coefZ: coefficient in front of impedance matrix
    - bx & bxy: detuning coefficients for the transverse tunespread (in units of sigma),
      in the plane of the coherent motion (bx) and in the other plane (bxy)
    - omega0: revolution angular frequency
    - omegas: synchrotron angular frequency
    - distribution: kind of transverse distribution ('gaussian' or 'parabolic')
    '''
    
    Mat=np.zeros((2*lmax+1,nmax+1,2*lmax+1,nmax+1),dtype=complex);

    for l in range(lmax,-lmax-1,-1):

    	for lprime in range(lmax,-lmax-1,-1):

    	    for n in range(nmax+1):

    	    	for nprime in range(nmax+1):

		    Mat[l+lmax,n,lprime+lmax,nprime] = coefdamper*matdamper[l+lmax,n,lprime+lmax,nprime];
		    Mat[l+lmax,n,lprime+lmax,nprime]+= coefZ*matZ[l+lmax,n,lprime+lmax,nprime];

		    if ( (l==lprime)and(n==nprime) ):

		    	# tune offset from azimuthal mode number (synchrotron sideband)
			Qls=omegas*l/omega0;
			# dispersion integral (for transverse Landau damping)
			if np.abs(bx*bxy)>0.:
			    I=dispersion_integral_oct_2D(dQc-Qls,bx,bxy,distribution=distribution);
		    	    # final matrix with Landau damping included
			    Mat[l+lmax,n,lprime+lmax,nprime] += omega0/I;
			else:
		    	    # final matrix without Landau damping
			    Mat[l+lmax,n,lprime+lmax,nprime] += omega0*(Qls-dQc);

    det=li.det(Mat.reshape(((2*lmax+1)*(nmax+1),(2*lmax+1)*(nmax+1))));
    #print "dQc=",dQc,", det=",det;
  
    return det;


def determinantDELPHI_tunespread_array(dQc,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'):

    ''' compute the determinant of the final matrix (including impedance, damper and transverse Landau damping) 
    as a function of coherent tune shift dQc. Version where dQc can be an array.

    - dQc: coherent tune shift (from unperturbed tune Q0). Can be a scalar or an array.
    - lmax: number of azimuthal modes in matrices (we go from -lmax to +lmax)
    - nmax: number of radial modes in matrices, minus one (we go from 0 to +nmax)
    - matdamper: precomputed damper matrix,
    - matZ: precomputed impedance matrix,
    - coefdamper: coefficient in front of damper matrix
    - coefZ: coefficient in front of impedance matrix
    - bx & bxy: detuning coefficients for the transverse tunespread (in units of sigma),
      in the plane of the coherent motion (bx) and in the other plane (bxy)
    - omega0: revolution angular frequency
    - omegas: synchrotron angular frequency
    - distribution: kind of transverse distribution ('gaussian' or 'parabolic')
    '''
    
    dQc_list=create_list(dQc,n=1);
    det=np.zeros(len(dQc_list),dtype=complex);
    
    for idQc,dQc in enumerate(dQc_list):
    
    	det[idQc]=determinantDELPHI_tunespread(dQc,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian')

    return det;
    

def solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,
	omega0,omegas,distribution='gaussian',kini=12):#,minimag=None,maximag=None,npts_all=None,npts_zoom=100,kini=10):

    ''' solve the equation determinant(final matrix)=0 (matrix includes impedance, damper and transverse Landau damping) 
    as a function of coherent tune shift, to find the eigenmodes.

    - lmax: number of azimuthal modes in matrices (we go from -lmax to +lmax)
    - nmax: number of radial modes in matrices, minus one (we go from 0 to +nmax)
    - matdamper: precomputed damper matrix,
    - matZ: precomputed impedance matrix,
    - coefdamper: coefficient in front of damper matrix,
    - coefZ: coefficient in front of impedance matrix,
    - bx & bxy: detuning coefficients for the transverse tunespread (in units of sigma),
      in the plane of the coherent motion (bx) and in the other plane (bxy),
    - omega0: revolution angular frequency,
    - omegas: synchrotron angular frequency,
    - distribution: kind of transverse distribution ('gaussian' or 'parabolic'),
    - kini: number of eigenvalues of the linear system (without Landau damping)
    to initialize the solving algorithm.
    '''
    
    from solve_lib import potential_roots_complex,roots_complex_list;

    #if minimag==None: minimag=-omegas/omega0;
    #if maximag==None: maximag=0;
    #if npts_all==None: npts_all=2*lmax+1;

    # f=function to solve is the matrix determinant (normalized to the value for zero coherent tuneshift)
    factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
    f=(lambda dQc: determinantDELPHI_tunespread(dQc,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution)/factnorm);
    
    # first solve the eigenvalue problem without Landau damping (provide initial guesses)
    freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=False);
    ind=np.argsort(np.imag(freqshift));
    # kini cannot be more than the number of eigenvalues with negative imag. part
    kini=max(kini,len(np.where(np.imag(freqshift[ind])<0)));

    list_pot_roots=[];

#    # OLD pre-solving algorithm (find all potential roots) (perform very poorly)
#    list_pot_roots=potential_roots_complex(f,-lmax*(1.+1./(npts_all-1.))*omegas/omega0,lmax*(1.+1./(npts_all-1.))*omegas/omega0,minimag,maximag,
#    	npts=npts_all,meshkind='lin');
#    print "pre-solving (global) -> done",len(list_pot_roots);
#
#    # 1st version: zoom around synchrotron sidebands
#    for l in range(-lmax,lmax+1):
#    	
#	#list_tmp=potential_roots_complex(f,(l-0.5)*omegas/omega0,(l+0.5)*omegas/omega0,minimag,maximag,
#    	#	npts=npts_zoom,meshkind='log',offset_log=l*omegas/omega0+1j*(minimag+maximag)/2.);
#	list_tmp=potential_roots_complex(f,(l-0.5)*omegas/omega0,(l+0.5)*omegas/omega0,minimag,maximag,
#    		npts=npts_zoom,meshkind='lin',offset_log=l*omegas/omega0);
#	
#    	print "pre-solving (l=",l,") -> done",len(list_tmp),(l-0.5)*omegas/omega0,(l+0.5)*omegas/omega0,minimag,maximag;
#	
#	for el in list_tmp:
#	    if not(el in list_pot_roots): list_pot_roots.append(el);
#   
#    # 2nd version: zoom around eigenvalues of unperturbed system
#    for i in ind[:kini]:
#    	
#	minimag=np.imag(freqshift[i])/omega0;
#	minreal=min(np.real(freqshift[i])/omega0,np.real(freqshift[i])/omega0+2*(2*bx+bxy));
#	maxreal=max(np.real(freqshift[i])/omega0,np.real(freqshift[i])/omega0+2*(2*bx+bxy));
#	list_tmp=potential_roots_complex(f,minreal,maxreal,minimag,maximag,
#    		npts=npts_zoom,meshkind='lin',offset_log=np.real(freqshift[i])/omega0);
#	
#    	print "pre-solving (i=",i,") -> done",len(list_tmp),2*(2*bx+bxy),minreal,maxreal,minimag,maximag;
#	
#	for el in list_tmp:
#	    if not(el in list_pot_roots): list_pot_roots.append(el);

    # NEW algorithm, much simpler (and works): take unperturbed eigenvalues
    # as initial guesses, plus the same with a real offset from the tunespread
    for fr in freqshift[ind[:kini]]:
	list_pot_roots.append((np.real(fr)/omega0,np.imag(fr)/omega0));
    	if (2*abs(2*bx+bxy)>1e-7): list_pot_roots.append((np.real(fr)/omega0+2*(2*bx+bxy),np.imag(fr)/omega0));
	
	# then zoom around eigenvalues of unperturbed system
	minimag=np.imag(fr)/omega0;maximag=0.;
	minreal=min(np.real(fr)/omega0,np.real(fr)/omega0+2*(2*bx+bxy));
	maxreal=max(np.real(fr)/omega0,np.real(fr)/omega0+2*(2*bx+bxy));
	list_tmp=potential_roots_complex(f,minreal,maxreal,minimag,maximag,
    		npts=20,meshkind='lin',offset_log=np.real(fr)/omega0);
	
    	print "pre-solving (unperturbed dQ=",fr/omega0,") -> done",len(list_tmp),list_tmp;
	
	for el in list_tmp:
	    if not(el in list_pot_roots): list_pot_roots.append(el);


    print list_pot_roots;
    #arr_pot_roots=np.array(list_pot_roots)*omega0;
    
    # actual solving algorithm (note: roots translated from tune shifts to angular frequency shits)
    eigenval=omega0*roots_complex_list(f,list_pot_roots,tolf=1e-3,tolx=1e-10);
    
    return eigenval;#,arr_pot_roots;


def solve_stability_diagram(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,
	omega0,omegas,distribution='gaussian',kini=12):#,minimag=None,maximag=None,npts_all=None,npts_zoom=100,kini=10):

    ''' solve the stability diagram equation dQu*dispersion_integral(Qc)=-1 
    with dQu the unperturbed coherent tune shifts from DELPHI (using impedance & damper
    but not Landau damping) as a function of coherent tune shift Qc, to 
    find the eigenmodes with transverse Landau damping (stability diagram approx.).

    - lmax: number of azimuthal modes in matrices (we go from -lmax to +lmax)
    - nmax: number of radial modes in matrices, minus one (we go from 0 to +nmax)
    - matdamper: precomputed damper matrix,
    - matZ: precomputed impedance matrix,
    - coefdamper: coefficient in front of damper matrix,
    - coefZ: coefficient in front of impedance matrix,
    - bx & bxy: detuning coefficients for the transverse tunespread (in units of sigma),
      in the plane of the coherent motion (bx) and in the other plane (bxy),
    - omega0: revolution angular frequency,
    - omegas: synchrotron angular frequency,
    - distribution: kind of transverse distribution ('gaussian' or 'parabolic'),
    - kini: number of eigenvalues of the linear system (without Landau damping)
    to initialize the solving algorithm.
    '''
    
    from solve_lib import potential_roots_complex,roots_complex_list;

    # first solve the eigenvalue problem without Landau damping (provide initial 
    # guesses & unperturbed coherent tune shifts)
    freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=False);
    ind=np.argsort(np.imag(freqshift));
    # kini cannot be more than the number of eigenvalues with negative imag. part
    kini=max(kini,len(np.where(np.imag(freqshift[ind])<0)));

    eigenval=[];
    
    # solve for each unperturbed eigenvalue
    for fr in freqshift[ind[:kini]]:

	# evaluate azimuthal mode number (used for stability diagram)
	l=np.int(np.ceil(fr.real/omegas));
	#print "fr=",fr/omega0,", l=",l;
	
	# solve equation
	Qls=omegas*l/omega0;
	# f: function to solve is the stability diagram eq. : dQu*dispersion_integral(Qc)=-1
	f=(lambda dQc: (fr/omega0)*dispersion_integral_oct_2D(dQc-Qls,bx,bxy,distribution=distribution)+1);

	# take unperturbed eigenvalues as initial guesses, plus the same with a 
	# real offset from the tunespread
	list_pot_roots=[(np.real(fr)/omega0-Qls,np.imag(fr)/omega0)];
    	if (2*abs(2*bx+bxy)>1e-7): list_pot_roots.append((np.real(fr)/omega0+2*(2*bx+bxy),np.imag(fr)/omega0));

	#print list_pot_roots;

	# actual solving algorithm (note: roots translated from tune shifts to angular frequency shits)
	eigenvaltmp=omega0*roots_complex_list(f,list_pot_roots,tolf=1e-3,tolx=1e-10);
    	#print eigenvaltmp;
	if len(eigenvaltmp)>1:
	    print "Pb in solve_stability_diagram: too many solutions",len(eigenvaltmp),eigenvaltmp;

	if len(eigenvaltmp)>=1:
	    itmp=np.argmin(np.imag(eigenvaltmp));
	    eigenval.append(eigenvaltmp[itmp]+Qls*omega0);# keep here most critical one (arbitrary choice...)
	else:
	    print "Pb in solve_stability_diagram: did not find any solution",fr/omega0,l,len(eigenvaltmp);
	    
    return np.array(eigenval);


def eigenmodesDELPHI_tunespread_converged(nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,Z,freqZ,coefdamper,coefZ,
	omegas,bx,bxy,flag_trapz=0,flagdamperimp=0,d=None,freqd=None,kmax=5,crit=5.e-2,abseps=1.e-3,
	lmaxold=-1,nmaxold=-1,matdamperold=None,matZold=None,distribution='gaussian',kini=12):#,minimag=None,maximag=None,npts_all=1000,npts_zoom=1000):

    ''' computes eigenmodes for increasing matrix size, until we get convergence on 
    the imaginary part of eigenvalues. We include tunespread.
    
    - nx: coupled-bunch mode considered (from 0 to M-1),
    - M: number of equidistant bunches,
    - omegaksi: chromatic angular frequency,
    - omega0: revolution angular frequency,
    - tunefrac: fractional part of the tune,
    - a,b: parameters for Laguerre polynomial decomposition in DELPHI,
    - taub: total bunch length (seconds) or 4*RMS for Gaussian bunches,
    - g: array with coef. of decomposition of initial distribution on Laguerre polynomials,
    - Z, freqZ: Z=damper impedance at the frequencies (NOT ANGULAR) freqZ
    - coefdamper, coefZ: coefficients in fromt resp. of damper and impedance matrix (precomputed beforehand),
    - omegas: synchrotron angular frequency.
    - bx & bxy: detuning coefficients for the transverse tunespread (in units of sigma),
      in the plane of the coherent motion (bx) and in the other plane (bxy)
    - flag_trapz: flag to use trapz method for computation of impedance matrix (if 1),
    - flagdamperimp: flag to use frequency dependent damper gain (if 1)(with d and freqd arrays),
    - d, freqd: d=damper impedance at the frequencies (NOT ANGULAR) freqd - used only if freqdamperimp==1,
    - kmax: number of eigenvalues to make converge,
    - crit: relative error tolerated on imaginary part of eigenvalues,
    - abseps: ~absolute error tolerated on imaginary part of eigenvalues (also used for absolute error determination
          on impedance and damper sums),
    - lmaxold, nmaxold: max azimuthal mode number of max radial mode number (radial mode begins at 0)
          already computed (in matdamperold and matZold),
    - matdamperold ,matZold: damper and impedance matrices already computed,
    - distribution: kind of transverse distribution ('gaussian' or 'parabolic'),
    - kini: number of eigenvalues of the linear system (without Landau damping)
    to initialize the solving algorithm.
    '''

    eigenvalold=np.ones(kmax)*(1.e50+1j*1.e50);
    err=1e50;lmax=1;nmax=1;
    
    # absolute errors tolerated in damper and impedance sums computations 
    absepsd=abseps/max(abs(coefdamper),10);#print "absepsd=",absepsd
    absepsZ=abseps/max(abs(coefZ.real),10);#print "absepsZ=",absepsZ

    flagcrit=True;
    
    while flagcrit:
    
        #t1=ti.clock()
    	
	matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,
		flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=absepsd,lmaxold=lmaxold,
		nmaxold=nmaxold,damperold=matdamperold);

    	matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,tunefrac,a,b,taub,g,Z,freqZ,
		flag_trapz=flag_trapz,abseps=absepsZ,lmaxold=lmaxold,nmaxold=nmaxold,couplold=matZold);

	#t2=ti.clock();
	#print "lmax=",lmax,", nmax=",nmax,", time for computing matrices [seconds]: ",t2-t1;
	
	if ((len(np.where(np.isinf(matdamper))[0])==0)and(len(np.where(np.isnan(matdamper))[0])==0))and((len(np.where(np.isinf(matZ))[0])==0)and(len(np.where(np.isnan(matZ))[0])==0)):

	    if (lmax<16)and(nmax<16):
	    	eigenval=solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,
			bx,bxy,omega0,omegas,distribution=distribution,kini=kini);
	    else:
	    	print "	warning: eigenmodesDELPHI_tunespread_converged: too big matrices, convergence stopped, error=",err,", lmax=",lmax,", nmax=",nmax;
		
	else:
	    print "	warning: eigenmodesDELPHI_tunespread_converged: some nan or inf in matrix, convergence stopped, error=",err,", lmax=",lmax,", nmax=",nmax;
	    break;
	
	ind=np.argsort(np.imag(eigenval));
	
	kmaxtmp=min(kmax,len(eigenval));
	
	err=np.max(np.abs(np.imag(eigenval[ind[:kmaxtmp]])
		-np.imag(eigenvalold[:kmaxtmp]))
		/(np.abs(np.imag(eigenvalold[:kmaxtmp]))+abseps));
	
	flagcrit=(err>crit)or(kmaxtmp<kmax);
		
	for k in range(min(kmax,len(eigenval))): eigenvalold[k]=eigenval[ind[k]];
	
	if (lmax>lmaxold)and(nmax>nmaxold):
	    lmaxold=lmax;nmaxold=nmax;
	    matZold=matZ;matdamperold=matdamper;

	lmax+=1;nmax+=1;


    return eigenval[ind],lmaxold,nmaxold,matdamperold,matZold;
	


def eigenmodesDELPHI_tunespread_converged_scan(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,
	bxscan,bxyscan,M,omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle='proton',flagnorm=0,
	flag_trapz=0,flagdamperimp=0,d=None,freqd=None,kmax=1,kmaxplot=10,
	crit=5.e-2,abseps=1.e-3,flagm0=False,distribution='gaussian',kini=12):
	
    ''' encapsulate DELPHI calculations (with tunespread), with scans on coupled-bunch modes, 
    damper gain, nb of particles, synchrotron tune, damper phase, and detuning coefficients.
    return tuneshifts (complex) for all these parameters scanned
    and also the tuneshifts of the most unstable coupled-bunch modes
    Note: bxscan & bxyscan are done simultaneously if they are of same length
    (then in the resulting table of tuneshifts, the length of bxyscan is 
    replaced by 1)
    '''

    f0=omega0/(2.*np.pi);
    Qfrac=Q-np.floor(Q);
    
    if (len(bxscan)==len(bxyscan)): bxyscanbis=np.array([0.]);# 1D scan simultaneously along bxscan and bxyscan; only the length of bxyscanbis is used then
    else: bxyscanbis=bxyscan; # 2D scan of bxscan and bxyscan

    if flagnorm==0:
	# normalization factor for damper
	dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,Qfrac,a,b,taub,g,
    	    flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=abseps);
	dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
    tuneshift_most=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
    if flagm0:
        tuneshiftm0nx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis)),dtype=complex);
	tuneshiftm0=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis)),dtype=complex);

    for iQp,Qp in enumerate(Qpscan):
    
   	omegaksi=Qp*omega0/eta;

	for inx,nx in enumerate(nxscan):

	    if flagnorm:
		# normalization factor for damper
		dnormfactor=compute_damper_matrix(0,0,nx,M,omegaksi,omega0,Qfrac,a,b,taub,g,
	    	    flagdamperimp=flagdamperimp,d=d,freqd=freqd,abseps=abseps);
		dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

	    lmaxi=-1;nmaxi=-1; # to keep track of maximum lmax and nmax for all modes
   	    lmax=-1;nmax=-1;matZ=None;matdamper=None;

	    for idamp,damp in enumerate(dampscan):

		for iNb,Nb in enumerate(Nbscan):

		    for iomegas,omegas in enumerate(omegasscan):

			for idphase,dphase in enumerate(dphasescan):

			    #print "Nb=",Nb,',lmax=',lmax,', nmax=',nmax;

			    coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,
				    dphase,M,Nb,gamma,Q,particle=particle);

			    for ibx,bx in enumerate(bxscan):
			    
			    	for ibxy,bxy in enumerate(bxyscanbis):
			    
				    # when bxscan and bxyscan have same length, we scan these
				    # 2 tables simultaneously
				    if (len(bxscan)==len(bxyscan)): bxy=bxyscan[ibx];
				    
				    freqshift,lmax,nmax,matdamper,matZ=eigenmodesDELPHI_tunespread_converged(nx,
					    M,omegaksi,omega0,Qfrac,a,b,taub,g,Z,freq,bx,bxy,coefdamper,
					    coefZ,omegas,flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,
					    freqd=freqd,kmax=kmax,crit=crit,abseps=abseps,lmaxold=lmax,
					    nmaxold=nmax,matdamperold=matdamper,matZold=matZ,
					    distribution=distribution,kini=kini);
			    
				    if flagm0:
					# extract mode m=0 tuneshift
					 tuneshiftm0nx[iQp,inx,idamp,iNb,iomegas,idphase,ibx,ibxy]=extract_between_bounds(freqshift/omega0,-omegas/(2.*omega0),omegas/(2.*omega0));

				    if (len(freqshift)>=kmaxplot):
			        	tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,ibx,ibxy,:]=freqshift[:kmaxplot]/omega0;
				    else:
			        	tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,ibx,ibxy,:len(freqshift)]=freqshift[:]/omega0;
			        
				    #lambdax[iNb,:]=freqshift[:kmaxplot]/omegas;
				    lmaxi=max(lmax,lmaxi);nmaxi=max(nmax,nmaxi);	

	# find the most unstable coupled-bunch mode
	for idamp,damp in enumerate(dampscan):

	    for iNb,Nb in enumerate(Nbscan):

		for iomegas,omegas in enumerate(omegasscan):

		    for idphase,dphase in enumerate(dphasescan):

			for ibx,bx in enumerate(bxscan):

			    for ibxy,bxy in enumerate(bxyscanbis):
			    
				if (len(bxscan)==len(bxyscan)): bxy=bxyscan[ibx];
				
				for kmode in range(kmaxplot):

				    inx=np.argmin(np.imag(tuneshiftnx[iQp,:,idamp,iNb,iomegas,idphase,ibx,ibxy,kmode]))
				    if kmode<kmax:
					print "Qp=",Qp,", M=",M,", d=",damp,", Nb=",Nb,", omegas=",omegas,", dphase=",dphase,", bx=",bx,", bxy=",bxy;
					print "   lmaxi=",lmaxi,", nmaxi=",nmaxi,", kmode=",kmode,", Most unstable coupled-bunch mode: ",nxscan[inx];
				    tuneshift_most[iQp,idamp,iNb,iomegas,idphase,ibx,ibxy,kmode]=tuneshiftnx[iQp,inx,idamp,iNb,iomegas,idphase,ibx,ibxy,kmode];

				if flagm0:

				    inx=np.argmax(np.abs(np.real(tuneshiftm0nx[iQp,:,idamp,iNb,iomegas,idphase,ibx,ibxy])))
				    if kmode<kmax:
					print "Qp=",Qp,", M=",M,", d=",damp,", Nb=",Nb,", omegas=",omegas,", dphase=",dphase,", bx=",bx,", bxy=",bxy;
					print "   lmaxi=",lmaxi,", nmaxi=",nmaxi,", kmode=",kmode,", Most unstable coupled-bunch mode: ",nxscan[inx];
				    tuneshiftm0[iQp,idamp,iNb,iomegas,idphase,ibx,ibxy]=tuneshiftm0nx[iQp,inx,idamp,iNb,iomegas,idphase,ibx,ibxy];

    if flagm0: return tuneshift_most,tuneshiftnx,tuneshiftm0;
    else: return tuneshift_most,tuneshiftnx;
    

def eigenmodesDELPHI_tunespread_converged_scan_lxplus(Qpscan,nxscan,dampscan,Nbscan,omegasscan,dphasescan,
	bxscan,bxyscan,M,omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle='proton',
	flagnorm=0,flag_trapz=0,flagdamperimp=0,d=None,freqd=None,
	kmax=1,kmaxplot=10,crit=5.e-2,abseps=1.e-3,flagm0=False,
	distribution='gaussian',kini=12,
	lxplusbatch=None,comment='',queue='1nh',dire=''):
	
    ''' same as eigenmodesDELPHI_tunespread_converged_scan with possibility to launch on lxplus
    lxplusbatch: if None, no use of lxplus batch system
                 if 'launch' -> launch calculation on lxplus on queue 'queue'
                 if 'retrieve' -> retrieve outputs
    comment is used to identify the batch job name and the pickle filename
    dire is the directory where to put/find the result; it should be a path 
    relative to the current directory (e.g. '../') '''
    
    
    import pickle as pick;
    import commands;
    if (len(bxscan)==len(bxyscan)): bxyscanbis=np.array([0.]);# 1D scan simultaneously along bxscan and bxyscan; only the length of bxyscanbis is used then
    else: bxyscanbis=bxyscan; # 2D scan of bxscan and bxyscan

    tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
    tuneshift_most=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
    if flagm0: tuneshiftm0=np.zeros((len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis)),dtype=complex);
    
    if (lxplusbatch==None):
        if flagm0: tuneshift_most,tuneshiftnx,tuneshiftm0=eigenmodesDELPHI_tunespread_converged_scan(Qpscan,
		nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,
		omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=freqd,
		kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
		distribution=distribution,kini=kini);
	
	else: tuneshift_most,tuneshiftnx=eigenmodesDELPHI_tunespread_converged_scan(Qpscan,
		nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,
		omega0,Q,gamma,eta,a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
		flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=freqd,
		kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
		distribution=distribution,kini=kini);

    else:
    
	filename="DELPHI_pickle_"+comment+".dat";
	status,lsres=commands.getstatusoutput('ls '+dire+'out'+filename);
	#print status,lsres
        
	if lxplusbatch.startswith('launch')and(status!=0): # test if output file not already there
    
    	    # pickle file
    	    fileall=open(filename,'w');
    	    listnames=['Qpscan','nxscan','dampscan','Nbscan','omegasscan','dphasescan',
	    	'bxscan','bxyscan','M','omega0','Q','gamma','eta','a','b','taub','g','Z','freq',
		'particle','flagnorm','flag_trapz','flagdamperimp','d','freqd','kmax',
		'kmaxplot','crit','abseps','flagm0','distribution','kini']
    	    for name in listnames: eval('pick.dump('+name+',fileall)');
	    fileall.close();

            # launch calculation on lxplus batch system
            filejob=open('batch'+comment+'.job','w');
            here=os.getcwd()+'/';#print here;
            print >> filejob, "mv "+here+"/"+filename+" .";
            #print >> filejob, "cp "+here+"/DELPHI.py .";
            #print >> filejob, "cp "+here+"/DELPHI_tunespread_script.py .";
	    #print >> filejob, "cp "+here+"/plot_lib.py .";
	    #print >> filejob, "cp "+here+"/io_lib.py .";
	    #print >> filejob, "cp "+here+"/solve_lib.py .";
	    #print >> filejob, "cp "+here+"/string_lib.py .";
	    #print >> filejob, "cp "+here+"/tables_lib.py .";
	    #print >> filejob, "cp "+here+"/C_complex.py .";
 	    #print >> filejob, "cp "+here+"/particle_param.py .";
            #print >> filejob, "chmod +x DELPHI_tunespread_script.py";
            print >> filejob, "DELPHI_tunespread_script.py "+filename+" > out_"+comment;
            print >> filejob, "cp out"+filename+" "+here+dire;
            print >> filejob, "cp out_"+comment+" "+here+dire;
            filejob.close();
            os.system("chmod 744 batch"+comment+".job");
            os.system("bsub"+user_option+" -e error1.out -q "+queue+" batch"+comment+".job");

    	else: # 'retrieve' case
    	    
	    status,lsres=commands.getstatusoutput('ls '+dire+'out'+filename);
	    if status==0:
		# pickle output file
    		fileall=open(dire+'out'+filename,'r');
    		tuneshift_most=pick.load(fileall);
		tuneshiftnx=pick.load(fileall);
		if flagm0: tuneshiftm0=pick.load(fileall);
		fileall.close();
	    
	    else:
	    	# fill tables with nans
		tuneshift_most.fill(np.nan+1j*np.nan);
		tuneshiftnx.fill(np.nan+1j*np.nan);
		if flagm0: tuneshiftm0.fill(np.nan+1j*np.nan);
		print " WARNING !!! no file",filename;
	    
	    os.system("rm -rf LSFJOB_* error1.out batch"+comment+".job");
	
    if flagm0: return tuneshift_most,tuneshiftnx,tuneshiftm0;
    else: return tuneshift_most,tuneshiftnx;
		

def DELPHI_tunespread_wrapper(imp_mod,Mscan,Qpscan,dampscan,Nbscan,omegasscan,dphasescan,
	bxscan,bxyscan,omega0,Qx,Qy,gamma,eta,a,b,taub,g,planes=['x','y'],nevery=1,particle='proton',flagnorm=0,
	flagdamperimp=0,d=None,freqd=None,
	kmax=1,kmaxplot=10,crit=5.e-2,abseps=1.e-3,flagm0=False,
	distribution='gaussian',kini=12,
	lxplusbatch=None,comment='',queue='1nh',dire='',flagQpscan_outside=True):
	
    ''' wrapper of eigenmodesDELPHI_tunespread_converged_scan_lxplus, with more scans.
    imp_mod is an impedance or wake model (see Impedance.py)
    from which one extract the planes given in 'planes'
    Mscan is the scan in number of bunches
    lxplusbatch: if None, no use of lxplus batch system
                 if 'launch' -> launch calculation on lxplus on queue 'queue'
                 if 'retrieve' -> retrieve outputs
    comment is used to identify the batch job name and the pickle filename
    dire is the directory where to put/find the result; it should be a path 
    relative to the current directory (e.g. '../')
    nevery indicates the downsampling (we take 1 frequency every "nevery" 
    frequencies), for the impedance.
    flagQpscan_outside: True to put the Qpscan outside the lxplus batch job,
    False so that it is inside the lxplus batch job.
    '''
    
    from Impedance import test_impedance_wake_comp,impedance_wake
    from copy import deepcopy

    if (len(bxscan)==len(bxyscan)): bxyscanbis=np.array([0.]);# 1D scan simultaneously along bxscan and bxyscan; only the length of bxyscanbis is used then
    else: bxyscanbis=bxyscan; # 2D scan of bxscan and bxyscan

    tuneshiftQp=np.zeros((len(planes),len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
    if flagm0: tuneshiftm0Qp=np.zeros((len(planes),len(Mscan),len(Qpscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis)),dtype=complex);

    for iplane,plane in enumerate(planes):
    	
	# select Zxdip or Zydip
	flagx=int(plane=='x');
	for iw in imp_mod:
	    if test_impedance_wake_comp(iw,flagx,1-flagx,0,0,plane): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);

	for iM,M in enumerate(Mscan):

	    flag_trapz=0; # by default no trapz method

	    if (M==1): nxscan=np.array([0]);flag_trapz=1;
	    else: nxscan=sort_and_delete_duplicates(np.concatenate((np.arange(0,M,M/20),np.arange(M/2-10,M/2+11),
		np.arange(M-10,M),np.arange(0,10))));print "number of coupled-bunch modes=",len(nxscan);

	    if flagQpscan_outside:
		# put loop on Qpscan outside the lxplus job
		for iQp,Qp in enumerate(Qpscan):
		    tuneshiftnx=np.zeros((1,len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscan),kmaxplot),dtype=complex);
		    if flagm0:
			tuneshiftQp[iplane,iM,iQp,:,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[iplane,iM,iQp,:,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			    kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			    distribution=distribution,kini=kini,
			    lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_Qp'+str(Qp)+'_'+plane,
			    queue=queue,dire=dire);

		    else:
			tuneshiftQp[iplane,iM,iQp,:,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus([Qp],
			    nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,omega0,eval('Q'+plane),gamma,eta,
			    a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			    flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			    kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			    distribution=distribution,kini=kini,
			    lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_Qp'+str(Qp)+strnorm[flagnorm]+'_'+plane,
			    queue=queue,dire=dire);
	    
	    else:
		# put loop on Qpscan inside the lxplus job
		tuneshiftnx=np.zeros((len(Qpscan),len(nxscan),len(dampscan),len(Nbscan),len(omegasscan),len(dphasescan),len(bxscan),len(bxyscanbis),kmaxplot),dtype=complex);
		if flagm0:
		    tuneshiftQp[iplane,iM,:,:,:,:,:,:,:,:],tuneshiftnx,tuneshiftm0Qp[iplane,iM,:,:,:,:,:,:,:]=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,omega0,eval('Q'+plane),gamma,eta,
			a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			distribution=distribution,kini=kini,
			lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_'+plane,
			queue=queue,dire=dire);

		else:
		    tuneshiftQp[iplane,iM,:,:,:,:,:,:,:,:],tuneshiftnx=eigenmodesDELPHI_converged_scan_lxplus(Qpscan,
			nxscan,dampscan,Nbscan,omegasscan,dphasescan,bxscan,bxyscan,M,omega0,eval('Q'+plane),gamma,eta,
			a,b,taub,g,Z,freq,particle=particle,flagnorm=flagnorm,
			flag_trapz=flag_trapz,flagdamperimp=flagdamperimp,d=d,freqd=d,
			kmax=kmax,kmaxplot=kmaxplot,crit=crit,abseps=abseps,flagm0=flagm0,
			distribution=distribution,kini=kini,
			lxplusbatch=lxplusbatch,comment=comment+'_'+str(M)+'b_'+strnorm[flagnorm]+'_'+plane,
			queue=queue,dire=dire);
	    

    if flagm0: return tuneshiftQp,tuneshiftm0Qp;
    else: return tuneshiftQp;
    

def extract_between_bounds(tuneshifts,lower,upper):
    ''' extract in a table of complex tuneshifts the "biggest" mode
    that has a real tuneshift between 'lower' and 'upper' '''
    
    ind=np.where((tuneshifts.real>lower)*(tuneshifts.real<upper))[0];
    #print ind,tuneshifts
    
    if (len(ind)>0):
    	i=np.argmax(np.abs(np.real(tuneshifts[ind])));
    	#print i,ind[i],tuneshifts[ind[i]]
	return tuneshifts[ind[i]];
	
    else:
    	print "Extract_between_bounds: no modes between bounds !";
	return 0;


def longdistribution_decomp(taub,typelong='Gaussian'):

    ''' decomposition over Laguerre polynomials of the longitudinal distribution
    for a bunch of length taub (seconds) (4*RMS for Gaussian)
    typelong can be "Gaussian", "parabolicamp", "parabolicline"
    '''
    
    from plot_lib import plot,init_figure,end_figure

    if typelong.startswith('Gaussian'):
    	a=8.;b=8.;g=np.ones(1)*8/(np.pi*taub*taub);
    
    else:

	ntau=10000;deltatau=taub/(2.*ntau);
	tau=np.arange(0,taub/2+deltatau,deltatau); # sampling for integration
	if typelong.startswith('parabolicline'):
            g0=6./(np.pi*taub**2) * np.sqrt(1.-(2.*tau/taub)**2);
            a=9.;b=9.;
        elif typelong.startswith('parabolicamp'):
            g0=8./(np.pi*taub**2) * (1.-(2.*tau/taub)**2);
            a=4.;b=4.;
    
	x=a*(tau/taub)**2;gnorm=real(g0*np.exp(b*(tau/taub)**2));
	tau2=np.arange(0,taub+taub/1000.,taub/1000.);
	gsum=np.zeros((1,len(tau2)));
	i=0;glag=[];glag.append(0.);
	tabtau=np.concatenate(tau,np.array([2.*taub]));
	tabg=np.concatenate(g0,np.array([0.]));

	while (np.trapz(np.abs(np.interp(tau2,tabtau,tabg)-gsum),tau2)
    	    /np.trapz(tabg,tabtau)>5.e-2)and(not(np.isnan(glag[-1]))):

            glag.append(np.trapz(gnorm*np.exp(-x)*laguerre(i,0,x)),x);
            if (not(isnan(glag(i+1)))):
        	gsum=gsum+glag(i+1)*laguerre(i,0,a*(tau2/taub)**2)*np.exp(-b*(tau2/taub)**2);
        	i+=i;
	i-=1;
	if (np.isnan(glag[-1])): glag.pop();
	# check the normalization
	s=np.sum(glag*(1.-a/b)**np.arange(i+1))*taub**2/(2.*b);
	print "check normalization (should be 1): "(s-1./(2.*np.pi))*2.*np.pi
	# plots
	fig,ax=init_figure();
	plot(np.arange(len(glag)),glag,'','b','Laguerre poly. coeff.',ax,0,xlab='Coeff. number');
	fig,ax=init_figure();
	plot(tau,g0,'g0','xb','Long. distribution',ax,0,xlab=r" $ \tau $");
	plot(tau2,gsum,'Sum of Laguerre decomposition','r','Long. distribution',ax,0,xlab=r" $ \tau $ ");

    return g,a,b;


def plot_TMCI(Nbscan,lambdax,ax,part='real',leg='',patcol='xb',
	xlab='Nb particles per bunch',title='',ms=15.,ylim=[-5,5]):

    ''' TMCI plot from lambdax=(Q-Q0)/Qs of all (radial+azimuthal) coherent 
    modes vs intensity Nbscan
    part = 'real' or 'imag'
    '''
    
    from plot_lib import plot,init_figure,end_figure
    
    # make a TMCI kind of plot
    sgn=1;sgnstr='';
    if (part.startswith('real')): strpart='Re';
    if (part.startswith('imag')): strpart='Im';sgn=-1;sgnstr='-'; # invert sign of imaginary part
    
    if (lambdax.ndim>=2):
    	lamsort=sgn*np.sort(getattr(lambdax,part));
    	plot(Nbscan,lamsort[:,:-1],'',patcol,sgnstr+" $ "+strpart+"(Q-Q_0)/Q_s $ ",ax,0,xlab=xlab,ms=ms);
    else:
    	lamsort=np.zeros((len(lambdax),1),dtype=complex);
	lamsort[:,0]=sgn*getattr(lambdax,part);
    
    plot(Nbscan,lamsort[:,-1],leg,patcol,sgnstr+" $ "+strpart+"(Q-Q_0)/Q_s $ ",ax,0,xlab=xlab,ms=ms);
    
    ax.set_title(title);
    ax.set_ylim(ylim);


def find_intensity_threshold(Nbscan,freqshift,thresgrowth=0.,ninterp=1e4):

    ''' find intensity threshold of instability
    when growth rate becomes more than thresgrowth'''
    
    # reinterpolate on finer mesh
    delta=(Nbscan[-1]-Nbscan[0])/ninterp;
    x=np.arange(Nbscan[0],Nbscan[-1]+delta,delta);
    y=np.interp(x,Nbscan,np.imag(freqshift));
    # find threshold
    ind=np.where(y<-thresgrowth)[0];
    if (len(ind)==0): Nbthres=Nbscan[-1];
    else: Nbthres=x[ind[0]];
    
    return Nbthres;
    
    
	
def MOSES_wrapper(rootname,Rres,fres,Qres,Qpscan,Nbscan,omegasscan,omega0,E,alphap,sigmaz,avbeta,
	lmax,nmax,mmin=-3,mmax=3,taumin=-0.5,taumax=0.5,firstline='MOSES input file #',action='launch',
	MOSES_exec='~/DFS/Documents/MOSES4W/MOSES\ application\ for\ Windows/MOSES4W.exe',
	direname_local='DFS',dirname_remote=r"//cernhomeN.cern.ch/"+user[0]+'/'+user,flag_win=True):

    ''' Wrapper to produce easily input files for MOSES (https://oraweb.cern.ch/pls/hhh/code_website.disp_code?code_name=MOSES), from scanned parameters.
    MOSES is a code to compute instabilities, written by Y. Chin (Y. H. Chin. User's Guide for New MOSES Version 2.0 ,
    CERN/LEP-TH/88-05, 1988) that works in a similar way as DELPHI (single-bunch, resonator impedance only, with Landau damping).
    You need to have installed the executable in the path 'MOSES_exec'. This is specially
    designed for a MOSES executable working under Windows, in a directory accessible from here (typically, on a mounted
    DFS file system).
    
    This wrapper, depending on 'action', either writes MOSES input files for scanned parameters and write a Windows batch file
    to launch all the calculations (you still have to launch calculations by hand), or 
    retrieve the MOSES output to get tuneshifts.
    Inputs:
     - rootname: root (including path) of the final input and batch filenames,
     - Rres, fres, Qres: transverse shunt impedance, resonance frequency and quality factor
       for the resonator impedance model. Each (or all) of them can be a list of values, or a scalar.
     - Qpscan: list of chromaticity (Q') to scan,
     - Nbscan: list of number of particles to scan (actually, each value is not used, but only
       the number of values, the first and the last ones, i.e. it makes a linear sampling between
       the first and last number of particles in Nbscan, with min(len(Nbscan),120) points),
     - omegasscan: list of omegas (synchrotron angular frequency in rad/s) to scan,
     - omega0: angular revolution frequency in rad/s,
     - E: energy in eV,
     - alphap: momentum compaction factor,
     - sigmaz: bunch length in m,
     - avbeta: average beta function in m (actually R/Q),
     - lmax: maximum azimuthal mode number to consider (we go from -lmax to +lmax),
     - nmax: number of radial mode -1 (0 means 1 radial mode),
     - mmin & mmax: minimum and maximum for the y-axis of MOSES TMCI plot for the 
       real part of the tune shift / Qs,
     - taumin & taumax: minimum and maximum for the y-axis of MOSES TMCI plot for the 
       imaginary part of the tune shift / Qs,
     - firstline: first line of MOSES input files (some comment),
     - action: 'launch' or 'retrieve': either write input file + Win. bat file to launch 
       all of them easily, or read output and give tuneshifts.
     - MOSES_exec: local path where to find MOSES' executable (Unix path).
     - direname_local: local directory name in which the remote file system is mounted (do not include
     the ~/ or /home/[username]/, just put the bare directory name or path from ~/).
     - dirname_remote: remote directory name to which direname_local corresponds (NOTE: you can keep
     the '/' like this, for Windows paths the replacement by '\' will be done automatically - see next flag).
     - flag_win: True if remote path is a Windows path; then all '/' are replaced by '\'.
    Outputs:
       all modes tuneshifts, tuneshifts of mode 0, and Iscan (arrays)'''
    
    e,m0,c,E0=proton_param();
    
    if (action=='launch'):
    	# Win. batch file name and directory name to put in it
	batname=rootname+'.bat';
	fidbat=open(batname,'wb');
	# create directory if not already there, and copy MOSES executable
	n=rootname.rfind('/');
	os.system("mkdir -p "+rootname[:n]);
	os.system("cp "+MOSES_exec+" "+rootname[:n]);
	# directory name to put in it at each line
	n1=rootname.find(dirname_local)+len(dirname_local);
	dirname=dirname_remote+rootname[n1:];
	if flag_win:
	    n2=dirname.rfind('/');
	    dirname=dirname[:n2].replace("/","\ ").replace(' ','');
    
    # generate list of resonators
    Rreslist=create_list(Rres);
    freslist=create_list(fres,n=len(Rreslist));
    Qreslist=create_list(Qres,n=len(Rreslist));
    if (len(Rreslist)!=len(freslist))or(len(Rreslist)!=len(Qreslist)):
    	print "Pb in MOSES_wrapper: length of resonator models list not all identical !";sys.exit();
    
    # revolution frequency (Hz and MHz)
    f0=omega0/(2.*np.pi);f0MHz=f0*1e-6;
    
    # intensity scan parameters for MOSES in mA, and number of intensities
    Iini=1e3*e*Nbscan[0]*f0;nstep=min(120,len(Nbscan));
    Istep=1e3*e*f0*max((Nbscan[-1]-Nbscan[0])/float(nstep-1),1e8);
    Iscan=np.arange(Iini,Iini+nstep*Istep,Istep);
    
    # bunch length in cm
    sigmazcm=100.*sigmaz;
    
    # energy in GeV
    EGeV=E*1e-9;
    
    tuneshift=np.zeros((len(Rreslist),len(Qpscan),nstep,len(omegasscan),(nmax+1)*(2*lmax+1)),dtype=complex);
    tuneshiftm0=np.zeros((len(Rreslist),len(Qpscan),nstep,len(omegasscan)),dtype=complex);
    
    for ires,R in enumerate(Rreslist):
        
	# shunt impedance in MOhm/m and resonance frequency in MHz
	RMOhm=R*1e-6;frMHz=freslist[ires]*1e-6;Q=Qreslist[ires];
	
    	for iQp,Qp in enumerate(Qpscan):
	
	    for iomegas,omegas in enumerate(omegasscan):
	
		# synchrotron tune
		Qs=omegas/omega0;
		# input/output suffix
		suffix='_R'+float_to_str(RMOhm)+'MOhmm_fr'+float_to_str(frMHz)+'MHz_Q'+float_to_str(Q)+'_Qp'+float_to_str(Qp)+'_Qs'+float_to_str(Qs);

	        if (action=='launch'):
		    # input filename
		    inputname=rootname+suffix+'.dat';
		    # write input
		    fid=open(inputname,'wb');
		    print >> fid, ' '+firstline;
		    print >> fid, ' &MPARM'
		    print >> fid, ' NUS='+fortran_str(Qs)+', ENGY='+fortran_str(EGeV)+', SGMZ='+fortran_str(sigmazcm)+',  BETAC='+fortran_str(avbeta)+',';
		    print >> fid, ' REVFRQ='+fortran_str(f0MHz)+', ALPHA='+fortran_str(alphap)+', CHROM='+fortran_str(Qp)+', SPRD=0.0';
		    print >> fid, ' &END'
		    print >> fid, ' &CPARM'
		    print >> fid, ' CRNT='+fortran_str(Iini)+', STPC='+fortran_str(Istep)+', NCR='+str(nstep)+', NMODF='+str(-lmax)+', NMODE='+str(lmax)+', KRAD='+str(nmax);
		    print >> fid, ' IPRINT=.TRUE., LPLE=.TRUE.';
		    print >> fid, ' &END'
		    print >> fid, ' &IPARM'
		    print >> fid, ' FREQ='+fortran_str(frMHz)+', RS='+fortran_str(RMOhm)+', QV='+fortran_str(Q);
		    print >> fid, ' &END'
		    print >> fid, ' &HPARM'
		    print >> fid, ' MMIN='+str(mmin)+', MMAX='+str(mmax)+', TAUMIN='+fortran_str(taumin)+', TAUMAX='+fortran_str(taumax);
		    print >> fid, ' &END'
		    fid.close();
		    # add one line to .bat file
		    print >> fidbat, 'start /d '+dirname+' MOSES4W.exe '+rootname[n+1:]+suffix+'.dat';
		
		elif (action=='retrieve'):
		    # output filename
		    outputname=rootname+suffix+'.top';
		    # read output
		    linecur=0; # number of initial header lines
		    
		    # real and imag. part of tuneshift
		    for ir,r in enumerate(['real','imag']):
			
			linecur += 29; # number of bla-bla lines
			
			for k in range((nmax+1)*(2*lmax+1)):
		            # read block (for a single mode)
			    s=read_ncol_file(outputname,ignored_rows=linecur);

			    if (np.max(np.abs(s[:,0]-Iscan))>1e-5): print "Pb in MOSES_wrapper: intensity scan incorrect";print np.max(np.abs(s[:,0]-Iscan));

			    tuneshift[ires,iQp,:,iomegas,k] += np.squeeze((1j**ir)*s[:,1]*Qs);
			    
			    # go to next mode block
			    linecur += len(s[:,0])+1;
			
    		    for iNb,I in enumerate(Iscan):
			# extract mode m=0 tuneshift
			tuneshiftm0[ires,iQp,iNb,iomegas]=extract_between_bounds(tuneshift[ires,iQp,iNb,iomegas,:],-Qs/2.,Qs/2.);
			# sort tuneshift (most unstable modes first) (Note: inverted sign for imag. part in MOSES)
			ind=np.argsort(np.imag(tuneshift[ires,iQp,iNb,iomegas,:]));
			tuneshiftnew=tuneshift[ires,iQp,iNb,iomegas,ind[::-1]];
			tuneshift[ires,iQp,iNb,iomegas,:]=tuneshiftnew;


    if (action=='launch'): fidbat.close();print " Now you can execute manually the file",batname;
    
    return tuneshift,tuneshiftm0,Iscan*1e-3/(e*f0);


def long_matching_from_sigmaz(sigmaz,gamma,eta,Qs,R,V,h,particle='proton',flaglinear=False):

    ''' computes delta_p/p0 (sigma) and longitudinal emittance (in eV.s) for a matched bunch
    from:
     - sigmaz = bunch length [m],
     - gamma = relativistic mass factor,
     - eta = slip factor = alphap - 1/gamma^2,
     - Qs = synchrotron tune (at zero amplitude),
     - R = machine radius (circumference / 2pi) [m],
     - V = RF voltage [V],
     - h = RF harmonic number,
     - particle -> 'proton' or 'electron',
     - flaglinear -> True for linear RF bucket (otherwise it does non-linear matching).
    '''

    e,m0,c,E0=eval(particle+'_param()');
    beta=np.sqrt(1.-1./(gamma**2))
    p0=m0*beta*gamma*c;
    E=gamma*E0/e; # total energy in eV

    
    if flaglinear:
    	# linear matching (based on B. Salvant Mathematica notebook)
	delta=Qs*sigmaz/(R*eta); # sigma(delta_p/p0)
	eps0=4*np.pi*p0*beta*sigmaz*delta/e; # long. emittance (eV.s)

    else:
    	# non-linear matching (based on E. Metral Excel sheet)
	taub=4*sigmaz/(beta*c); # total bunch length in seconds (4 times RMS)
	print "total bunch length in ns:",taub*1e9;
    	f0=c*beta/(2*np.pi*R); # rev. frequency
	# long. emittance (eV.s)
	eps0=np.sqrt(V*256*E*R**2/(2*np.pi*np.abs(eta)* h**3 * c**2))*(np.sin(np.pi*h*f0*taub/2.255)/np.sin(np.pi/2.255))**2;
	# RF bucket acceptance (eV.s)
	accept=eps0*(np.sin(np.pi/2.255)/np.sin(np.pi*h*f0*taub/2.255))**2;
	# sigma(delta_p/p0)
	delta=np.pi*h*f0*accept*np.sin(np.pi*h*f0*taub/2)/(8*E*beta**2);
	
    
    return delta,eps0


def Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton'):

    ''' computes Qs (at zero amplitude) from RF parameters:
     - V = RF voltage [V],
     - h = RF harmonic number,
     - gamma = relativistic mass factor,
     - eta = slip factor = alphap - 1/gamma^2,
     - phis = synchrotron phase [rad],
     - particle -> 'proton' or 'electron'.
    '''

    e,m0,c,E0=eval(particle+'_param()');
    
    beta=np.sqrt(1.-1./(gamma**2))
    p0=m0*beta*gamma*c;

    Qs=np.sqrt(e*V*eta*h*np.cos(phis)/(2*np.pi*beta*c*p0));
    
    return Qs;
    

def eta_from_Qs_RF_param(Qs,V,h,gamma,phis=0.,particle='proton'):

    ''' computes eta (slip factor) from RF parameters:
     - Qs = synchrotron tune,
     - V = RF voltage [V],
     - h = RF harmonic number,
     - gamma = relativistic mass factor,
     - eta = slip factor = alphap - 1/gamma^2,
     - phis = synchrotron phase [rad],
     - particle -> 'proton' or 'electron'.
    '''

    e,m0,c,E0=eval(particle+'_param()');
    
    beta=np.sqrt(1.-1./(gamma**2))
    p0=m0*beta*gamma*c;

    eta=Qs**2*2*np.pi*beta*c*p0/(e*V*h*np.cos(phis));
    
    return eta;
    

def write_CFG_HEADTAIL(cfgname,particle='proton',Nb=1.e11,betax=65.976,betay=71.526,
	sigmaz=0.0937,emitx=2e-6,emity=2e-6,delta=1.4435e-4,Qs=2.34243e-3,alphap=3.225e-4,
	circ=26658.883,gamma=4263.16,nkick=1,nturns=200000,pipex=0.05,pipey=0.05,
	nsigmaz=2,Qx=64.31,Qy=59.32,Qpx=0,Qpy=0,isyn=4,nMP_per_sl=5000,nbunch=1,
	nsl=200,spacing=20,ibunchtb=0,iwake=1,ipipe=8,ntwake=20,fresT=1e9,QT=1.,
	RT=0.,fresL=200e6,QL=140.,RL=0.,condRW=1.4e6,lenRW=26658.883,ibeta=0,
	iwaketb=6,ispace=0,smooth=3,kickswitch=1,xkick=1,ykick=1.,zkick=0.,iamp=1,
	icoupl=0,coupl=0.0015,dispx=0,isext=0,sext_str=-0.254564,disp_sext=2.24,
	iloss=2,Qsecx=0,Qsecy=0,nturns_betn=1,start_turn=199000,end_turn=199100,
	VRF=12e6,h=35640,VRF2_start=0.,VRF2_end=0.,h2=18480,phase=0.,start_RF2=2000,
	end_RF2=3000,prate=0.,alphap_sec=0.,phase_shift_max=1,octfoc=0.,octdefoc=0.,
	dratex=0.02,dratey=0.02,nMP_prb=500,ipre=1,idampbeam=0):
	
    ''' Write an HEADTAIL input file (.cfg)
    Parameters are the same (and in the same order) as in the .cfg file, except that:
    - units are all SI except prate in GeV/c/s (but nsigmaz, xkick and ykick are still in number of sigmas,
    and spacing in RF buckets),
    - number of macroparticles=nMP_per_sl*nsl.
    cfgname contains the name of the input file (with extension).
    Default values are typical of 2012 LHC operation at 4TeV (except detuning and Qsec).
    '''
    
    nMP=nMP_per_sl*nsl;
    emitx *= 1e6;
    emity *= 1e6;
    fresT /= 1e9;
    fresL /= 1e6;
    RT /= 1e6;
    RL /= 1e6;
    if particle.startswith('proton'): ipart=1;
    elif particle.startswith('electron'): ipart=2;
    
    fid=open(cfgname,'wb');
    print >> fid, 'Flag_for_bunch_particles_(1->protons_2->positrons_3&4->ions):   ',ipart;
    print >> fid, 'Number_of_particles_per_bunch:                                  ',Nb;
    print >> fid, 'Horizontal_beta_function_at_the_kick_sections_[m]:              ',betax;
    print >> fid, 'Vertical_beta_function_at_the_kick_sections_[m]:                ',betay;
    print >> fid, 'Bunch_length_(rms_value)_[m]:                                   ',sigmaz;
    print >> fid, 'Normalized_horizontal_emittance_(rms_value)_[um]:               ',emitx;
    print >> fid, 'Normalized_vertical_emittance_(rms_value)_[um]:                 ',emity;
    print >> fid, 'Longitudinal_momentum_spread:                                   ',delta;
    print >> fid, 'Synchrotron_tune:                                               ',Qs;
    print >> fid, 'Momentum_compaction_factor:                                     ',alphap;
    print >> fid, 'Ring_circumference_length_[m]:                                  ',circ;
    print >> fid, 'Relativistic_gamma:                                             ',gamma;
    print >> fid, 'Number_of_kick_sections:                                        ',nkick;
    print >> fid, 'Number_of_turns:                                                ',nturns;
    print >> fid, 'Horizontal_semiaxis_of_beam_pipe_[m]                            ',pipex;
    print >> fid, 'Vertical_semiaxis_of_beam_pipe_[m]                              ',pipey;
    print >> fid, 'Longitud_extension_of_the_bunch_(+/-N*sigma_z)                  ',nsigmaz;
    print >> fid, 'Horizontal_tune:                                                ',Qx;
    print >> fid, 'Vertical_tune:                                                  ',Qy;
    print >> fid, "Horizontal_chromaticity_[Q'x]:                                  ",Qpx;
    print >> fid, "Vertical_chromaticity_[Q'y]:                                    ",Qpy;
    print >> fid, 'Flag_for_synchrotron_motion:                                    ',isyn;
    print >> fid, 'Number_of_macroparticles_per_bunch:                             ',nMP;
    print >> fid, 'Number_of_bunches:                                              ',nbunch;
    print >> fid, 'Number_of_slices_in_each_bunch:                                 ',nsl;
    print >> fid, 'Spacing_between_consecutive_bunches_centroids_[RF_bkts]:        ',spacing;
    print >> fid, 'Switch_for_bunch_table:                                         ',ibunchtb;
    print >> fid, 'Switch_for_wake_fields:                                         ',iwake;
    print >> fid, 'Switch_for_pipe_geometry_(0->round_1->flat):                    ',ipipe;
    print >> fid, 'Number_of_turns_for_the_wake:                                   ',ntwake;
    print >> fid, 'Res_frequency_of_broad_band_resonator_[GHz]:                    ',fresT;
    print >> fid, 'Transverse_quality_factor:                                      ',QT;
    print >> fid, 'Transverse_shunt_impedance_[MOhm/m]:                            ',RT;
    print >> fid, 'Res_frequency_of_longitudinal_resonator_[MHz]:                  ',fresL;
    print >> fid, 'Longitudinal_quality_factor:                                    ',QL;
    print >> fid, 'Longitudinal_shunt_impedance_[MOhm]:                            ',RL;
    print >> fid, 'Conductivity_of_the_resistive_wall_[1/Ohm/m]:                   ',condRW;
    print >> fid, 'Length_of_the_resistive_wall_[m]:                               ',lenRW;
    print >> fid, 'Switch_for_beta:                                                ',ibeta;
    print >> fid, 'Switch_for_wake_table:                                          ',iwaketb;
    print >> fid, 'Flag_for_the_transverse_space_charge:                           ',ispace;
    print >> fid, 'Smoothing_order_for_longitudinal_space_charge:                  ',smooth;
    print >> fid, 'Switch_for_initial_kick:                                        ',kickswitch;
    print >> fid, 'x-kick_amplitude_at_t=0_[sigmas]:                               ',xkick;
    print >> fid, 'y-kick_amplitude_at_t=0_[sigmas]:                               ',ykick;
    print >> fid, 'z-kick_amplitude_at_t=0_[m]:                                    ',zkick;
    print >> fid, 'Switch_for_amplitude_detuning:                                  ',iamp;
    print >> fid, 'Linear_coupling_switch(1->on_0->off):                           ',icoupl;
    print >> fid, 'Linear_coupling_coefficient_[1/m]:                              ',coupl;
    print >> fid, 'Average_dispersion_function_in_the_ring_[m]:                    ',dispx;
    print >> fid, 'Sextupolar_kick_switch(1->on_0->off):                           ',isext;
    print >> fid, 'Sextupole_strength_[1/m^2]:                                     ',sext_str;
    print >> fid, 'Dispersion_at_the_sextupoles_[m]:                               ',disp_sext;
    print >> fid, 'Switch_for_losses_(0->no_losses_1->losses):                     ',iloss;
    print >> fid, "Second_order_horizontal_chromaticity_(Qx''):                    ",Qsecx;
    print >> fid, "Second_order_vertical_chromaticity_(Qy''):                      ",Qsecy;
    print >> fid, 'Number_of_turns_between_two_bunch_shape_acquisitions:           ',nturns_betn;
    print >> fid, 'Start_turn_for_bunch_shape_acquisitions:                        ',start_turn;
    print >> fid, 'Last_turn_for_bunch_shape_acquisitions:                         ',end_turn;
    print >> fid, 'Main_rf_voltage_[V]:                                            ',VRF;
    print >> fid, 'Main_rf_harmonic_number:                                        ',h;
    print >> fid, 'Initial_2nd_rf_voltage_[V]:                                     ',VRF2_start;
    print >> fid, 'Final_2nd_rf_cavity_voltage_[V]:                                ',VRF2_end;
    print >> fid, 'Harmonic_number_of_2nd_rf:                                      ',h2;
    print >> fid, 'Relative_phase_between_cavities:                                ',phase;
    print >> fid, 'Start_turn_for_2nd_rf_ramp:                                     ',start_RF2;
    print >> fid, 'End_turn_for_2nd_rf_ramp:                                       ',end_RF2;
    print >> fid, 'Linear_Rate_of_Change_of_Momentum_[GeV/c/sec]:                  ',prate;
    print >> fid, 'Second_Order_Momentum_Compaction_Factor:                        ',alphap_sec;
    print >> fid, 'Max_phase_shift_delay_after_transition_crossing_[turns]:        ',phase_shift_max;
    print >> fid, 'LHC_focusing_octupoles_current_[A]:                             ',octfoc;
    print >> fid, 'LHC_defocusing_octupoles_current_[A]:                           ',octdefoc;
    print >> fid, 'Horizontal_damper_rate_[inverse_of_nb_damping_turns]:           ',dratex;
    print >> fid, 'Vertical_damper_rate_[inverse_of_nb_damping_turns]:             ',dratey;
    print >> fid, 'Number_of_macroparticles_in_prb_file:                           ',nMP_prb;
    print >> fid, 'Flag_for_wake_pretreatment_(0->no_1->pretreatment)              ',ipre;
    print >> fid, 'Flag_for_damping_beam_average_(1->yes_0->bunch_by_bunch_damper) ',idampbeam;
    
    fid.close();
    
    return;
    

class DELPHI_calc(object):
    ''' class for DELPHI computations (NOT USED AND NOT FINISHED)'''
     
    def __init__(self,machine='LHC',M=1,nx=0,omegaksi=0.,circ=26658.883,Q=64.31,a=8.,b=8.,
     	taub=1.2e-9,gamma=7460.52,particle='proton',g=None,Z=None,freqZ=None,flag_trapz=0,
	flagdamperimp=0,d=None,freqd=None,kmax=5,crit=5.e-2,abseps=1.e-3,lmax=-1,nmax=-1,
	matdamper=None,matZ=None,flageigenvect=False,flagnormscan=np.array([0.]),
	dampscan=np.array([0.]),dphasescan=np.array([0.]),Nbscan=np.array([1.5e11]),
	Qsscan=np.array([2.e-3])):
	# default values typical of LHC (single-bunch) at 7 TeV
	
	self.machine=machine;
	self.M=M;
	self.nx=nx;
	self.omegaksi=omegaksi;
	self.circ=circ;
	self.Q=Q;
	self.a=a;
	self.b=b;
	self.taub=taub;
	self.gamma=gamma;
	self.particle=particle;
	self.g=g;
	self.Z=Z;
	self.freqZ=freqZ;
	self.flag_trapz=flag_trapz;
	self.flagdamperimp=flagdamperimp;
	self.d=d;
	self.freqd=freqd;
	self.kmax=kmax;
	self.crit=crit;
	self.abseps=abseps;
	self.lmax=lmax;
	self.nmax=nmax;
	self.matdamper=matdamper;
	self.matZ=matZ;
	self.flageigenvect=flageigenvect;
	self.flagnormscan=flagnormscan;
	self.dampscan=dampscan;
	self.dphasescan=dphasescan;
	self.Nbscan=Nbscan;
	self.Qsscan=Qsscan;
	
	


if __name__ == "__main__":


    machine='VEPP';
    e=1.602176487e-19; # elementary charge
    m0=9.10938e-31; # electron mass in kg
    c=299792458; # speed of light
    E0=m0*c**2 # rest energy

    E=1.8e9; # injection energy=1.8 GeV
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    sigmaz=7.5e-2; # RMS bunch length (m)
    bunchlength=4.*sigmaz; # full length in m
    VEPPlength=beta*c/8.377e5; # total circumference in m
    circ=VEPPlength;
    R=circ/(2.*np.pi) # machine radius
    taub=bunchlength/(beta*c);
    Qx=7.62;
    Qs=0.025;
    alphap=0.01645; # momentum compaction factor
    M=1; # number of bunches
    max_azim=2;nrad=2;kmax=10;nxmin=0;nxmax=-1;nxadded=[];
    f0=c*beta/circ # rev. frequency

    beta1=15;betaav=R/Qx

    Qpx=-7.5;

    # parameters for intensity scan
    deltaI=0.5;maxI=50;
    deltaI_0=0.2;maxI_0=20;

    dphase=0; # additional phase applied to the damper (rad)

    # one broad band model (NB: shunt impedances to be multiplied by beta(location)/(R/Q) )
    f=np.concatenate((10.**np.arange(-1,7),np.arange(1.e7,1.001e10,1.e7),10.**np.arange(10.1,13.1,0.1),
    	10.**np.arange(14,16)));
    Zx=np.zeros((len(f),2));
    #print f,len(f)
    R1=2.5e6*beta1/betaav;f1=R*f0*0.795/sigmaz;Q1=1.;
    Zt=(R1*f1/f)/(1.-1j*Q1*(f1/f-f/f1));
    Zx[:,0]=np.real(Zt);Zx[:,1]=np.imag(Zt);
    model='_1BB_Karliner-Popov';

    eta=alphap-1./(gamma*gamma);
    omegaksi=Qpx*2*np.pi*f0/eta;
    g=np.zeros(1);
    a=8.;b=8;g[0]=8./(np.pi*taub*taub);
    
    # damper model
    fd=np.concatenate((np.arange(2.e4,2.002e7,2.e4),np.arange(2.1e7,1.001e9,1e6),np.arange(1.1e9,1.01e10,1.e8)));
    fd=np.concatenate((np.flipud(-fd),np.array([0.]),fd));
    Zd=np.zeros((len(fd),2));
    #print fd,len(fd)
    L0=0.4;L1=0.1;L2=0.2;tauf=800*L1/c;print tauf*2*np.pi*f0 # set of parameters with bare tune (including integer part)
    dimp=damper_imp_Karliner_Popov(L0,L1,L2,tauf,R,Qx,f0,fd);
    Zd[:,0]=np.real(dimp);Zd[:,1]=np.imag(dimp);

    #mat=compute_impedance_matrix(2,1,0,M,omegaksi,2.*np.pi*f0,Qx-np.floor(Qx),a,b,taub,g,Zx,f,1);
    mat=compute_damper_matrix(2,1,0,M,omegaksi,2.*np.pi*f0,Qx-np.floor(Qx),a,b,taub,g,Zd,fd,flagdamperimp=1);

    print mat;

    #taub=1.0006923259e-9;
    #DArray1=c_double * 1;
    #g=DArray1(8/(np.pi*taub*taub));

    #prototype = CFUNCTYPE(c_double, c_int,c_int,c_double,c_double,c_double, c_long, DArray1)
    #Gln = prototype(('Gln', libDELPHI))

    #print Gln(0,0,-1.e9,8/(taub*taub), taub, 0, g);


# example:
# arr = (c_int * len(pyarr))(*pyarr)
