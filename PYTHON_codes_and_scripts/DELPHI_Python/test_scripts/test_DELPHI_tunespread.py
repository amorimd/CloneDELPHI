#!/usr/bin/python

import sys
if len(sys.argv)>2: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=str(sys.argv[2]);
elif len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);lxplusbatchDEL=None;
else: lxplusbatchImp=None;lxplusbatchDEL=None;
print lxplusbatchImp,lxplusbatchDEL;   

import commands
out=commands.getoutput("hostname")
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from string import *
import time
import numpy as np
from copy import deepcopy
import pylab,os,re
sys.path.append("../../LHC_impedance_and_scripts/")
from plot_lib import plot,init_figure,end_figure,plot2D
from io_lib import *
from tables_lib import select_in_table
from particle_param import *
from Impedance import *
from DELPHI import *
from LHC_param import LHC_param
from LHC_imp import *
from HLLHC_imp import *
from LHC_coll_imp import *
import time as ti  


def integ_line(a,b,cst,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'):

    # integrate along a line (horizontal or vertical) the determinant function 
    # (from a to b, the line being at the constant position cst - a, b and cst 
    # can be real or imaginary -> that determinaes if the line is hor. or ver.))
    
    from scipy import integrate as inte;

    if ((a.real==0)and(b.real==0))and(cst.imag==0):
    	fact_imag=1j;
    elif ((a.imag==0)and(b.imag==0))and(cst.real==0):
    	fact_imag=1;
    else:
    	print "Pb in integ_line: a,b and cst not consistent !";sys.exit();
    
    factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
    f_real=(lambda x: factnorm/determinantDELPHI_tunespread(fact_imag*x+cst,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real);
    f_imag=(lambda x: factnorm/determinantDELPHI_tunespread(fact_imag*x+cst,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag);

    int_real,err=inte.quadrature(f_real,fact_imag*np.conjugate(a),fact_imag*np.conjugate(b),tol=1.e-7,maxiter=1000,vec_func=False);
    int_imag,err=inte.quadrature(f_imag,fact_imag*np.conjugate(a),fact_imag*np.conjugate(b),tol=1.e-7,maxiter=1000,vec_func=False);

    integ=int_real+1j*int_imag;
    
    return integ;


def test_residue(ReQmin,ReQmax,ImQmin,ImQmax,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'):

    from scipy import integrate as inte;

    factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
    fimin_real=(lambda x: factnorm/determinantDELPHI_tunespread(x+1j*ImQmin,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real);
    fimax_real=(lambda x: factnorm/determinantDELPHI_tunespread(x+1j*ImQmax,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real);
    frmin_real=(lambda x: factnorm/determinantDELPHI_tunespread(ReQmin+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real);
    frmax_real=(lambda x: factnorm/determinantDELPHI_tunespread(ReQmax+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real);
    fimin_imag=(lambda x: factnorm/determinantDELPHI_tunespread(x+1j*ImQmin,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag);
    fimax_imag=(lambda x: factnorm/determinantDELPHI_tunespread(x+1j*ImQmax,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag);
    frmin_imag=(lambda x: factnorm/determinantDELPHI_tunespread(ReQmin+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag);
    frmax_imag=(lambda x: factnorm/determinantDELPHI_tunespread(ReQmax+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag);

    A_real,err=inte.quadrature(fimin_real,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    A_imag,err=inte.quadrature(fimin_imag,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    B_real,err=inte.quadrature(frmax_real,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    B_imag,err=inte.quadrature(frmax_imag,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    C_real,err=inte.quadrature(fimax_real,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    C_imag,err=inte.quadrature(fimax_imag,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    D_real,err=inte.quadrature(frmin_real,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
    D_imag,err=inte.quadrature(frmin_imag,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);

    integ=A_real+1j*A_imag+B_real+1j*B_imag-C_real-1j*C_imag-D_real-1j*D_imag;
    
    return integ;


if __name__ == "__main__":

    e,m0,c,E0=proton_param();
    root_LHC='../../LHC_impedance_and_scripts';
    results_dir='../../../DELPHI_results/LHC/test';
    os.system("mkdir -p "+results_dir);
    
    # fixed parameters
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h=LHC_param(E0,E=4e12);

    bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-100.,100.);
    print bx,bxy;
    
    if False:
	stabx_MATLAB=read_ncol_file('/home/nmounet/Documents/Stability_diagrams_LHC/stabx_m100A_Gaussian_eps2_noQsec_axx_m4p1e-5_axy_2p9e-5.dat',ignored_rows=1);

	Qscan_gau=10.**np.arange(-9,-2,0.05);
	Qscan_gau=np.hstack((-Qscan_gau[::-1],Qscan_gau));
	Qscan_par=np.arange(-1e-3,1e-3,1e-7);

	stabx_Python=np.zeros(len(Qscan_gau),dtype=complex);
	for iQ,Q in enumerate(Qscan_gau):
	    stabx_Python[iQ]=-1./dispersion_integral_oct_2D(Q,bx,bxy,distribution='gaussian')
	    #print Q,stabx_Python[iQ];

	fig,ax=init_figure();
	ax.plot(stabx_MATLAB[:,0],stabx_MATLAB[:,1],'b');
	ax.plot(np.real(stabx_Python),-np.imag(stabx_Python),'xr');
	end_figure(fig,ax);
	# nice !

	stabx_par=np.zeros(len(Qscan_par),dtype=complex);
	#stabx_par2=np.zeros(len(Qscan_par),dtype=complex);
	for iQ,Q in enumerate(Qscan_par):
	    stabx_par[iQ]=-1./dispersion_integral_oct_2D(Q,bx,bxy,distribution='parabolic')
	    #stabx_par2[iQ]=-1./dispersion_integral_oct_2D2(Q,bx,bxy,distribution='parabolic')

	bxm,bxym=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,100.,-100.);

	stabx_gaum=np.zeros(len(Qscan_gau),dtype=complex);
	for iQ,Q in enumerate(Qscan_gau):
	    stabx_gaum[iQ]=-1./dispersion_integral_oct_2D(Q,bxm,bxym,distribution='gaussian')

	stabx_parm=np.zeros(len(Qscan_par),dtype=complex);
	#stabx_parm2=np.zeros(len(Qscan_par),dtype=complex);
	for iQ,Q in enumerate(Qscan_par):
	    stabx_parm[iQ]=-1./dispersion_integral_oct_2D(Q,bxm,bxym,distribution='parabolic')
	    #stabx_parm2[iQ]=-1./dispersion_integral_oct_2D2(Q,bxm,bxym,distribution='parabolic')

	fig,ax=init_figure();
	ax.plot(np.real(stabx_Python),-np.imag(stabx_Python),'--r',label='Gaussian, -100 A');
	ax.plot(np.real(stabx_par),-np.imag(stabx_par),'r',label='Parabolic, -100 A');
	#ax.plot(np.real(stabx_par2),-np.imag(stabx_par2),'xr',label='New Parabolic, -100 A');
	ax.plot(np.real(stabx_gaum),-np.imag(stabx_gaum),'--b',label='Gaussian, 100 A');
	ax.plot(np.real(stabx_parm),-np.imag(stabx_parm),'b',label='Parabolic, 100 A');
	#ax.plot(np.real(stabx_parm2),-np.imag(stabx_parm2),'xb',label='New Parabolic, 100 A');
	end_figure(fig,ax);
	
	pylab.show()
	# all is OK
    
    if False:
	# test computation of dispersion integral in some limit cases
	Ioct=np.hstack((-10.**np.arange(2,-4.1,-0.1),10.**np.arange(-4,2.1,0.1)));
	distribution='parabolic';Q=3e-4-1e-5*1j;
	bx=np.zeros(len(Ioct));bxy=np.zeros(len(Ioct));Idisp=np.zeros(len(Ioct),dtype=complex);
	for io,o in enumerate(Ioct):
	    bx[io],bxy[io]=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-o,o);
	    Idisp[io]=dispersion_integral_oct_2D(Q,bx[io],bxy[io],distribution=distribution)
	
	fig,ax=init_figure();
	ax.plot(Ioct,np.real(1./Idisp),'b',label='Inverse of dispersion integral');
	ax.plot(0,-Q.real,'xr',label='Limit value',ms=12);
	end_figure(fig,ax);
	fig,ax=init_figure();
	ax.plot(Ioct,np.imag(1./Idisp),'b',label='Inverse of dispersion integral');
	ax.plot(0,-Q.imag,'xr',label='Limit value',ms=12);
	end_figure(fig,ax);
	pylab.show()
	# it works fine (all distributions, all kind of signs for Q)
	
    if False:
	# 2D plot of 1/dispersion integral
	distribution='gaussian';
	ReQmin=-Qs;ReQmax=Qs;
	ImQmin=-Qs/10.;ImQmax=0.;
	npts=1000;
	bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-100.,100.);
    	real_range=np.linspace(ReQmin,ReQmax,npts+1);
    	imag_range=np.linspace(ImQmin,ImQmax,npts+1);
	ftable_real=np.zeros((npts+1,npts+1));
	ftable_imag=np.zeros((npts+1,npts+1));

	for ir,r in enumerate(real_range):
	    for im,m in enumerate(imag_range):
		ftable_real[im,ir]=np.real(1./dispersion_integral_oct_2D(r+1j*m,bx,bxy,distribution=distribution));
		ftable_imag[im,ir]=np.imag(1./dispersion_integral_oct_2D(r+1j*m,bx,bxy,distribution=distribution));
	
	fig,ax=init_figure();
	plot2D(ftable_real,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='');
	ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	end_figure(fig,ax);
	
	fig,ax=init_figure();
	plot2D(ftable_imag,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='');
	ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	end_figure(fig,ax);

	pylab.show()
	
    
    if False:
	# test DELPHI with tunespread
	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

	beam='1';
	param_filename_coll=root_LHC+'/Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
	settings_filename_coll=root_LHC+'/Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';

	imp_mod,wake_mod=LHC_imp_model_v1(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll=None,dire=root_LHC+"/LHC_elements/",commentcoll='_2012_v2',direcoll='Coll_2012_v2/',lxplusbatch='retrieve',
	    beam=beam,squeeze='0p6m_3m_0p6m_3m',wake_calc=False,ftypescan=0,nflog=100,zpar=z_param(),
	    flagplot=False,root_result=results_dir,commentsave='');

	# longitudinal distribution initialization
	g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

	nevery=1;lmax=2;nmax=4;M=1;nx=0;
	Qp=5;damp=0.;
	Qp=0;damp=0.02;
	Qp=15;damp=0.02;
	omegaksi=Qp*omega0/eta;Nb=3e11;maximag=Qs;
	dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,Qxfrac,a,b,taub,g);
	dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

	bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-0.,0.);
	#print bx,bxy;

	# select Zxdip
	for iw in imp_mod:
	    if test_impedance_wake_comp(iw,1,0,0,0,'x'): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);

	#fig,ax=init_figure();
	#ax.loglog(freq,Z[:,0],'b',label='Real');
	#ax.loglog(freq,Z[:,1],'r',label='Imag');
	#end_figure(fig,ax);
	#pylab.show()

	t1=ti.clock()
	coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='proton');

	matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g);

	matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,Z,freq,flag_trapz=1);
	t2=ti.clock();
	print "time for matrix computation [seconds]: ",t2-t1;

	t1=ti.clock()
	freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=False);
	t2=ti.clock();
	print "time for solving eigenvalue problem [seconds]: ",t2-t1;

	ind=np.argsort(np.imag(freqshift));
	print freqshift[ind[:12]]/omega0;
	
    if False:
	from scipy import integrate as inte;
	ReQmin=-2.5*Qs;ReQmax=2.5*Qs;ReQmax=-2*Qs-3e-5;print ReQmin,ReQmax;
	ImQmin=-Qs;ImQmax=2*freqshift[ind[0]].imag/omega0;ImQmax=0.;print ImQmin,ImQmax;
	factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'));
	fimin_real=(lambda x: factnorm/determinantDELPHI_tunespread_array(x+1j*ImQmin,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').real);
	fimax_real=(lambda x: factnorm/determinantDELPHI_tunespread_array(x+1j*ImQmax,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').real);
	frmin_real=(lambda x: factnorm/determinantDELPHI_tunespread_array(ReQmin+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').real);
	frmax_real=(lambda x: factnorm/determinantDELPHI_tunespread_array(ReQmax+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').real);
	fimin_imag=(lambda x: factnorm/determinantDELPHI_tunespread_array(x+1j*ImQmin,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').imag);
	fimax_imag=(lambda x: factnorm/determinantDELPHI_tunespread_array(x+1j*ImQmax,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').imag);
	frmin_imag=(lambda x: factnorm/determinantDELPHI_tunespread_array(ReQmin+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').imag);
	frmax_imag=(lambda x: factnorm/determinantDELPHI_tunespread_array(ReQmax+1j*x,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').imag);
	
	npts=1000;
    	real_range=np.linspace(ReQmin,ReQmax,npts+1);
    	imag_range=np.linspace(ImQmin,ImQmax,npts+1);
	fig,ax=init_figure();
	plot(real_range,fimin_real(real_range),'Re[1/f], along Im='+str(ImQmin),'b','1/f',ax,0,xlab='Re[Q]');
	plot(real_range,fimin_imag(real_range),'Im[1/f], along Im='+str(ImQmin),'--b','1/f',ax,0,xlab='Re[Q]');
	plot(real_range,fimax_real(real_range),'Re[1/f], along Im='+str(ImQmax),'r','1/f',ax,0,xlab='Re[Q]');
	plot(real_range,fimax_imag(real_range),'Im[1/f], along Im='+str(ImQmax),'--r','1/f',ax,0,xlab='Re[Q]');
	ax.set_xlim([ReQmin,ReQmax]);
	end_figure(fig,ax);
	
	fig,ax=init_figure();
	plot(imag_range,frmin_real(imag_range),'Re[1/f], along Re='+str(ReQmin),'b','1/f',ax,0,xlab='Im[Q]');
	plot(imag_range,frmin_imag(imag_range),'Im[1/f], along Re='+str(ReQmin),'--b','1/f',ax,0,xlab='Im[Q]');
	plot(imag_range,frmax_real(imag_range),'Re[1/f], along Re='+str(ReQmax),'r','1/f',ax,0,xlab='Im[Q]');
	plot(imag_range,frmax_imag(imag_range),'Im[1/f], along Re='+str(ReQmax),'--r','1/f',ax,0,xlab='Im[Q]');
	ax.set_xlim([ImQmin,ImQmax]);
	end_figure(fig,ax);
	
	pylab.show()
	
	print "Calcul A";
	A_real,err=inte.quadrature(fimin_real,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	A_imag,err=inte.quadrature(fimin_imag,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	print "Calcul B";
	B_real,err=inte.quadrature(frmax_real,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	B_imag,err=inte.quadrature(frmax_imag,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	print "Calcul C";
	C_real,err=inte.quadrature(fimax_real,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	C_imag,err=inte.quadrature(fimax_imag,ReQmin,ReQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	print "Calcul D";
	D_real,err=inte.quadrature(frmin_real,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	D_imag,err=inte.quadrature(frmin_imag,ImQmin,ImQmax,tol=1.e-7,maxiter=1000,vec_func=False);
	print A_real+1j*A_imag+B_real+1j*B_imag-C_real-1j*C_imag-D_real-1j*D_imag;
	print A_real,A_imag,B_real,B_imag,C_real,C_imag,D_real,D_imag;
	
	# it's slow and it does not work...

    if False:
	kmax=6;
	#na=1;nz=1;errmax=1.;
	#frold=np.ones(kmax,dtype=complex)*1e50;fr=np.zeros(kmax,dtype=complex);
	#minimag=-1.5e-5;maximag=0;
	#while (errmax>0.01)and(False):
	#    freqshift_spread,freqshift_spread_pot=solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,
	#    	minimag=minimag,maximag=maximag,distribution='gaussian',npts_all=na*2*lmax,npts_zoom=nz*10);
	#    if nz==1: frtmp=deepcopy(freqshift_spread);
	#    else: frtmp=sort_and_delete_duplicates(np.hstack((fr[:le],freqshift_spread)),tolerance=1e-4);
	#    le=min(kmax,len(frtmp));
	#    indtmp=np.argsort(np.imag(frtmp));
	#    fr[:le]=frtmp[indtmp[:le]];
	#    errmax=np.max(np.abs(np.imag(fr-frold)/np.imag(fr)));
	#    print "na=",na,",nz=",nz,",minimag=",minimag,",maximag=",maximag,",tuneshifts:",fr/omega0;
	#    frold[:le]=fr[:le];
	#    na+=1;nz+=3;
	#    if fr[0].imag<0: minimag=10*fr[0].imag/omega0;
	#    if le<kmax: maximag=max(fr[le-1].imag/2.,2*fr[le-1].imag)/omega0;
	#    else: maximag=fr[le-1].imag/omega0;
	#print fr/omega0;
	
	Ioctscan=-np.array([0.,5.,10.,30.,50.,100.]);
	npts_zoom=100;
	distribution='gaussian';kini=12;
	#distribution='parabolic';
	for Ioct in Ioctscan:
	    bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-Ioct,Ioct);
	    print bx,bxy;
	    t1=ti.clock()
	    #freqshift_spread,freqshift_spread_pot=solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,
	    #	minimag=None,maximag=None,distribution=distribution,npts_all=2*lmax+1,npts_zoom=npts_zoom,kini=12);
	    freqshift_spread=solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,
		distribution=distribution,kini=kini);
	    t2=ti.clock();
	    print "Ioct=",-Ioct,", time for solving [seconds]: ",t2-t1;
	    ind_spread=np.argsort(np.imag(freqshift_spread));
	    print freqshift_spread[ind_spread[:12]]/omega0;
	    #print freqshift_spread_pot/omega0

	    npts=100;factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
	    ReQminscan=np.array([-2.1,-1.1,-0.1,0.9,1.9])*Qs;
	    ReQmaxscan=np.array([-2.01,-1.01,-0.01,0.99,1.99])*Qs;
	    ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0.;

	    for iRe,ReQmin in enumerate(ReQminscan):
    		ReQmax=ReQmaxscan[iRe];
		real_range=np.linspace(ReQmin,ReQmax,npts+1);
    		imag_range=np.linspace(ImQmin,ImQmax,npts+1);
		ftable=np.zeros((npts+1,npts+1));
		for ir,r in enumerate(real_range):
		    for im,m in enumerate(imag_range):
			ftable[im,ir]=1./np.abs(determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution)/factnorm);

		fig1,ax1=init_figure();
		plot2D(ftable,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax1,colorlabel='');
		plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros (no Landau damping)','xk','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
		plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros (Landau damping)','ow','Im[Q]',ax1,0,xlab='Re[Q]',ms=10);
		#plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'pot. zeros (Landau damping)','ok','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
		ax1.set_xlim([ReQmin,ReQmax]);ax1.set_ylim([ImQmin,ImQmax]);
		end_figure(fig1,ax1,save=results_dir+'/plot2D_lmax'+str(lmax)+'_nmax'+str(nmax)+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1e11)+'e11_'+distribution+'_oct'+str(-Ioct)+'_l'+str(iRe-lmax)+'_kini'+str(kini));
	
	#pylab.show()


    if False:
	# 2D plot of matrix determinant
	ReQmin=-2.03*Qs;ReQmax=2*Qs;
	ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0;
	npts=200;
    	real_range=np.linspace(ReQmin,ReQmax,npts+1);
    	imag_range=np.linspace(ImQmin,ImQmax,npts+1);
	ftable=np.zeros((npts+1,npts+1));
	ftable_real=np.zeros((npts+1,npts+1));
	ftable_imag=np.zeros((npts+1,npts+1));
	factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian'));
	print factnorm

	for ir,r in enumerate(real_range):
	    for im,m in enumerate(imag_range):
		ftable_real[im,ir]=determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').real/factnorm;
		ftable_imag[im,ir]=determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian').imag/factnorm;
		ftable[im,ir]=1./np.abs(determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution='gaussian')/factnorm);
	
	# plot real part
	fig,ax=init_figure();
	plot2D(ftable_real,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='',colorlim=[-1,1]);
	plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'exact zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros from solve_lib','ow','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'potential zeros','ok','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	end_figure(fig,ax);
	
	# plot imag part
	fig,ax=init_figure();
	plot2D(ftable_imag,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='',colorlim=[-1,1]);
	plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros from solve_lib','ow','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'potential zeros','ok','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	end_figure(fig,ax);
	
	# plot 1/abs(f)
	fig,ax=init_figure();
	plot2D(ftable,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='');
	plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros from solve_lib','ow','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	#plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'potential zeros','ok','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	end_figure(fig,ax);
	pylab.show()
	
    if False:
	# 2D plot of matrix determinant, with Landau damping
	npts=100;
	#bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-30.,30.);
	print bx,bxy;
	factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
	print factnorm

	if False:
	    ReQmin=-1.06*Qs;ReQmax=-1.01*Qs;#ReQmin=-0.2e-4;ReQmax=-0.;
	    ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0.;
	    #ImQmin=-1e-6;ImQmax=-0.;
    	    real_range=np.linspace(ReQmin,ReQmax,npts+1);
    	    imag_range=np.linspace(ImQmin,ImQmax,npts+1);
	    ftable=np.zeros((npts+1,npts+1));
	    ftable_real=np.zeros((npts+1,npts+1));
	    ftable_imag=np.zeros((npts+1,npts+1));
	    for ir,r in enumerate(real_range):
		for im,m in enumerate(imag_range):
		    ftable_real[im,ir]=determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).real/factnorm;
		    ftable_imag[im,ir]=determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution).imag/factnorm;
		    ftable[im,ir]=1./np.abs(determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution)/factnorm);

    	    real_range2=np.linspace(freqshift[ind[0]].real/omega0+2*(2*bx+bxy),freqshift[ind[0]].real/omega0,npts_zoom+1);
    	    imag_range2=np.linspace(freqshift[ind[0]].imag/omega0,0.,npts_zoom+1);
	    gridR,gridI=np.meshgrid(real_range2,imag_range2);

	    fig,ax=init_figure();
	    plot2D(ftable_real,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='',colorlim=[-1,1]);
	    plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	    plot(gridR,gridI,'zoomed grid','+k','Im[Q]',ax,0,xlab='Re[Q]',ms=4);
	    ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	    end_figure(fig,ax);
	    fig,ax=init_figure();
	    plot2D(ftable_imag,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='',colorlim=[-1,1]);
	    plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	    plot(gridR,gridI,'zoomed grid','+k','Im[Q]',ax,0,xlab='Re[Q]',ms=4);
	    ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	    end_figure(fig,ax);
	    fig,ax=init_figure();
	    plot2D(ftable,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax,colorlabel='');
	    plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros (no Landau damping)','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	    plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros (Landau damping)','ow','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	    plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'pot. zeros (Landau damping)','ok','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	    plot(gridR,gridI,'zoomed grid','+k','Im[Q]',ax,0,xlab='Re[Q]',ms=4);
	    ax.set_xlim([ReQmin,ReQmax]);ax.set_ylim([ImQmin,ImQmax]);
	    end_figure(fig,ax);
	
	ReQminscan=np.array([-2.1,-1.1,-0.1,0.9,1.9])*Qs;
	ReQmaxscan=np.array([-2.01,-1.01,-0.01,0.99,1.99])*Qs;
	ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0.;
	
	for iRe,ReQmin in enumerate(ReQminscan):
    	    ReQmax=ReQmaxscan[iRe];
	    real_range=np.linspace(ReQmin,ReQmax,npts+1);
    	    imag_range=np.linspace(ImQmin,ImQmax,npts+1);
	    ftable=np.zeros((npts+1,npts+1));
	    for ir,r in enumerate(real_range):
		for im,m in enumerate(imag_range):
		    ftable[im,ir]=1./np.abs(determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution)/factnorm);
	
	    fig1,ax1=init_figure();
	    plot2D(ftable,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax1,colorlabel='');
	    plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros (no Landau damping)','xk','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
	    plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros (Landau damping)','ow','Im[Q]',ax1,0,xlab='Re[Q]',ms=10);
	    #plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'pot. zeros (Landau damping)','ok','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
	    ax1.set_xlim([ReQmin,ReQmax]);ax1.set_ylim([ImQmin,ImQmax]);
	    end_figure(fig1,ax1);

	# stability diagram plot
	Qscan_gau=10.**np.arange(-9,-2,0.05);
	Qscan_gau=np.hstack((-Qscan_gau[::-1],Qscan_gau));
	stabx_Python=np.zeros(len(Qscan_gau),dtype=complex);
	for iQ,Q in enumerate(Qscan_gau):
	    stabx_Python[iQ]=-1./dispersion_integral_oct_2D(Q,bx,bxy,distribution=distribution)
	    #print Q,stabx_Python[iQ];

	fig,ax=init_figure();
	plot(np.real(stabx_Python)+ReQmax,-np.imag(stabx_Python),'Stab. diagram','b','-Im[Q]',ax,0,xlab='Re[Q]');
	plot(np.real(freqshift)/omega0,-np.imag(freqshift)/omega0,'zeros','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
	ax.set_xlim([ReQmin,ReQmax]);
	end_figure(fig,ax);
	
	pylab.show()
	
    if True:
    
	Ioctscan=np.array([-100.,10.,100.]);
	Qpscan=np.arange(-10,21,2);
	dampscan=np.array([0,0.02]);#dampscan=np.array([0.]);
	nevery=1;lmax=2;nmax=4;M=1;nx=0;Nb=3e11;
	npts=100;col=['b','r','g','m','c','k','y'];pat=['x','o','+','^','v','d'];
	
	distribution='gaussian';kini=12;
    	# comparison with stability diagram
	avbetax=R/Qx;avbetay=R/Qy; # average beta functions used

	beam='1';
	param_filename_coll=root_LHC+'/Coll_settings/coll_ph1_beta_4000GeV_sq0p6_b'+beam+'_2012.txt';
	settings_filename_coll=root_LHC+'/Coll_settings/coll_settings_physics_fill_3265_B'+beam+'.txt';

	imp_mod,wake_mod=LHC_imp_model_v1(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	    beta_filename_coll=None,dire=root_LHC+"/LHC_elements/",commentcoll='_2012_v2',direcoll='Coll_2012_v2/',lxplusbatch='retrieve',
	    beam=beam,squeeze='0p6m_3m_0p6m_3m',wake_calc=False,ftypescan=0,nflog=100,zpar=z_param(),
	    flagplot=False,root_result=results_dir,commentsave='');

	# longitudinal distribution initialization
	g,a,b=longdistribution_decomp(taub,typelong="Gaussian");

	dnormfactor=compute_damper_matrix(0,0,0,M,0.,omega0,Qxfrac,a,b,taub,g);
	dnormfactor=2.*np.pi*dnormfactor[0,0,0,0];

	for damp in dampscan:
	
	    # select Zxdip
	    for iw in imp_mod:
		if test_impedance_wake_comp(iw,1,0,0,0,'x'): Z=deepcopy(iw.func[::nevery,:]);freq=deepcopy(iw.var[::nevery]);

	    Qdet=np.zeros((len(Ioctscan),len(Qpscan)),dtype=complex);
	    Qstab=np.zeros((len(Ioctscan),len(Qpscan)),dtype=complex);

	    for iQp,Qp in enumerate(Qpscan):

		t1=ti.clock()
		omegaksi=Qp*omega0/eta;
		coefdamper,coefZ=computes_coef(f0,damp,b,g[0],dnormfactor,taub,dphase,M,Nb,gamma,Qx,particle='proton');
		matdamper=compute_damper_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g);
		matZ=compute_impedance_matrix(lmax,nmax,nx,M,omegaksi,omega0,Qxfrac,a,b,taub,g,Z,freq,flag_trapz=1);
		t2=ti.clock();
		print "time for matrix computation [seconds]: ",t2-t1;

		t1=ti.clock()
		freqshift,v=eigenmodesDELPHI(lmax,nmax,matdamper,matZ,coefdamper,coefZ,omegas,flageigenvect=False);
		t2=ti.clock();
		print "time for solving eigenvalue problem [seconds]: ",t2-t1;

		ind=np.argsort(np.imag(freqshift));
		print freqshift[ind[:12]]/omega0;

		for iIoct,Ioct in enumerate(Ioctscan):
		    bx,bxy=detuning_coef_oct_LHC('x',gamma,2e-6,2e-6,-Ioct,Ioct);
		    print bx,bxy;

		    t1=ti.clock()
		    freqshift_spread=solve_determinantDELPHI_tunespread(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,
			distribution=distribution,kini=kini);
		    t2=ti.clock();
		    print "d=",damp,", Qp=",Qp,", Ioct=",-Ioct,", time for solving det. eq. [seconds]: ",t2-t1;

		    t1=ti.clock()
		    freqshift_stab=solve_stability_diagram(lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,
			distribution=distribution,kini=kini);
		    t2=ti.clock();
		    print "d=",damp,", Qp=",Qp,", Ioct=",-Ioct,", time for solving stab. diagram [seconds]: ",t2-t1;

		    if len(freqshift_spread)==0: freqshift_spread=np.array([0.]);
		    if len(freqshift_stab)==0: freqshift_stab=np.array([0.]);
		    
		    ind_spread=np.argsort(np.imag(freqshift_spread));
		    print freqshift_spread[ind_spread[:12]]/omega0;
		    ind_stab=np.argsort(np.imag(freqshift_stab));
		    print freqshift_stab[ind_stab[:12]]/omega0;

		    Qdet[iIoct,iQp]=freqshift_spread[ind_spread[0]]/omega0;
		    Qstab[iIoct,iQp]=freqshift_stab[ind_stab[0]]/omega0;

		    if True:
			# stability diagram plot
			ReQminscan=np.array([-2.1,-1.1,-0.1,0.9,1.9])*Qs;
			ReQmaxscan=np.array([-2.01,-1.01,-0.01,0.99,1.99])*Qs;
			ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0.;

			Qscan_gau=10.**np.arange(-9,-2,0.05);
			Qscan_gau=np.hstack((-Qscan_gau[::-1],Qscan_gau));
			stabx_Python=np.zeros(len(Qscan_gau),dtype=complex);

			for iRe,ReQmin in enumerate(ReQminscan):
    			    ReQmax=ReQmaxscan[iRe];
			    
			    for iQ,Q in enumerate(Qscan_gau):
				stabx_Python[iQ]=-1./dispersion_integral_oct_2D(Q,bx,bxy,distribution=distribution)
				#print Q,stabx_Python[iQ];

			    fig,ax=init_figure();
			    plot(np.real(stabx_Python)+ReQmax,-np.imag(stabx_Python),'Stab. diagram','b','-Im[Q]',ax,0,xlab='Re[Q]');
			    plot(np.real(freqshift)/omega0,-np.imag(freqshift)/omega0,'zeros (unperturbed)','xk','Im[Q]',ax,0,xlab='Re[Q]',ms=15);
			    ax.set_xlim([ReQmin,ReQmax]);
			    end_figure(fig,ax,save=results_dir+'/plot_stab_diagram_lmax'+str(lmax)+'_nmax'+str(nmax)+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1e11)+'e11_'+distribution+'_oct'+str(-Ioct)+'_l'+str(iRe-lmax));
		    
		    if True:
			# 2D plots
			factnorm=np.abs(determinantDELPHI_tunespread(0,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution));
			ReQminscan=np.array([-2.1,-1.1,-0.1,0.9,1.9])*Qs;
			ReQmaxscan=np.array([-2.01,-1.01,-0.01,0.99,1.99])*Qs;
			ImQmin=1.1*freqshift[ind[0]].imag/omega0;ImQmax=0.;

			for iRe,ReQmin in enumerate(ReQminscan):
    			    ReQmax=ReQmaxscan[iRe];
			    real_range=np.linspace(ReQmin,ReQmax,npts+1);
    			    imag_range=np.linspace(ImQmin,ImQmax,npts+1);
			    ftable=np.zeros((npts+1,npts+1));
			    for ir,r in enumerate(real_range):
				for im,m in enumerate(imag_range):
				    ftable[im,ir]=1./np.abs(determinantDELPHI_tunespread(r+1j*m,lmax,nmax,matdamper,matZ,coefdamper,coefZ,bx,bxy,omega0,omegas,distribution=distribution)/factnorm);

			    fig1,ax1=init_figure();
			    plot2D(ftable,ReQmin,ReQmax,ImQmin,ImQmax,'Re[Q]','Im[Q]','',ax1,colorlabel='');
			    plot(np.real(freqshift)/omega0,np.imag(freqshift)/omega0,'zeros (no Landau damping)','xk','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
			    plot(np.real(freqshift_spread)/omega0,np.imag(freqshift_spread)/omega0,'zeros (Landau damping)','ow','Im[Q]',ax1,0,xlab='Re[Q]',ms=10);
			    #plot(freqshift_spread_pot[:,0]/omega0,freqshift_spread_pot[:,1]/omega0,'pot. zeros (Landau damping)','ok','Im[Q]',ax1,0,xlab='Re[Q]',ms=15);
			    ax1.set_xlim([ReQmin,ReQmax]);ax1.set_ylim([ImQmin,ImQmax]);
			    end_figure(fig1,ax1,save=results_dir+'/plot2D_lmax'+str(lmax)+'_nmax'+str(nmax)+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1e11)+'e11_'+distribution+'_oct'+str(-Ioct)+'_l'+str(iRe-lmax)+'_kini'+str(kini));
			    os.system("rm -f "+results_dir+"/plot2D_lmax"+str(lmax)+'_nmax'+str(nmax)+'_d'+float_to_str(damp)+'_Qp'+float_to_str(Qp)+'_Nb'+float_to_str(Nb/1e11)+'e11_'+distribution+'_oct'+str(-Ioct)+'_l'+str(iRe-lmax)+'_kini'+str(kini)+'.eps');

	    fig,ax=init_figure();
	    for iIoct,Ioct in enumerate(Ioctscan):
		plot(Qpscan,-np.imag(Qdet[iIoct,:]),'DELPHI, Ioct='+str(-Ioct)+' A','-'+col[iIoct],'-Im[Q]',ax,0,xlab=" $ Q' $ ");
		plot(Qpscan,-np.imag(Qstab[iIoct,:]),'Stab. diagram, Ioct='+str(-Ioct)+' A',pat[iIoct]+col[iIoct],'-Im[Q]',ax,0,xlab=" $ Q' $ ",ms=15);

	    end_figure(fig,ax,save=results_dir+'/plot_vs_Qp_lmax'+str(lmax)+'_nmax'+str(nmax)+'_d'+float_to_str(damp)+'_Nb'+float_to_str(Nb/1e11)+'e11_'+distribution+'_kini'+str(kini));    
	    
	#pylab.show()
