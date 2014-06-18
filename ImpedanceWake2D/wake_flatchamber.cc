/*
 *  wake_flatchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the wakes in a flat chamber (see CERN note by N. Mounet and E. Metral, 
 "Electromagnetic fields and beam coupling impedances in a multilayer flat chamber", 2010)
 
In input : typical input file is

Machine:	LHC
Relativistic Gamma:	479.6
Impedance Length in m:	1
Number of upper layers in the chamber wall:	2
Layer 1 inner half gap in mm:	18.375
Layer 1 DC resistivity (Ohm.m):	2e-10
Layer 1 relaxation time for resistivity (ps):	2.1
Layer 1 real part of dielectric constant:	1
Layer 1 magnetic susceptibility:	0
Layer 1 relaxation frequency of permeability (MHz):	Infinity
Layer 1 thickness in mm:	0.05
Layer 2 DC resistivity (Ohm.m):	7.2e-7
Layer 2 relaxation time for resistivity (ps):	0
Layer 2 real part of dielectric constant:	1.43
Layer 2 magnetic susceptibility:	0.02
Layer 2 relaxation frequency of permeability (MHz):	Infinity
Layer 2 thickness in mm:	Infinity
Top bottom symmetry (yes or no):	yes
Number of lower layers in the chamber wall:	2
Layer -1 inner half gap in mm:	18.375
Layer -1 DC resistivity (Ohm.m):	2e-10
Layer -1 relaxation time for resistivity (ps):	2.1
Layer -1 real part of dielectric constant:	1
Layer -1 magnetic susceptibility:	0
Layer -1 relaxation frequency of permeability (MHz):	Infinity
Layer -1 thickness in mm:	0.05
Layer -2 DC resistivity (Ohm.m):	7.2e-7
Layer -2 relaxation time for resistivity (ps):	0
Layer -2 real part of dielectric constant:	1.43
Layer -2 magnetic susceptibility:	0.02
Layer -2 relaxation frequency of permeability (MHz):	Infinity
Layer -2 thickness in mm:	Infinity
linear (1) or logarithmic (0) or both (2) scan in z for the wake:	0
sampling distance in m for the linear sampling:	0.5e-5
zmin in m of the linear sampling:	0.5e-5
zmax in m of the linear sampling:	0.01
Number of points per decade for the logarithmic sampling:	100
exponent (10^) of zmin (in m) of the logarithmic sampling:	-1.99
exponent (10^) of zmax (in m) of the logarithmic sampling:	6
added z [m]:	5e6 1e7
factor weighting the longitudinal impedance error:	100.
tolerance (in wake units) to achieve:	1.e9
frequency above which the mesh bisecting is linear [Hz]:	1.e11
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If top-bottom symmetry is set (with "yes" or "y" or "1") the lower layers (with a
minus sign) are ignored. Also if there are more layers than indicated by the number of upper (lower) layers,
the additional one are ignored. The last layer is always assumed to go to infinity.

In output one gives six files with the wakes (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant term). Each has 2 columns : distance (behind source) and
wake (SI units).
In output one also gives six files with the impedances computed on the final mesh used to 
compute the wakes (files are for longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant, impedances). Each file has 3 columns : frequency, real part
and imaginary part of the impedance (SI units).

 */

#include <iostream>
#include <sstream>
#include <vector>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sys/time.h>

#include <complex>
#include <cmath>
#include <stdlib.h>
#include <ablas.h>
//#include <matinv.h>
//#include <densesolver.h>
//#include <trfac.h>
#include <amp.h>
#include <mpfr.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort_double.h>

/* to use cmatrixgemm (ALGLIB routine -- 16.12.2009 Bochkanov Sergey)
INPUT PARAMETERS
  M - matrix size, M>0
  N - matrix size, N>0
  K - matrix size, K>0
  Alpha - coefficient 
  A - matrix 
  IA - submatrix offset 
  JA - submatrix offset 
  OpTypeA - transformation type: * 0 - no transformation * 1 - transposition * 2 - conjugate transposition 
  B - matrix 
  IB - submatrix offset 
  JB - submatrix offset 
  OpTypeB - transformation type: * 0 - no transformation * 1 - transposition * 2 - conjugate transposition 
  Beta - coefficient 
  C - matrix 
  IC - submatrix offset 
  JC - submatrix offset 

template<unsigned int Precision> void cmatrixgemm(int m, int n, int k, 
	amp::campf<Precision> alpha, const ap::template_2d_array< amp::campf<Precision> >& a,
	int ia, int ja, int optypea, const ap::template_2d_array< amp::campf<Precision> >& b, 
	int ib, int jb, int optypeb, amp::campf<Precision> beta, 
	ap::template_2d_array< amp::campf<Precision> >& c, int ic, int jc);
*/


#define  Precision	200  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension
#define  MAXMEM		150000 // Maximum number of elements of arrays with impedances, etas and chis

using std::complex;
using std::log;
using std::exp;
using std::abs;
using std::max;
using std::min;
using std::cout;

const  amp::ampf<Precision> C=299792458;    // velocity of light [m/s]

struct params {unsigned int m; unsigned int n;
  	unsigned int N; unsigned int M;
	ap::template_1d_array< amp::campf<Precision> > eps1;
  	ap::template_1d_array< amp::campf<Precision> > mu1;
  	ap::template_1d_array< amp::campf<Precision> > nu2;
	ap::template_1d_array< amp::campf<Precision> > eps1ratio;
  	ap::template_1d_array< amp::campf<Precision> > mu1ratio;
  	ap::template_1d_array< amp::campf<Precision> > nu2ratio;
	ap::template_1d_array< amp::ampf<Precision> > b;
	ap::template_1d_array< amp::campf<Precision> > eps1m;
  	ap::template_1d_array< amp::campf<Precision> > mu1m;
  	ap::template_1d_array< amp::campf<Precision> > nu2m;
	ap::template_1d_array< amp::campf<Precision> > eps1mratio;
  	ap::template_1d_array< amp::campf<Precision> > mu1mratio;
  	ap::template_1d_array< amp::campf<Precision> > nu2mratio;
	ap::template_1d_array< amp::ampf<Precision> > bm;
	amp::ampf<Precision> beta;
	amp::ampf<Precision> kovergamma;};
// eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.


struct params_diff {unsigned int N; unsigned int M;
        ap::template_1d_array< amp::ampf<Precision> > rho;
        ap::template_1d_array< amp::ampf<Precision> > tau;
        ap::template_1d_array< amp::ampf<Precision> > epsb;
        ap::template_1d_array< amp::ampf<Precision> > chi;
        ap::template_1d_array< amp::ampf<Precision> > fmu;
        ap::template_1d_array< amp::ampf<Precision> > b;
        ap::template_1d_array< amp::ampf<Precision> > rhom;
        ap::template_1d_array< amp::ampf<Precision> > taum;
        ap::template_1d_array< amp::ampf<Precision> > epsbm;
        ap::template_1d_array< amp::ampf<Precision> > chim;
        ap::template_1d_array< amp::ampf<Precision> > fmum;
        ap::template_1d_array< amp::ampf<Precision> > bm;
        amp::ampf<Precision> beta; amp::ampf<Precision> gamma;
	int flag_topbotsym; double L; double *freqi; unsigned int interp_type; unsigned long nf;
	complex<double> *Zxdipfi; complex<double> *Zxdipdi;
	complex<double> *Zydipfi; complex<double> *Zydipdi;
	complex<double> *Zyquadfi; complex<double> *Zyquaddi;
        complex<double> *Zlongfi; complex<double> *Zlongdi;
	complex<double> *Zycstfi; complex<double> *Zycstdi;
	complex<double> *x;};


// global arrays with memory of etas and chis
ap::template_1d_array< amp::ampf<Precision> > kxmem;
ap::template_1d_array< amp::campf<Precision> > eta1mem,eta2mem,chi1mem,chi2mem;
unsigned long mem; // current number of elements of kxmem, eta1mem, etc.

// global arrays with memory of impedances
double *freqmem;
complex<double> *Zxdipmem,*Zydipmem,*Zlongmem,*Zyquadmem,*Zycstmem;
complex<double> *tolintabs00mem, *tolintabs01mem, *tolintabs02mem, *tolintabs11mem;
unsigned long impmem; // current number of elements of freqmem, Zxdipmem, etc.

double factlong=100.; // multiplication factor for longitudinal impedance (precision should be better for long. wake) (default value)

  
/******************************************************************************
 *** locateMP: search a multiprecision table ordered in ascending order	    ***
 *** Effect         : Function that gives the position lprov (integer) in   ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
 *** Parameters     : table, z, n (table is indexed from 0 to n)            ***
 ******************************************************************************/

unsigned long locateMP (ap::template_1d_array< amp::ampf<Precision> >& table, 
	amp::ampf<Precision> z, unsigned long n)

{ unsigned long il, iu, im, lprov;

 il=0;
 iu=n+1;
 while (iu-il >1) {
   im=(iu+il)/2; // midpoint
   if (z >= table(im))
     il=im;
   else
     iu=im;
   }
 if (z==table(0)) lprov=0;
 else if (z==table(n)) lprov=n; 
 else if (z<table(0)) lprov=0; 
 else if (z>table(n)) lprov=n+1; 
 else lprov=iu;

 return lprov;

}

/******************************************************************************
 *** locate: search a table of doubles ordered in ascending order	    ***
 *** Effect         : Function that gives the position lprov (integer) in   ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
 *** Parameters     : table, z, n (table is indexed from 0 to n)            ***
 ******************************************************************************/

unsigned long locate (double *table, double z, unsigned long n)

{ unsigned long il, iu, im, lprov;

 il=0;
 iu=n+1;
 while (iu-il >1) {
   im=(iu+il)/2; // midpoint
   if (z >= table[im])
     il=im;
   else
     iu=im;
   }
 if (z==table[0]) lprov=0;
 else if (z==table[n]) lprov=n; 
 else if (z<table[0]) lprov=0; 
 else if (z>table[n]) lprov=n+1; 
 else lprov=iu;

 return lprov;

}

/**************************************************
*** csqrt: Complex square root in multiprecision
***
***************************************************/

  amp::campf<Precision> csqrt(amp::campf<Precision> z){
  
    /* Complex square root of z in multiprecision. Convention specified in CERN note mentioned at the
    beginning. */
  
    amp::campf<Precision> result;
    amp::ampf<Precision> rho,phi;
    
    rho=amp::sqrt(amp::abscomplex(z));
    phi=amp::atan2(z.y,z.x)/2;
    
    result.x=rho*amp::cos(phi);
    result.y=rho*amp::sin(phi);
    
    return result;
  
  }
  

/**************************************************
*** cexp: Complex exponential in multiprecision
***
***************************************************/

  amp::campf<Precision> cexp(amp::campf<Precision> z){
  
    /* Complex exponential of z in multiprecision. */
  
    amp::campf<Precision> result;
    amp::ampf<Precision> rho,phi;
    
    rho=amp::exp(z.x);
    phi=z.y;
    
    result.x=rho*amp::cos(phi);
    result.y=rho*amp::sin(phi);
    
    return result;
  
  }
  
  
/**************************************************
*** multilayerm: computes matrix for field matching (multiprecision)
***
***************************************************/

  void multilayerm(ap::template_2d_array< amp::campf<Precision> >& m,
  	amp::campf<Precision>& kyn,
  	unsigned int N,
	ap::template_1d_array< amp::campf<Precision> >& eps1, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2, 
	ap::template_1d_array< amp::campf<Precision> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	amp::ampf<Precision> kx, amp::ampf<Precision> beta){
	
    /* computes the matrix M for the multilayer field matching, from mu1, eps1, b of each of the N layers
    and the horizontal wave number kx, angular frequency omega, relativistic velocity factor beta
    and wave number k.
    It also gives in output kyn, which is ky for the last layer. */

    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.
    
    amp::campf<Precision> kyp,kyp1,exp1,exp2,exp3,exp4;
    amp::campf<Precision> tmp1,tmp2,tmpplus,tmpminus,fac,fac2;
    amp::ampf<Precision> kx2,fac3;
    ap::template_2d_array< amp::campf<Precision> > mp1p;
    ap::template_2d_array< amp::campf<Precision> > mold;
   
    mp1p.setbounds(1,4,1,4);
    mold.setbounds(1,4,1,4);
  
    for (int i=1; i<=4; i++) {
      for (int j=1; j<=4; j++) {
        if (i==j) mold(i,j)=1;
	else mold(i,j)=0;
      }
    }
    
    kx2=amp::sqr(kx);
    kyp=csqrt(kx2+nu2(1));
    fac3=kx/(2*beta);

    m=mold;
    
    for (unsigned int p=1; p<=N-1; p++) {

      kyp1=csqrt(kx2+nu2(p+1));
      
      tmp1=kyp*nu2ratio(p)/kyp1;
      exp1=cexp((kyp-kyp1)*b(p));
      exp2=cexp(-(kyp+kyp1)*b(p));
      exp3=1/exp2;
      exp4=1/exp1;
      
      tmp2=tmp1*eps1ratio(p);
      tmpplus=(1+tmp2)/2;
      tmpminus=(1-tmp2)/2;
      mp1p(1,1)=tmpplus*exp1;
      mp1p(1,2)=tmpminus*exp2;
      mp1p(2,1)=tmpminus*exp3;
      mp1p(2,2)=tmpplus*exp4;
     
      fac=fac3*(nu2ratio(p)-1)/kyp1;
      fac2=fac*eps1(p+1);
      mp1p(1,3)=-fac2*exp1;
      mp1p(1,4)=-fac2*exp2;
      mp1p(2,3)=fac2*exp3;
      mp1p(2,4)=fac2*exp4;
      
      tmp2=tmp1*mu1ratio(p);
      tmpplus=(1+tmp2)/2;
      tmpminus=(1-tmp2)/2;
      mp1p(3,3)=tmpplus*exp1;
      mp1p(3,4)=tmpminus*exp2;
      mp1p(4,3)=tmpminus*exp3;
      mp1p(4,4)=tmpplus*exp4;
     
      fac2=fac*mu1(p+1);
      mp1p(3,1)=-fac2*exp1;
      mp1p(3,2)=-fac2*exp2;
      mp1p(4,1)=fac2*exp3;
      mp1p(4,2)=fac2*exp4;
      
      ablas::cmatrixgemm<Precision>(4,4,4,1,mp1p,1,1,0,mold,1,1,0,0,m,1,1);
      
      kyp=kyp1;
      mold=m;
    }
    
    kyn=kyp1;
  
    return;
  }
  
  
/**************************************************
*** matinv4: computes the 2 first lines of coefficients
*** of the inverse of a 4x4 matrix (multiprecision).
*** This is optimized.
***************************************************/

  void matinv4(ap::template_2d_array< amp::campf<Precision> >& mat) {
  
    /* input= matrix mat (bounds 0..3,0..3). output is also mat, with the first two lines
       filled with the coefficients of its inverse */
    
    // We need the two first columns of the cofactor matrix, and the total determinant
    ap::template_2d_array< amp::campf<Precision> > cof; // cofactor matrix (two first columns only)
    ap::template_1d_array< amp::campf<Precision> > det2; // array with the 2x2 determinants of the 2 last columns
    amp::campf<Precision> det4; //determinant of the initial 4x4 matrix mat
    
    cof.setbounds(0,3,0,1);
    det2.setbounds(0,5);
    
    // det2 is ordered from top to bottom
    det2(0)=mat(0,2)*mat(1,3)-mat(0,3)*mat(1,2);
    det2(1)=mat(0,2)*mat(2,3)-mat(0,3)*mat(2,2);
    det2(2)=mat(0,2)*mat(3,3)-mat(0,3)*mat(3,2);
    det2(3)=mat(1,2)*mat(2,3)-mat(1,3)*mat(2,2);
    det2(4)=mat(1,2)*mat(3,3)-mat(1,3)*mat(3,2);
    det2(5)=mat(2,2)*mat(3,3)-mat(2,3)*mat(3,2);
    
    // cofactors of the first column
    cof(0,0)=  mat(1,1)*det2(5)-mat(2,1)*det2(4)+mat(3,1)*det2(3);
    cof(1,0)=-(mat(0,1)*det2(5)-mat(2,1)*det2(2)+mat(3,1)*det2(1));
    cof(2,0)=  mat(0,1)*det2(4)-mat(1,1)*det2(2)+mat(3,1)*det2(0);
    cof(3,0)=-(mat(0,1)*det2(3)-mat(1,1)*det2(1)+mat(2,1)*det2(0));
    
    // total determinant
    det4=mat(0,0)*cof(0,0)+mat(1,0)*cof(1,0)+mat(2,0)*cof(2,0)+mat(3,0)*cof(3,0);
    
    // cofactors of the second column
    cof(0,1)=-(mat(1,0)*det2(5)-mat(2,0)*det2(4)+mat(3,0)*det2(3));
    cof(1,1)=  mat(0,0)*det2(5)-mat(2,0)*det2(2)+mat(3,0)*det2(1);
    cof(2,1)=-(mat(0,0)*det2(4)-mat(1,0)*det2(2)+mat(3,0)*det2(0));
    cof(3,1)=  mat(0,0)*det2(3)-mat(1,0)*det2(1)+mat(2,0)*det2(0);
    
    // final coefficients sought for (we transpose the cofactors and divide by determinant)
    mat(0,0)=cof(0,0)/det4;
    mat(0,1)=cof(1,0)/det4;
    mat(0,2)=cof(2,0)/det4;
    mat(0,3)=cof(3,0)/det4;
    mat(1,0)=cof(0,1)/det4;
    mat(1,1)=cof(1,1)/det4;
    mat(1,2)=cof(2,1)/det4;
    mat(1,3)=cof(3,1)/det4;
    
    
  }
  
/**************************************************
*** etachi: computes eta1, eta2, chi1 and chi2 (multiprecision)
***
***************************************************/

  void etachi(amp::campf<Precision>& chi1,
  	amp::campf<Precision>& chi2,
  	amp::campf<Precision>& eta1,
  	amp::campf<Precision>& eta2,
  	unsigned int N, unsigned int M,
	ap::template_1d_array< amp::campf<Precision> >& eps1, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2, 
	ap::template_1d_array< amp::campf<Precision> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	ap::template_1d_array< amp::campf<Precision> >& eps1m, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1m, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2m, 
	ap::template_1d_array< amp::campf<Precision> >& eps1mratio, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1mratio, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2mratio, 
	ap::template_1d_array< amp::ampf<Precision> >& bm,
	amp::ampf<Precision> kx, amp::ampf<Precision> beta){
	
    /* function that computes chi1, chi2, eta1 and eta2 for a given horizontal wave number kx, from mu1, eps1, b 
    of each of the N upper layers and mu1m, eps1m, bm of each of the M lower layers, 
    and from the angular frequency omega, relativistic velocity factor beta and wave number k. */
    
    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.

    ap::template_2d_array< amp::campf<Precision> > mat,matprime; // 4*4 field matching matrices for upper and lower layers
    ap::template_2d_array< amp::campf<Precision> > pcal; // 4*4 final matrix to be inverted, and then its inverse
    amp::campf<Precision> kyn,kym; // ky for the last layers
    int i,l; // flags for last layers
    int info;
    unsigned int lprov;
    //matinv::matinvreport<Precision> repi;
    timeval c1,c2; // clock ticks


    // try to find kx in the table kxmem
    if (mem==0) lprov=0;
    else lprov=locateMP(kxmem,kx,mem-1);
 
    if ( (mem!=0)&&(kx==kxmem(lprov)) ) {
      eta1=eta1mem(lprov);eta2=eta2mem(lprov);
      chi1=chi1mem(lprov);chi2=chi2mem(lprov);
    }
    else if ( ( (mem!=0)&&(lprov>0)) && (kx==kxmem(lprov-1)) ) {
      eta1=eta1mem(lprov-1);eta2=eta2mem(lprov-1);
      chi1=chi1mem(lprov-1);chi2=chi2mem(lprov-1);
    }
    else {

      /* setting bounds for matrices (beware, pcal first index has to be set to 0 instead of 1, because of the
      function matinv::cmatrixinverse */
      mat.setbounds(1,4,1,4);
      matprime.setbounds(1,4,1,4);
      pcal.setbounds(0,3,0,3);

      // compute the field matching 4*4 matrices (upper and lower layers)
      //gettimeofday(&c1,0);
      multilayerm(mat, kyn, N, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, kx, beta);
      multilayerm(matprime, kym, M, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, kx, beta);
      //gettimeofday(&c2,0);
      //std::cout << "multilayer: time= " << (c2.tv_sec-c1.tv_sec)*1.e6+(c2.tv_usec-c1.tv_usec) << " usec\n";

      i=1;l=2;

      // compute the final 4*4 matrix P and invert it
      for (int j=0; j<=3; j++) {
	pcal(0,j)=mat(i,j+1);
	pcal(1,j)=mat(i+2,j+1);
	pcal(2,j)=matprime(l,j+1);
	pcal(3,j)=matprime(l+2,j+1);
      }

      /*for (int k=1; k<=4; k++) {
	for (int j=1; j<=4; j++) {
          printf("%d %d : %s %s\n", k,j,mat(k,j).x.toDec().c_str(),mat(k,j).y.toDec().c_str());
	}
	printf("\n");
      }*/

      //matinv::cmatrixinverse(pcal,4,info,repi);
      matinv4(pcal); // ~50 times quicker !
      /*amp::ampf<Precision> sum=0;
      for (unsigned int p=0;p<=1;p++) {
        for (unsigned int q=0;q<=3;q++) sum+=amp::abscomplex((pcal(p,q)-pcal2(p,q))/pcal(p,q));
      }
      std::cout << "difference: " << sum.toDec().c_str() << "\n";*/

      /*for (int k=0; k<=3; k++) {
	for (int j=0; j<=3; j++) {
          printf("%d %d : %s %s\n", k+1,j+1,pcal(k,j).x.toDec().c_str(),pcal(k,j).y.toDec().c_str());
	}
	printf("\n");
      }*/

      // computes chi1, chi2, eta1 and eta2 at kx
      // 4 next lines: matinv version
      chi1=pcal(0,0)*mat(i,2)+pcal(0,1)*mat(i+2,2);
      chi2=pcal(1,0)*mat(i,2)+pcal(1,1)*mat(i+2,2);
      eta1=pcal(0,2)*matprime(l,1)+pcal(0,3)*matprime(l+2,1);
      eta2=pcal(1,2)*matprime(l,1)+pcal(1,3)*matprime(l+2,1);
      //chi1=sol(0,0);chi2=sol(1,0);eta1=sol(0,1);eta2=sol(1,1);

      /*printf("kx : %s\n", kx.toDec().c_str());
      printf("chi1 : %s %s\n", chi1.x.toDec().c_str(),chi1.y.toDec().c_str());
      printf("chi2 : %s %s\n", chi2.x.toDec().c_str(),chi2.y.toDec().c_str());
      printf("eta1 : %s %s\n", eta1.x.toDec().c_str(),eta1.y.toDec().c_str());
      printf("eta2 : %s %s\n", eta2.x.toDec().c_str(),eta2.y.toDec().c_str());*/
     
     for(unsigned int k=mem; k>=lprov+1; k--) {
       kxmem(k)=kxmem(k-1);
       eta1mem(k)=eta1mem(k-1);eta2mem(k)=eta2mem(k-1);
       chi1mem(k)=chi1mem(k-1);chi2mem(k)=chi2mem(k-1);
     }
     kxmem(lprov)=kx;
     eta1mem(lprov)=eta1;eta2mem(lprov)=eta2;
     chi1mem(lprov)=chi1;chi2mem(lprov)=chi2;
     mem++;
       
    }

    return;
    
  }
  
/**************************************************
*** integrand: computes the complex integrand in u (std precision)
***
***************************************************/

  std::complex<double> integrand(unsigned int m, unsigned int n,
  	unsigned int N, unsigned int M,
	ap::template_1d_array< amp::campf<Precision> >& eps1, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2, 
	ap::template_1d_array< amp::campf<Precision> >& eps1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1ratio, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2ratio, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	ap::template_1d_array< amp::campf<Precision> >& eps1m, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1m, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2m, 
	ap::template_1d_array< amp::campf<Precision> >& eps1mratio, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1mratio, 
  	ap::template_1d_array< amp::campf<Precision> >& nu2mratio, 
	ap::template_1d_array< amp::ampf<Precision> >& bm,
	amp::ampf<Precision> u, amp::ampf<Precision> beta,
	amp::ampf<Precision> kovergamma){
	
    /* function that computes the integrand in alphamn for a given u (kx=k sinh(u)/gamma) 
    and given azimuthal mode numbers m and n, from mu1, eps1, b 
    of each of the N upper layers and mu1m, eps1m, bm of each of the M lower layers, 
    and from the angular frequency omega, relativistic velocity factor beta, wave number k and k/gamma. */
    
    // NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.
    
    amp::campf<Precision> chi1,chi2,eta1,eta2; // chi and eta functions of kx
    amp::ampf<Precision> kx;
    amp::campf<Precision> inte; // result in multiprecision
    int m1powm,m1pown,m1powmn;
    double x,y;
    std::complex<double> result; // result in standard double precision

    if (m % 2 ==0) m1powm=1;
    else m1powm=-1;
    if (n % 2 ==0) m1pown=1;
    else m1pown=-1;
    
    m1powmn=m1powm*m1pown;
    
    // computes kx
    kx=kovergamma*amp::sinh(u);
    
    // computes eta and chi functions at kx
    etachi(chi1, chi2, eta1, eta2, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, kx, beta);
    
    // computes integrand
    inte=amp::cosh(m*u)*amp::cosh(n*u)*( chi1+m1powm*eta1+m1pown*chi2+m1powmn*eta2 );
    
    // check if real or imaginary part of inte is NaN. In that case, replace it by zero.
    x=double(amp::ampf<Precision>(inte.x).toDouble());
    y=double(amp::ampf<Precision>(inte.y).toDouble());
    if (x != x) x=0; // this is the way to check if it's a NaN (all comparisons give false with a NaN)
    if (y != y) y=0; // this is the way to check if it's a NaN (all comparisons give false with a NaN)

    result=std::complex<double>(x,y);
    //printf("%13.8e %13.8e %13.8e\n",double(amp::ampf<Precision>(u).toDouble()),x,y);
    return result;
    
  }
  
/**************************************************
*** integrand_real: computes the real part of the integrand in u (std precision)
***
***************************************************/

  double integrand_real(double x, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<Precision> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<Precision> > b,bm;
    amp::ampf<Precision> beta,gamma,kovergamma,u;
    std::complex<double> inte;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    
    u=x;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma);

    return inte.real();
    
  }

/**************************************************
*** integrand_imag: computes the imag. part of the integrand in u (std precision)
***
***************************************************/

  double integrand_imag(double x, void *p){
	
    /* function encapsulating the computation of the imaginary part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<Precision> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<Precision> > b,bm;
    amp::ampf<Precision> beta,gamma,kovergamma,u;
    std::complex<double> inte;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    
    
    u=x;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma);

    return inte.imag();
    
  }


/**************************************************
*** integrand_real_modif: computes the real part of the integrand in t (std precision)
*** instead of u (change of variable u=(1-t)/t )
***************************************************/

  double integrand_real_modif(double t, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<Precision> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<Precision> > b,bm;
    amp::ampf<Precision> beta,gamma,kovergamma,u;
    std::complex<double> inte;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    
    
    if (t==0) u=1.e4; // integrand should be very close to zero anyway
    else u=(1.-t)/t;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma);

    if ((t==0)&&(std::abs(inte)>1.e-10)) printf ("Warning: approx for t=0 too rough, %13.8e\n", std::abs(inte));

    return inte.real()/(t*t);
    
  }


/**************************************************
*** integrand_imag_modif: computes the imag part of the integrand in t (std precision)
*** instead of u (change of variable u=(1-t)/t )
***************************************************/

  double integrand_imag_modif(double t, void *p){
	
    /* function encapsulating the computation of the real part of the integrand of alphamn, 
    for gsl integration, see also function integrand */
    
    struct params *param=(struct params *)p;
    unsigned int N,M; // number of upper and lower layers
    unsigned int m,n; // indices of alphamn (azimuthal mode numbers)
    ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,eps1ratio,eps1mratio,mu1,mu1m,mu1ratio,mu1mratio;
    ap::template_1d_array< amp::campf<Precision> > nu2,nu2m,nu2ratio,nu2mratio;
    ap::template_1d_array< amp::ampf<Precision> > b,bm;
    amp::ampf<Precision> beta,gamma,kovergamma,u;
    std::complex<double> inte;

    m=(param->m);
    n=(param->n);
    M=(param->M);
    N=(param->N);
    eps1=(param->eps1);
    eps1m=(param->eps1m);
    eps1ratio=(param->eps1ratio);
    eps1mratio=(param->eps1mratio);
    mu1=(param->mu1);
    mu1m=(param->mu1m);
    mu1ratio=(param->mu1ratio);
    mu1mratio=(param->mu1mratio);
    nu2=(param->nu2);
    nu2m=(param->nu2m);
    nu2ratio=(param->nu2ratio);
    nu2mratio=(param->nu2mratio);
    b=(param->b);
    bm=(param->bm);
    beta=(param->beta);
    kovergamma=(param->kovergamma);
    
    
    if (t==0) u=1.e4; // integrand should be very close to zero anyway
    else u=(1.-t)/t;

    inte=integrand(m, n, N, M, eps1, mu1, nu2, eps1ratio, mu1ratio, nu2ratio, b, eps1m, mu1m, nu2m, eps1mratio, mu1mratio, nu2mratio, bm, u, beta, kovergamma);
    
    if ((t==0)&&(std::abs(inte)>1.e-10)) printf ("Warning: approx for t=0 too rough, %13.8e\n", std::abs(inte));

    return inte.imag()/(t*t);
    
  }


/**************************************************
*** integrate: performs the integration (using GSL)
*** of integrand_real or integrand_imag
***************************************************/

  double integrate(int flagreal, unsigned int M, unsigned int N, 
  	ap::template_1d_array< amp::ampf<Precision> > b, 
	ap::template_1d_array< amp::ampf<Precision> > bm, 
	amp::ampf<Precision> beta, ap::template_1d_array< amp::campf<Precision> > eps1,
	ap::template_1d_array< amp::campf<Precision> > eps1m,
	ap::template_1d_array< amp::campf<Precision> > mu1,
	ap::template_1d_array< amp::campf<Precision> > mu1m,
	amp::ampf<Precision> omega, amp::ampf<Precision> k, amp::ampf<Precision> kovergamma,
	unsigned int m, unsigned int n, double tolintabs, size_t limit, gsl_integration_workspace *w){
	
    /* In input: flagreal : 1 to integrate integrand_real, otherwise integrate integrand_imag
       The final m and n are the azimuthal mode numbers (e.g. m=0, n=0 will give alpha_00 ), 
       and tolintabs is the absolute error permitted
       The rest are the parameters of the multilayer computation (see etachi and integrand
       functions)*/

    
    struct params param; // input parameters for the integrand functions
    ap::template_1d_array< amp::campf<Precision> > eps1ratio,eps1mratio,mu1ratio,mu1mratio,nu2,nu2m,nu2ratio,nu2mratio;
    amp::ampf<Precision> beta2,k2;
    gsl_function F;
    double tolint=1.e-6; // relative error permitted for gsl adaptative integration
    double x,err; // result and error
    int status;
    
    if ((m==0)&&(n==0)) tolint*=1e-2; // higher precision for alpha00
 
    // precompute various ratio
    eps1ratio.setbounds(1,N);mu1ratio.setbounds(1,N);nu2.setbounds(1,N+1);nu2ratio.setbounds(1,N);
    eps1mratio.setbounds(1,M);mu1mratio.setbounds(1,M);nu2m.setbounds(1,M+1);nu2mratio.setbounds(1,M);
    beta2=amp::sqr(beta);
    k2=amp::sqr(k);
    //upper layers
    nu2(1)=k2*(1-beta2*eps1(1)*mu1(1));
    for (unsigned int p=1; p<=N; p++) {
      nu2(p+1)=k2*(1-beta2*eps1(p+1)*mu1(p+1));
      eps1ratio(p)=eps1(p)/eps1(p+1);
      mu1ratio(p)=mu1(p)/mu1(p+1);
      nu2ratio(p)=nu2(p+1)/nu2(p);
      eps1(p)=1/eps1(p);
      mu1(p)=1/mu1(p);
    }
    eps1(N+1)=1/eps1(N+1);
    mu1(N+1)=1/mu1(N+1);
    //lower layers
    nu2m(1)=k2*(1-beta2*eps1m(1)*mu1m(1));
    for (unsigned int p=1; p<=M; p++) {
      nu2m(p+1)=k2*(1-beta2*eps1m(p+1)*mu1m(p+1));
      eps1mratio(p)=eps1m(p)/eps1m(p+1);
      mu1mratio(p)=mu1m(p)/mu1m(p+1);
      nu2mratio(p)=nu2m(p+1)/nu2m(p);
      eps1m(p)=1/eps1m(p);
      mu1m(p)=1/mu1m(p);
    }
    eps1m(M+1)=1/eps1m(M+1);
    mu1m(M+1)=1/mu1m(M+1);
    
        
    // parameters
    param.M=M+1;
    param.N=N+1;
    param.b=b;
    param.bm=bm;
    param.beta=beta;
    param.eps1=eps1;
    param.eps1m=eps1m;
    param.eps1ratio=eps1ratio;
    param.eps1mratio=eps1mratio;
    param.mu1=mu1;
    param.mu1m=mu1m;
    param.mu1ratio=mu1ratio;
    param.mu1mratio=mu1mratio;
    param.nu2=nu2;
    param.nu2m=nu2m;
    param.nu2ratio=nu2ratio;
    param.nu2mratio=nu2mratio;
    param.kovergamma=kovergamma;
    param.m=m;
    param.n=n;

    // integration with adaptative integration on infinite interval (QAGIU GSL algorithm)
    F.params=&param;

    if (flagreal) {
      // real part
      F.function=&integrand_real;
    } else {
      // imaginary part
      F.function=&integrand_imag;
    }

    status=gsl_integration_qagiu(&F, 0, 0., tolint, 15, w, &x, &err);
    //status=gsl_integration_qagiu(&F, 0, tolintabs, tolint, 15, w, &x, &err);

    // deal with GSL errors
    if (status) {
      //printf("alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));

      if ( ( (status==GSL_EMAXITER)||(status==GSL_EDIVERGE) )||(status==GSL_EROUND) ) {
	if (flagreal) F.function=&integrand_real_modif;
	else F.function=&integrand_imag_modif;
        status=gsl_integration_qag(&F, 0., 1., 0., tolint, limit, 1, w, &x, &err);
        //status=gsl_integration_qag(&F, 0., 1., tolintabs, tolint, limit, 1, w, &x, &err);
        //printf("GSL_EDIVERGE: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
	if (status) printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
      /*else if (status==GSL_EROUND) {
	while ((std::abs(err/x)>2*tolint)&&((status==GSL_EROUND)||(status==GSL_EDIVERGE))) {
          tolint*=2;tolintabs*=2;
	  status=gsl_integration_qagiu(&F, 0, tolintabs, tolint, limit, w, &x, &err);
	}
	printf("GSL_EROUND: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }*/
      else {
        printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
    } 
    
    /*if (std::abs(x)<tolintabs) {
      status=gsl_integration_qagiu(&F, 0, std::abs(x), tolint, limit, w, &x, &err);
      if (status) {
        printf("Warning: alpha%d%d, omega= %s, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",m,n,omega.toDec().c_str(),flagreal,x,std::abs(err/x));
      }
    }*/

    return x;
     
  }



/**************************************************
*** read_input: read an input file line
***
***************************************************/

  void read_input(std::string line, std::string description, unsigned int& param0, double& param1, 
  	std::string& param2, amp::ampf<Precision>& param3, int type){
	
    /* function reading the number or string at the end of an input line, identifying the description string.
    type=0 if output is an unsigned int, 1 if output is a double, 2 if output is a char array, and 3 if output is an
    amp high precision number.
    Output is in "param[type_number]". */
    
    size_t found;

    if (line.find(description) != std::string::npos){
      found=line.find("\t");

      if (found != std::string::npos){
         
        switch(type){
	  
	  case 0:
    
	    param0=atoi(line.substr(found+1,std::string::npos).c_str());
	    break;
	    
	  case 1:
	  
	    param1=strtod(line.substr(found+1,std::string::npos).c_str(),NULL);
	    break;
	    
	  case 2:
	  
	    param2=line.substr(found+1,std::string::npos);
	    break;
	    
	  case 3:
	  
	    param3=line.substr(found+1,std::string::npos).c_str();
	    break;
	
	}

      }

    }
    
    return;

  }


/**************************************************
*** read_input_layer: read an input file line
*** for a specific "layer property" line
***************************************************/

  void read_input_layer(std::string line, std::string description, int p, amp::ampf<Precision>& param){
	
    /* function reading the number or string at the end of an input line for a line containing a property of layer p.
    Output is an amp high precision number, in param. */
    
    unsigned int dummy0;
    double dummy1;
    std::string dummy2;
    
    std::ostringstream des;  
    des << "Layer " << p << " " << description;
    //cout << des.str() << "\n";
    //cout << line << "\n";
    read_input(line,des.str(),dummy0,dummy1,dummy2,param,3);
    
    return;
    
  }


/**************************************************
*** impedance: computes the impedances (std precision)
***
***************************************************/

  void impedance(complex<double>& Zxdip, complex<double>& Zydip,
  	complex<double>& Zyquad,
	complex<double>& Zlong, complex<double>& Zycst,
  	unsigned int N, unsigned int M,
	ap::template_1d_array< amp::ampf<Precision> >& rho, 
  	ap::template_1d_array< amp::ampf<Precision> >& tau, 
  	ap::template_1d_array< amp::ampf<Precision> >& epsb, 
  	ap::template_1d_array< amp::ampf<Precision> >& chi, 
  	ap::template_1d_array< amp::ampf<Precision> >& fmu, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	ap::template_1d_array< amp::ampf<Precision> >& rhom, 
  	ap::template_1d_array< amp::ampf<Precision> >& taum, 
  	ap::template_1d_array< amp::ampf<Precision> >& epsbm, 
  	ap::template_1d_array< amp::ampf<Precision> >& chim, 
  	ap::template_1d_array< amp::ampf<Precision> >& fmum, 
	ap::template_1d_array< amp::ampf<Precision> >& bm,
	amp::ampf<Precision> beta, amp::ampf<Precision> gamma, 
	int flag_topbotsym, double L, double freq) {
	
    /* function that computes the impedances at a given frequency freq, from rho, tau, epsb, chi, fmu, b 
    of each of the N upper layers and rhom, taum, epsbm, chim, fmum, bm of each of the M lower layers, 
    and from the relativistic velocity factor beta. */

    amp::ampf<Precision> omega, k, kovergamma, mu0, eps0, Z0;
    ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,mu1,mu1m; /* eps1 and mu1 for upper (without m at
 					the end) and lower layers (with m at the end) */
    amp::campf<Precision> jimagMP; // imaginary constant in multiprecision
    //gsl_integration_workspace *w;
    double x,y,errx,erry;
    double tolintabs=1.e-15,tolint=1.e-6; // absolute and relative error permitted for gsl adaptative integration
    complex<double> cst,alpha00,alpha01,alpha02,alpha11,jimag; // some constant, alphamn constants and imaginary constant
    complex<double> tolintabs00, tolintabs01, tolintabs02, tolintabs11;
    size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
    unsigned long lprov;
    
    // try to find freq in the table freqmem
    //cout.precision(18);
    //cout << "impedance: Frequency " << freq << "\n";
    //cout.flush();
    if (impmem==0) lprov=0;
    else lprov=locate(freqmem,freq,impmem-1);
 
    if ( (impmem!=0)&&(freq==freqmem[lprov]) ) {
      Zxdip=Zxdipmem[lprov];Zydip=Zydipmem[lprov];
      Zyquad=Zyquadmem[lprov];Zlong=Zlongmem[lprov];
      Zycst=Zycstmem[lprov];
      //cout << "Frequency found in impedance memory " << lprov << "\n";
    }
    else if ( ( (impmem!=0)&&(lprov>0)) && (freq==freqmem[lprov-1]) ) {
      Zxdip=Zxdipmem[lprov-1];Zydip=Zydipmem[lprov-1];
      Zyquad=Zyquadmem[lprov-1];Zlong=Zlongmem[lprov-1];
      Zycst=Zycstmem[lprov-1];
      //cout << "Frequency found in impedance memory " << lprov-1 << "\n";
    }
    
    
    else {

      //cout << "Frequency not found in impedance memory\n";
      //cout.flush();

      if (impmem==0) {
	tolintabs00=complex<double>(0.,0.);
	tolintabs01=complex<double>(0.,0.);
	tolintabs02=complex<double>(0.,0.);
	tolintabs11=complex<double>(0.,0.);
      }
      else if ( (freqmem[lprov]-freq)<(freq-freqmem[lprov-1]) ) {
        tolintabs00=tolintabs00mem[lprov];
        tolintabs01=tolintabs01mem[lprov];
        tolintabs02=tolintabs02mem[lprov];
        tolintabs11=tolintabs11mem[lprov];
      }
      else {
        tolintabs00=tolintabs00mem[lprov-1];
        tolintabs01=tolintabs01mem[lprov-1];
        tolintabs02=tolintabs02mem[lprov-1];
        tolintabs11=tolintabs11mem[lprov-1];        
      }
      
      //cout << "lprov=" << lprov << "\n";
      //cout.flush();

      // constants
      mu0=4e-7*amp::pi<Precision>();
      eps0=1/(mu0*amp::sqr(C));
      Z0=mu0*C;
      jimagMP.x=0;jimagMP.y=1;
      jimag=complex<double>(0.,1.);

      // workspace allocation for gsl adaptative integration
      gsl_integration_workspace *w=gsl_integration_workspace_alloc(limit);
      //for (unsigned int i=0;i<500;i++) cout << "tata ";cout << "\n";cout.flush();

      // allocation for eps1, mu1, eps1m and mu1m
      eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);
      eps1m.setbounds(1,M+1);mu1m.setbounds(1,M+1);

      //for (unsigned int i=0;i<500;i++) cout << "tata\n";
      mem=0; // initialize memory (in kx) at each frequency
      omega=amp::twopi<Precision>()*freq;
      k=omega/(beta*C);
      kovergamma=k/gamma;

      // first layer (inside the chamber) is always vacuum
      eps1(1)=1;mu1(1)=1;
      eps1m(1)=1;mu1m(1)=1;
      
      // computes the layer properties for the angular freq. omega
      for (unsigned int p=2;p<=N+1; p++) {
        if (rho(p).isFiniteNumber()) {
	  eps1(p)=epsb(p)+1/(jimagMP*eps0*rho(p)*omega*(1+jimagMP*omega*tau(p)));
	} else {
	  eps1(p)=epsb(p);
	}
	mu1(p)=1+chi(p)/(1+jimagMP*omega/(amp::twopi<Precision>()*fmu(p)));
	//cout << p << "\n";
	//printf("%s %s %s\n%s %s\n",b(p).toDec().c_str(),eps1(p).x.toDec().c_str(),eps1(p).y.toDec().c_str(),
	//		mu1(p).x.toDec().c_str(),mu1(p).y.toDec().c_str());
      }
      for (unsigned int p=2;p<=M+1; p++) {
	if (rhom(p).isFiniteNumber()) {
	  eps1m(p)=epsbm(p)+1/(jimagMP*eps0*rhom(p)*omega*(1+jimagMP*omega*taum(p)));
	} else {
	  eps1m(p)=epsbm(p);
	}
	mu1m(p)=1+chim(p)/(1+jimagMP*omega/(amp::twopi<Precision>()*fmum(p)));
	//printf("%s %s %s\n%s %s\n",bm(p).toDec().c_str(),eps1m(p).x.toDec().c_str(),eps1m(p).y.toDec().c_str(),
	//		mu1m(p).x.toDec().c_str(),mu1m(p).y.toDec().c_str());
      }

      //printf("%s %s\n",b(1).toDec().c_str(),bm(1).toDec().c_str());
      //printf("%s %s %s\n",omega.toDec().c_str(),k.toDec().c_str(),kovergamma.toDec().c_str());
      //for (unsigned int i=0;i<500;i++) cout << "tata ";cout << "\n";cout.flush();

      // computes alpha00
      x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.real(),limit,w);
      y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.imag(),limit,w);
      alpha00=complex<double>(x,y);
      //printf("alpha00: %13.8e %13.8e\n",alpha00.real(),alpha00.imag());

      // computes alpha01
      if (flag_topbotsym==0) {
	// computes alpha01
	x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.real(),limit,w);
	y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.imag(),limit,w);
	alpha01=complex<double>(x,y);
	}
      else alpha01=complex<double>(0.,0.);
      //printf("alpha01: %13.8e %13.8e\n",alpha01.real(),alpha01.imag());

      // computes alpha02
      x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,tolintabs02.real(),limit,w);
      y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,tolintabs02.imag(),limit,w);
      alpha02=complex<double>(x,y);
      //printf("alpha02: %13.8e %13.8e\n",alpha02.real(),alpha02.imag());

      // computes alpha11
      x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.real(),limit,w);
      y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.imag(),limit,w);
      alpha11=complex<double>(x,y);
      //printf("alpha11: %13.8e %13.8e\n",alpha11.real(),alpha11.imag());


      // computes the impedances
      cst=jimag*L*double(amp::ampf<Precision>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<Precision>())).toDouble());
      Zlong=cst*alpha00;
      Zycst=cst*alpha01/double(amp::ampf<Precision>(gamma).toDouble());
      cst=cst*double(amp::ampf<Precision>(k/(2*amp::sqr(gamma))).toDouble());
      Zxdip=cst*(alpha02-alpha00);
      Zydip=2.*cst*alpha11;
      //Zxquad=-Zxdip; //unused
      Zyquad=cst*(alpha02+alpha00);

      //cout.precision(18);
      //cout << freq << " " << Zlong << " " << Zxdip << " " << Zydip << " " << Zyquad << " " << Zycst << '\n';
      //cout << freq << " " << Zxdip.real() << " " << Zxdip.imag() << '\n';
      //cout << cst << " " << alpha02.real() << " " << alpha02.imag() << " " << alpha00.real() << " " << alpha00.imag() << '\n';
      //cout.flush();

      gsl_integration_workspace_free(w);
 
      // add to the memory
      for(unsigned int k=impmem; k>=lprov+1; k--) {
	freqmem[k]=freqmem[k-1];
	Zxdipmem[k]=Zxdipmem[k-1];Zydipmem[k]=Zydipmem[k-1];
	Zyquadmem[k]=Zyquadmem[k-1];Zycstmem[k]=Zycstmem[k-1];
	Zlongmem[k]=Zlongmem[k-1];
	tolintabs00mem[k]=tolintabs00mem[k-1];tolintabs01mem[k]=tolintabs01mem[k-1];
	tolintabs02mem[k]=tolintabs02mem[k-1];tolintabs11mem[k]=tolintabs11mem[k-1];
      }
      freqmem[lprov]=freq;
      Zxdipmem[lprov]=Zxdip;Zydipmem[lprov]=Zydip;
      Zyquadmem[lprov]=Zyquad;Zycstmem[lprov]=Zycst;
      Zlongmem[lprov]=Zlong;
      tolintabs00mem[lprov]=tolintabs*complex<double>(abs(alpha00.real()),abs(alpha00.imag()));  
      tolintabs01mem[lprov]=tolintabs*complex<double>(abs(alpha01.real()),abs(alpha01.imag()));  
      tolintabs02mem[lprov]=tolintabs*complex<double>(abs(alpha02.real()),abs(alpha02.imag()));  
      tolintabs11mem[lprov]=tolintabs*complex<double>(abs(alpha11.real()),abs(alpha11.imag()));
      impmem++;
       
    }

    return;

  }


/**************************************************
*** pchip: interpolated value at z, given the function 
*** values fi at the points zi and the slopes of the 
*** interpolated polynomial di at zi. Tables have size nz.
***************************************************/

  complex<double> pchip(double z, double *zi, 
  	complex<double> *fi, complex<double> *di, unsigned long nz) {
	
    double delta,t,t2,t3,tm,tm2,tm3,h1,h2,h3,h4;
    complex<double> value;
    unsigned long l;
    
    l=locate(zi, z, nz-1);
    if (z==zi[0]) l=1; /* change value from 0 to 1 in this case (because of the slightly modifed 
    		version of locate we use here) */
		
    delta=zi[l]-zi[l-1];
    //cout << l << " " << delta << "\n";
    t=(zi[l]-z)/delta;
    tm=(z-zi[l-1])/delta;
    t2=t*t;t3=t*t2;
    tm2=tm*tm;tm3=tm*tm2;
    
    // cubic Hermite basis functions
    h1=3.*t2-2.*t3;
    h2=3.*tm2-2.*tm3;
    h3=-delta*(t3-t2);
    h4=delta*(tm3-tm2);
    
    // interpolated value
    value=fi[l-1]*h1+fi[l]*h2+di[l-1]*h3+di[l]*h4;
	
    return value;

  }


/**************************************************
*** interp: interpolated value at z, given the function 
*** values fi at the points zi and the slopes of the 
*** interpolated polynomial di at zi (in case of pchip interpolation).
*** Tables have size nz.
*** interp_type is the type of interpolation: 
*** 0 for pchip, 1 for linear, 2 for exponential, 3 for power law
***************************************************/

  complex<double> interp(double z, double *zi, complex<double> *fi, 
  	complex<double> *di, unsigned long nz, unsigned int interp_type) {
	
    complex<double> value;
    unsigned long l;
    
    switch (interp_type) {
      case 0:
        value=pchip(z, zi, fi, di, nz);
	return value;
    
      case 1:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; /* change value from 0 to 1 in this case (because of the slightly modifed 
    		version of locate we use here) */
        // linearly interpolated value
        value=fi[l-1]+(fi[l]-fi[l-1])*(z-zi[l-1])/(zi[l]-zi[l-1]);
	return value;
	
      /*case 2:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; // change value from 0 to 1 in this case (because of the slightly modifed 
    			   // version of locate we use here)
        // log of the function is linearly interpolated
        value=exp(log(fi[l-1])+(log(fi[l])-log(fi[l-1]))*(z-zi[l-1])/(zi[l]-zi[l-1]));
	return value;

      case 3:
        l=locate(zi, z, nz-1);
        if (z==zi[0]) l=1; // change value from 0 to 1 in this case (because of the slightly modifed 
    			   //version of locate we use here)
        // log of the function w.r.t to log of the variable is linearly interpolated
        value=exp(log(fi[l-1])+(log(fi[l])-log(fi[l-1]))*(log(z)-log(zi[l-1]))/(log(zi[l])-log(zi[l-1])));
	return value;
      */
    }
    
  }


/**************************************************
*** pchipslope: slopes of the pchip interpolating
*** polynomial (Source: Monotone Piecewise Cubic Interpolation, 
*** F. N. Fritsch & R. E. Carlson, SIAM Journal on Numerical Analysis,
*** Vol. 17, No. 2, 1980, pp. 238-246 )
***************************************************/

  void pchipslope(complex<double>* d, double* x, complex<double> * y, 
	  unsigned long length){
    /* Derivative values for monotonic Piecewise Cubic Hermite Interpolation (pchip)
      d = pchipslope(x,y,length) computes the first derivatives, d[k] = P'(x[k]), of the interpolating polynomial*/

    // computes first the "finite differences" derivatives of y with respect to x
    complex<double>* del=new complex<double>[length-1];
    for (unsigned long i=0; i<length-1; i++) {
      del[i] = (y[i+1] - y[i])/(x[i+1]-x[i]);
    }

    //  Special case n=2, use linear interpolation.
    if (length==2) {
      //complex<double>* d = new complex<double>[2];
      d[0]=del[0];
      d[1]=del[0];
      return;
    }

    //  Slopes at interior points.
    //  d(k) = weighted average of del(k-1) and del(k) when they have the same sign.
    //  d(k) = 0 when del(k-1) and del(k) have opposites signs or either is zero.

    //complex<double>* d = new complex<double>[length];
    for (unsigned long i=0; i<length; i++) {
      d[i]=0.;
    }

    std::vector<unsigned long> k;
    for (unsigned long i=0; i<length-2; i++) {
      if (abs(del[i])!=0. || abs(del[i+1])!=0.) {
	k.push_back(i);
      }
    }

    double h[length-1];
    for (unsigned long i=0; i<length-1; i++){
      h[i]=x[i+1]-x[i];
    }

    double hs[k.size()];
    double w1[k.size()];
    double w2[k.size()];
    for (unsigned long i=0; i<k.size(); i++){
      hs[i] = h[k[i]] + h[k[i]+1];
      w1[i] = ( h[ k[i] ] + hs[i]) / (3.*hs[i]);
      w2[i] = ( h[ k[i]+1 ] + hs[i]) / (3.*hs[i]);
    }

    double dmax[k.size()];
    double dmin[k.size()];

    for (unsigned long i=0; i<k.size(); i++){
      dmax[i] = max ( abs(del[k[i]]), abs(del[k[i]+1] )); 
      dmin[i] = min ( abs(del[k[i]]), abs(del[k[i]+1] )); 
      d[k[i]+1] = dmin[i] / std::conj( w1[i] * del[k[i]] / dmax[i] + w2[i]*del[k[i]+1]/ dmax[i]);
      /*cout << "i= " << i << " ; dmax= " << dmax[i] << " ; dmin= " << dmin[i] << '\n' ;
      cout << "k[i]= " << k[i] << " ; w1= " << w1[i] << " ; w2= " << w2[i] << '\n' ;
      cout << "del[k[i]]= " << del[k[i]] << " ; del[k[i]+1]= " << del[k[i]+1] << " ; d[k[i]+1]= " << d[k[i]+1] << '\n' ;
      cout << '\n' ;*/
    }

    //  Slopes at end points.
    //  Set d(1) and d(n) via non-centered, shape-preserving three-point formulae.

    d[0] = (((double)2.*h[0]+h[1])*del[0] - h[0]*del[1])/(h[0]+h[1]);

    //if isreal(d) && (sign(d(1)) ~= sign(del(1)))
    //d(1) = 0;
    if (del[0]/abs(del[0]) != del[1]/abs(del[1]) && (abs(d[0]) > (double)3.*abs(del[0]) )) {
      d[0] = (double)3.*del[0];
    }

    d[length-1] = (((double)2.*h[length-2]+h[length-3])*del[length-2] - h[length-2]*del[length-3])/(h[length-2]+h[length-3]);

    if (del[length-2]/abs(del[length-2]) != del[length-3]/abs(del[length-3]) && abs(d[length-1]) > abs((double)3.*del[length-2])) {
      d[length-1] = (double)3.*del[length-2];
    }
    
    delete[] del;
    return;
  }

/***********************************************************
*** integrand_diff: function computing
*** the difference in norm between
*** the impedances at freq=exp(u) and their interpolation
***********************************************************/

  double integrand_diff(double u, void *p){

    /* function encapsulating the computation of the difference in norm between the 
    impedances and their pchip interpolation */

    struct params_diff *param=(struct params_diff *)p;
    unsigned int N,M,interp_type; // number of upper and lower layers, and interpolation type
    unsigned long nf; // number of frequencies in interpolation 
    // multilayer parameters
    ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,chi,fmu,b,rhom,taum,epsbm,chim,fmum,bm;
    amp::ampf<Precision> beta,gamma; // relativistic velocity factor
    double *freqi; // frequencies of the interpolation
    complex<double> *Zxdipfi,*Zxdipdi; // values and derivatives of the interpolation at the frequencies freqi (Zxdip)
    complex<double> *Zydipfi,*Zydipdi; // values and derivatives of the interpolation at the frequencies freqi (Zydip)
    complex<double> *Zyquadfi,*Zyquaddi; // values and derivatives of the interpolation at the frequencies freqi (Zyquad)
    complex<double> *Zlongfi,*Zlongdi; // values and derivatives of the interpolation at the frequencies freqi (Zlong)
    complex<double> *Zycstfi,*Zycstdi; // values and derivatives of the interpolation at the frequencies freqi (Zycst)
    complex<double> Zxdip,Zydip,Zyquad,Zlong,Zycst; // exact impedances
    complex<double> pZxdip,pZydip,pZyquad,pZlong,pZycst; // interpolated impedances
    double L,result,freq;
    int flag_topbotsym;
    
    freq=exp(u);

    M=(param->M);N=(param->N);
    rho=(param->rho);rhom=(param->rhom);
    tau=(param->tau);taum=(param->taum);
    epsb=(param->epsb);epsbm=(param->epsbm);
    chi=(param->chi);chim=(param->chim);
    fmu=(param->fmu);fmum=(param->fmum);
    b=(param->b);bm=(param->bm);
    beta=(param->beta);
    gamma=(param->gamma);
    flag_topbotsym=(param->flag_topbotsym);
    L=(param->L);

    freqi=(param->freqi);
    interp_type=(param->interp_type);
    nf=(param->nf);
    Zxdipfi=(param->Zxdipfi);Zxdipdi=(param->Zxdipdi);
    Zydipfi=(param->Zydipfi);Zydipdi=(param->Zydipdi);
    Zyquadfi=(param->Zyquadfi);Zyquaddi=(param->Zyquaddi);
    Zycstfi=(param->Zycstfi);Zycstdi=(param->Zycstdi);
    Zlongfi=(param->Zlongfi);Zlongdi=(param->Zlongdi);
    

    // computes the impedances at freq
    impedance(Zxdip, Zydip, Zyquad, Zlong, Zycst, N, M, rho, tau, epsb, chi, fmu, b, 
	rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, freq);

    pZxdip=interp(freq,freqi,Zxdipfi,Zxdipdi,nf,interp_type);
    pZydip=interp(freq,freqi,Zydipfi,Zydipdi,nf,interp_type);
    pZyquad=interp(freq,freqi,Zyquadfi,Zyquaddi,nf,interp_type);
    pZlong=interp(freq,freqi,Zlongfi,Zlongdi,nf,interp_type);
    pZycst=interp(freq,freqi,Zycstfi,Zycstdi,nf,interp_type);

    /*result=max(abs(Zxdip-pZxdip),abs(Zydip-pZydip));
    result=max(result,max(abs(Zxquad-pZxquad),abs(Zyquad-pZyquad)));
    result=max(result,max(factlong*abs(Zlong-pZlong),abs(Zycst-pZycst)));*/
    result=abs(Zxdip-pZxdip)*abs(Zxdip-pZxdip)+abs(Zydip-pZydip)*abs(Zydip-pZydip);
    result+=abs(Zyquad-pZyquad)*abs(Zyquad-pZyquad)+abs(Zycst-pZycst)*abs(Zycst-pZycst);
    result+=factlong*factlong*abs(Zlong-pZlong)*abs(Zlong-pZlong);
    /*result=abs(Zxdip.real()-pZxdip.real())+abs(Zxdip.imag()-pZxdip.imag());
    result=max(result,abs(Zydip.real()-pZydip.real())+abs(Zydip.imag()-pZydip.imag()));
    result=max(result,abs(Zyquad.real()-pZyquad.real())+abs(Zyquad.imag()-pZyquad.imag()));
    result=max(result,abs(Zycst.real()-pZycst.real())+abs(Zycst.imag()-pZycst.imag()));
    result=max(result,factlong*(abs(Zlong.real()-pZlong.real())+abs(Zlong.imag()-pZlong.imag())));*/

    //cout.precision(18);
    //cout << u << "\t" << std::sqrt(result)*freq << '\n';
    //cout << freq << "\t" << Zlong.real() << " " << Zlong.imag() << " " << Zxdip.real() << " "<< Zxdip.imag() << " " << Zydip.real() << " "<< Zydip.imag() << " " << Zyquad.real()<< " " << Zyquad.imag() << " " << Zycst.real() << " " << Zycst.imag() << " " << '\n';
    //cout << freq << "\t" << pZlong.real() << " " << pZlong.imag() << " " << pZxdip.real() << " "<< pZxdip.imag() << " " << pZydip.real() << " "<< pZydip.imag() << " " << pZyquad.real()<< " " << pZyquad.imag() << " " << pZycst.real() << " " << pZycst.imag() << " " << '\n';
    //cout.flush();
    
    return std::sqrt(result)*freq;
    
  }

/***********************************************************
*** integrand_diff_freq: function computing
*** the difference in norm between
*** the impedances at freq and their interpolation.
*** variable is not u=log(freq), it is here freq directly
***********************************************************/

  double integrand_diff_freq(double freq, void *p){

    /* function encapsulating the computation of the difference in norm between the 
    impedances and their pchip interpolation */
    double u;
    
    u=log(freq);
    return integrand_diff(u,p)/freq;
    
  }

/***********************************************************
*** integrand_diff2: function computing
*** the difference in norm between
*** the impedances at freq=exp(u) and their interpolation, minus
*** a certain vector x
***********************************************************/

  double integrand_diff2(double u, void *p){

    /* function encapsulating the computation of the difference in norm between the 
    impedances and their pchip interpolation */

    struct params_diff *param=(struct params_diff *)p;
    unsigned int N,M,interp_type; // number of upper and lower layers, and interpolation type
    unsigned long nf; // number of frequencies in interpolation 
    // multilayer parameters
    ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,chi,fmu,b,rhom,taum,epsbm,chim,fmum,bm;
    amp::ampf<Precision> beta,gamma; // relativistic velocity factor
    double *freqi; // frequencies of the interpolation
    complex<double> *Zxdipfi,*Zxdipdi; // values and derivatives of the interpolation at the frequencies freqi (Zxdip)
    complex<double> *Zydipfi,*Zydipdi; // values and derivatives of the interpolation at the frequencies freqi (Zydip)
    complex<double> *Zyquadfi,*Zyquaddi; // values and derivatives of the interpolation at the frequencies freqi (Zyquad)
    complex<double> *Zlongfi,*Zlongdi; // values and derivatives of the interpolation at the frequencies freqi (Zlong)
    complex<double> *Zycstfi,*Zycstdi; // values and derivatives of the interpolation at the frequencies freqi (Zycst)
    complex<double> Zxdip,Zydip,Zyquad,Zlong,Zycst; // exact impedances
    complex<double> pZxdip,pZydip,pZyquad,pZlong,pZycst; // interpolated impedances
    complex<double> *x; // vector to take away from the difference
    double L,result,freq;
    int flag_topbotsym;

    x=new complex<double>[5];
    
    freq=exp(u);

    M=(param->M);N=(param->N);
    rho=(param->rho);rhom=(param->rhom);
    tau=(param->tau);taum=(param->taum);
    epsb=(param->epsb);epsbm=(param->epsbm);
    chi=(param->chi);chim=(param->chim);
    fmu=(param->fmu);fmum=(param->fmum);
    b=(param->b);bm=(param->bm);
    beta=(param->beta);
    gamma=(param->gamma);
    flag_topbotsym=(param->flag_topbotsym);
    L=(param->L);

    freqi=(param->freqi);
    interp_type=(param->interp_type);
    nf=(param->nf);
    Zxdipfi=(param->Zxdipfi);Zxdipdi=(param->Zxdipdi);
    Zydipfi=(param->Zydipfi);Zydipdi=(param->Zydipdi);
    Zyquadfi=(param->Zyquadfi);Zyquaddi=(param->Zyquaddi);
    Zycstfi=(param->Zycstfi);Zycstdi=(param->Zycstdi);
    Zlongfi=(param->Zlongfi);Zlongdi=(param->Zlongdi);
    
    x=(param->x);
    

    // computes the impedances at freq
    impedance(Zxdip, Zydip, Zyquad, Zlong, Zycst, N, M, rho, tau, epsb, chi, fmu, b, 
	rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, freq);

    pZxdip=interp(freq,freqi,Zxdipfi,Zxdipdi,nf,interp_type);
    pZydip=interp(freq,freqi,Zydipfi,Zydipdi,nf,interp_type);
    pZyquad=interp(freq,freqi,Zyquadfi,Zyquaddi,nf,interp_type);
    pZlong=interp(freq,freqi,Zlongfi,Zlongdi,nf,interp_type);
    pZycst=interp(freq,freqi,Zycstfi,Zycstdi,nf,interp_type);

    /*result=max(abs(Zxdip-pZxdip),abs(Zydip-pZydip));
    result=max(result,max(abs(Zxquad-pZxquad),abs(Zyquad-pZyquad)));
    result=max(result,max(factlong*abs(Zlong-pZlong),abs(Zycst-pZycst)));*/
    result=abs(Zxdip-pZxdip-x[0])*abs(Zxdip-pZxdip-x[0])+abs(Zydip-pZydip-x[1])*abs(Zydip-pZydip-x[1]);
    result+=abs(Zyquad-pZyquad-x[2])*abs(Zyquad-pZyquad-x[2])+abs(Zycst-pZycst-x[3])*abs(Zycst-pZycst-x[3]);
    result+=factlong*abs(Zlong-pZlong-x[4])*abs(Zlong-pZlong-x[4]);

    //cout.precision(18);
    //cout << u << "\t" << std::sqrt(result)*freq << '\n';
    //cout << freq << "\t" << Zlong.real() << " " << Zlong.imag() << " " << Zxdip.real() << " "<< Zxdip.imag() << " " << Zydip.real() << " "<< Zydip.imag() << " " << Zyquad.real()<< " " << Zyquad.imag() << " " << Zycst.real() << " " << Zycst.imag() << " " << '\n';
    //cout << freq << "\t" << pZlong.real() << " " << pZlong.imag() << " " << pZxdip.real() << " "<< pZxdip.imag() << " " << pZydip.real() << " "<< pZydip.imag() << " " << pZyquad.real()<< " " << pZyquad.imag() << " " << pZycst.real() << " " << pZycst.imag() << " " << '\n';
    //cout.flush();
    
    return std::sqrt(result)*freq;
    
  }

/***********************************************************
*** mean_diff: function computing
*** the average values of the differences between
*** the impedances in the memory between lprov and lprov2,
*** and their pchip interpolation.
***********************************************************/

  complex<double>* mean_diff(unsigned long lprov, unsigned long lprov2, void *p){
  
    struct params_diff *param=(struct params_diff *)p;
    unsigned long nf; // number of frequencies in interpolation 
    double *freqi; // frequencies of the interpolation
    complex<double> *x;
    complex<double> *Zxdipfi,*Zxdipdi; // values and derivatives of the interpolation at the frequencies freqi (Zxdip)
    complex<double> *Zydipfi,*Zydipdi; // values and derivatives of the interpolation at the frequencies freqi (Zydip)
    complex<double> *Zyquadfi,*Zyquaddi; // values and derivatives of the interpolation at the frequencies freqi (Zyquad)
    complex<double> *Zlongfi,*Zlongdi; // values and derivatives of the interpolation at the frequencies freqi (Zlong)
    complex<double> *Zycstfi,*Zycstdi; // values and derivatives of the interpolation at the frequencies freqi (Zycst)

    x=new complex<double>[5];
    
    freqi=(param->freqi);
    nf=(param->nf);
    Zxdipfi=(param->Zxdipfi);Zxdipdi=(param->Zxdipdi);
    Zydipfi=(param->Zydipfi);Zydipdi=(param->Zydipdi);
    Zyquadfi=(param->Zyquadfi);Zyquaddi=(param->Zyquaddi);
    Zycstfi=(param->Zycstfi);Zycstdi=(param->Zycstdi);
    Zlongfi=(param->Zlongfi);Zlongdi=(param->Zlongdi);
    
    for (unsigned long k=0; k<=4; k++) x[k]=complex<double>(0.,0.);
    
    for (unsigned long k=lprov; k<=lprov2; k++) {
      x[0]+=Zxdipmem[k]-pchip(freqmem[k],freqi,Zxdipfi,Zxdipdi,nf);
      x[1]+=Zydipmem[k]-pchip(freqmem[k],freqi,Zydipfi,Zydipdi,nf);
      x[2]+=Zyquadmem[k]-pchip(freqmem[k],freqi,Zyquadfi,Zyquaddi,nf);
      x[3]+=Zlongmem[k]-pchip(freqmem[k],freqi,Zlongfi,Zlongdi,nf);
      x[4]+=Zycstmem[k]-pchip(freqmem[k],freqi,Zycstfi,Zycstdi,nf);
    }
    for (unsigned long k=0; k<=4; k++) x[k]/=(double)(lprov2-lprov+1);
    
    return x;
    
  }

/***********************************************************
*** Function to compute Phi(x)
***********************************************************/

  complex<long double> Phi(complex<long double> x,double eps){
    
    complex<long double> Phi(0.L,0.L),power(1.L,0.L);
    int condition;
    int i=0;
    long double fact=1.L,fact2,abspower,il,expoabs;
    
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = (il+6.L)/((il+3.L)*(il+4.L));
            Phi+= power * fact * fact2;
	    power *= x;
	    abspower = abs(power);
	    fact = fact/(il+1.L);
	    //cout << abs(std::real(Phi)) << " " << abs(std::imag(Phi)) <<
	    //		" " << min(abs(std::real(Phi)),abs(std::imag(Phi))) << '\n';
      	    condition =  (2.L* abspower * expoabs * fact / (il+5.L) > 
	    	(long double)eps*min(abs(std::real(Phi)),abs(std::imag(Phi))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Phi= (std::pow(x,3)*expo - 6.L*x*(expo+1.L) + 12.L*(expo-1.L))/(std::pow(x,4));
    }
    return Phi;
  }
  

/***********************************************************
*** Function to compute Psi(x)
***********************************************************/

  complex<long double> Psi(complex<long double> x,double eps){
    
    complex<long double> Psi(0.L,0.L),power(1.L,0.L);
    int condition;
    long double fact=1.L,fact2,abspower,il,expoabs;
    int i=0;
     
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = 1.L/((il+3.L)*(il+4.L));
	    Psi+= - power * fact * fact2;
	    power *= x;
	    fact = fact/(il+1.L);
	    abspower = abs(power);
      	    condition =  (abspower * expoabs * fact / ((il+4.L)*(il+5.L)) >
	    	(long double)eps*min(abs(std::real(Psi)),abs(std::imag(Psi))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Psi = (-std::pow(x,2)*expo + 2.L*x*(2.L*expo+1.L) - 6.L*(expo-1.L))/(std::pow(x,4));
    }
    return Psi;
  }
  
/***********************************************************
*** Function to compute Lambda(x)
***********************************************************/

  complex<long double> Lambda(complex<long double> x,double eps){
    
    complex<long double> Lambda(0.L,0.L),power(1.L,0.L);
    int condition;
    long double fact=1.L,fact2,abspower,il,expoabs;
    int i=0;
     
    if (abs(x)<1.L) {
      expoabs=exp(abs(x));
      do {  
            // fact is 1/factorial(i) and power is x^i
	    il=(long double)i;
	    fact2 = 1.L/(il+2.L);
	    Lambda+= power * fact * fact2;
	    power *= x;
	    abspower = abs(power);
      	    fact = fact/(il+1.L);
	    condition =  (abspower * expoabs * fact/(il+3.L) >
	    	(long double)eps*min(abs(std::real(Lambda)),abs(std::imag(Lambda))));
      	    i++;} while ( condition ) ;
    }
    else {
      complex<long double> expo = exp(x);
      Lambda = (x*expo - expo+1.L)/(x*x);
    }
    return Lambda;
  }

/***********************************************************
*** Function to compute Fourier integral on semi-infinite domain
***********************************************************/

complex<double> fourier_integral_inf(complex<double>* fi,complex<double>* d,complex<double>
	df, long double t, long double* omegai, long double* delta, unsigned long length, double eps,
	unsigned int* interp_type, int flaginf) {

// Computes the Fourier integral at t from omegai(1) to omegai(end) of f(omega)*exp(j*omega*t)
// where f is a function. Add a correcting term to take into account the
// rest of the integral between omegai(end) and infinity (or -infinity and
// -omegai(1) if omegai(end)<0)

// In input (omegai,fi) are the abscissae where f is interpolated to perform the 
// calculation, and the corresponding values of f. df is the derivative of f at
// the freqi which is maximum in absolute value.
// We use a cubic interpolation between the points omegai (piecewise monotonic - pchip). 
// The slopes of the interpolation
// are in d (of size length). delta (of size length-1) are the differences
// between sucessive freqi. eps is the relative precision in the computation of the auxiliary
// functions Phi and Psi.
// THe type of interpolation on each interval is in interp_type (0 for pchip, 1 for linear)
 
// We use Filon's type method, plus an asymptotic method for the correcting term toward infinity if flaginf=1.

  long double s; 
  complex<long double> x,expo,fint(0.L,0.L),Lambdat,Lambdamt,Phit,Phimt,Psit,Psimt,dummy1,dummy2,aprime,bprime;

  // computes the fourier integral at t between omegai(1) and omegai(end)
  for (unsigned long i=0; i<length-1; i++) {
      s=delta[i]*t;
      x=complex<long double>(0.L,s);
      expo=exp(x);
      switch (interp_type[i]){
        case 0:
          // pchip interpolation
	  // computes some auxiliary variables
	  Phit=Phi(x,eps);
	  Phimt=Phi(-x,eps);
	  Psit=Psi(x,eps);
	  Psimt=Psi(-x,eps);      
	  //complex<double> dummy1 = exp(complex<double>(0.,omegai[i+1]*t) ) * (fi[i]*Phimt - d[i] * delta[i] * Psimt);
	  //complex<double> dummy2 = exp(complex<double>(0.,omegai[i]*t) ) * (fi[i+1]*Phit + d[i+1] * delta[i] * Psit);
	  //fint += delta[i]*(dummy1+dummy2);
	  dummy1 = delta[i]* (-expo * (complex<long double>)d[i]*Psimt +
	  	(complex<long double>)d[i+1]* Psit);
	  dummy2 = (complex<long double>)fi[i+1]*Phit + expo * (complex<long double>)fi[i]*Phimt;
	  fint += exp(complex<long double>(0.L,omegai[i]*t) )*delta[i]*(dummy1+dummy2);
	  break;
	  
        case 1:
          // linear interpolation
	  Lambdat=Lambda(x,eps);
	  Lambdamt=Lambda(-x,eps);
	  dummy1 = (complex<long double>)fi[i+1]*Lambdat + expo * (complex<long double>)fi[i]*Lambdamt;
	  fint += exp(complex<long double>(0.L,omegai[i]*t) )*delta[i]*dummy1;
	  break;
	  
        /*case 2:
          // exponential interpolation
	  break;
	*/
      }
    }
  //cout << std::real(d[length-1]) << "  " << std::imag(d[length-1]) << '\n'; // to check interpolation 
  //fi(1:n-1).*Phimt - d(1:n-1).*delta.*Psimt)

  // add the correcting term for the rest of the integral
  if (flaginf==1) {
    if (omegai[length-1]>=0) {
      fint+= exp(complex<long double>(0.L,t*omegai[length-1]))*
    	  (-(complex<long double>)fi[length-1]/complex<long double>(0.L,t) - 
	  (complex<long double>)df*(1.L/(t*t)));
      //cout << t << " " << exp(complex<double>(0.,t*omegai[length-1]))*(-fi[length-1]/complex<double>(0.,t) - df*(1./(t*t))) << "\n";
    }
    else {
      fint+=exp(complex<long double>(0.L,t*omegai[0]))*
    	  ((complex<long double>)fi[0]/complex<long double>(0.L,t) + 
	  (complex<long double>)df*(1.L/(t*t)));
    }
  }
  return (complex<double>)fint;
}

 
/**************************************************
 *** 			main program		***
 ***						***
 ***************************************************/

main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Zycstoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 char Zxdipoutput2[MAXCHARFILE+10],Zydipoutput2[MAXCHARFILE+10],Zlongoutput2[MAXCHARFILE+10];
 char Zxquadoutput2[MAXCHARFILE+10],Zyquadoutput2[MAXCHARFILE+10],Zycstoutput2[MAXCHARFILE+10];
 char Wxdipoutput[MAXCHARFILE+10],Wydipoutput[MAXCHARFILE+10],Wlongoutput[MAXCHARFILE+10];
 char Wxquadoutput[MAXCHARFILE+10],Wyquadoutput[MAXCHARFILE+10],Wycstoutput[MAXCHARFILE+10];
 char Wxdipoutput2[MAXCHARFILE+10],Wydipoutput2[MAXCHARFILE+10],Wlongoutput2[MAXCHARFILE+10];
 char Wxquadoutput2[MAXCHARFILE+10],Wyquadoutput2[MAXCHARFILE+10],Wycstoutput2[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,topbot,commentoutput,dummy2;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZycst, *filZlong, *filInput;
 FILE *filZxdip2, *filZydip2, *filZxquad2, *filZyquad2, *filZycst2, *filZlong2;
 FILE *filWxdip, *filWydip, *filWxquad, *filWyquad, *filWycst, *filWlong;
 FILE *filWxdip2, *filWydip2, *filWxquad2, *filWyquad2, *filWycst2, *filWlong2;
 unsigned int N,M,dummy0; // number of upper and lower layers, then dummy parameter
 unsigned int n_input,n_added; /* number of input lines, number of individually added frequencies */
 double dzlin; // dz in linear scan
 unsigned int flag_topbotsym,typescan,nzlog; /* flag for top-bottom symmetry (1 if such a symmetry), type of frequency scan,
 		number of z per decade in log. scan */
 unsigned long kmain,ind[20],lprov,lprov2,j,imax=ULONG_MAX,nz,nzdup,nf; // some temporary indices, total number of z in the scan, number of duplicate z, number of frequencies;
 ap::template_1d_array< amp::ampf<Precision> > b,bm,thick,thickm; // position of upper and lower boundaries, and thickness of the layers
 ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,chi,fmu,rhom,taum,epsbm,chim,fmum; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<Precision> beta,gamma,dummy3; // parameters
 // The next three are for gsl adaptative integration algorithm
 gsl_integration_workspace *w;
 gsl_function F;
 double tolintabs=1.e-2,tolintrel=1.e-2; // absolute and relative error permitted for gsl adaptative integration
 double freqlin=1.e11; // frequency limit above which we switch from a log mesh to a linear mesh (default value)
 double eps=1.e-20; // precision for the calculation of Phi, Psi and Lambda
 double x,y,xlin,maxi,err,L,zminlog,zmaxlog,zminlin,zmaxlin,zadded[15],dummy1;
 double *freq,*z,*t,dif,sum,*newfreq,*inte;
 long double *omegai,*delta;
 double freqmin0=1.e-5,freqmin,freqmax,freqmax0;
 complex<double> *Zxdipfi,*Zxdipdi,*Zydipfi,*Zydipdi,*Zyquadfi,*Zyquaddi,*Zlongfi,*Zlongdi,*Zycstfi,*Zycstdi; //impedances
 complex<double> *newZxdipdi,*newZydipdi,*newZyquaddi,*newZlongdi,*newZycstdi; //impedances slopes on subintervals
 complex<double> *newZxdipfi,*newZydipfi,*newZyquadfi,*newZlongfi,*newZycstfi; //impedances on subintervals
 complex<double> *xx;
 double *Wakexdip,*Wakeydip,*Wakeyquad,*Wakelong,*Wakeycst; //wakes
 double *Wakexdipold,*Wakeydipold,*Wakeyquadold,*Wakelongold,*Wakeycstold; //wakes (previous loop)
 complex<double> jimag; // imaginary constant
 unsigned int *interp_type; // kind of interpolation on each interval (0 for pchip, 1 for linear)
 struct params_diff param; // input parameters for the integrand functions (gsl integration)
 size_t limit=1000,found,neval; // limit (number of intervals) for gsl integration algorithm
 time_t start,end,time1,time2; // times
 clock_t c1,c2; // clock ticks
 bool condition_int,condition_freqmin,condition_freqmax;
 //std::vector<double> newfreq;
 static const long double pi = 3.141592653589793238462643383279502884197;
 double tol=1.e9; /* absolute error permitted on the integral of the difference (in norm)
 	betweens impedances and their interpolation (default value) */
 
 
 xx=new complex<double>[5];
 for (unsigned long k=0;k<=4;k++) xx[k]=complex<double>(0.,0.);

					
 // start time
 time(&start);
 
 // allocation for memory of etas and chis
 kxmem.setbounds(0,MAXMEM);
 eta1mem.setbounds(0,MAXMEM);eta2mem.setbounds(0,MAXMEM);
 chi1mem.setbounds(0,MAXMEM);chi2mem.setbounds(0,MAXMEM);
 mem=0;

 // allocation for memory of impedances
 freqmem=new double[MAXMEM];
 Zxdipmem=new complex<double>[MAXMEM];Zydipmem=new complex<double>[MAXMEM];
 Zyquadmem=new complex<double>[MAXMEM];Zlongmem=new complex<double>[MAXMEM];
 Zycstmem=new complex<double>[MAXMEM];
 tolintabs00mem=new complex<double>[MAXMEM];tolintabs01mem=new complex<double>[MAXMEM];
 tolintabs02mem=new complex<double>[MAXMEM];tolintabs11mem=new complex<double>[MAXMEM];
 impmem=0;

 // allocation of freq, impedances and slopes for the full interpolation
 freq=new double[MAXMEM];interp_type=new unsigned int[MAXMEM];inte=new double[MAXMEM];
 Zxdipfi=new complex<double>[MAXMEM];Zydipfi=new complex<double>[MAXMEM];
 Zyquadfi=new complex<double>[MAXMEM];Zlongfi=new complex<double>[MAXMEM];
 Zycstfi=new complex<double>[MAXMEM];
 Zxdipdi=new complex<double>[MAXMEM];Zydipdi=new complex<double>[MAXMEM];
 Zyquaddi=new complex<double>[MAXMEM];Zlongdi=new complex<double>[MAXMEM];
 Zycstdi=new complex<double>[MAXMEM];
 // allocation of freq, impedances and slopes for interpolation on a subinterval
 newfreq=new double[5];
 newZxdipfi=new complex<double>[5];newZydipfi=new complex<double>[5];
 newZyquadfi=new complex<double>[5];newZlongfi=new complex<double>[5];
 newZycstfi=new complex<double>[5];
 newZxdipdi=new complex<double>[5];newZydipdi=new complex<double>[5];
 newZyquaddi=new complex<double>[5];newZlongdi=new complex<double>[5];
 newZycstdi=new complex<double>[5];

 // default values of the parameters (in case)
 flag_topbotsym=1;
 typescan=0;
 zminlog=-2;zmaxlog=6;nzlog=10;n_added=0;
 zminlin=1.e-2;zmaxlin=1;dzlin=1.e-2;nz=0;
 N=2;M=2;L=1.;gamma="479.6";
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   /* next line is when the input file comes from windows or else and has some ^M characters in the
   end of each line */
   //data[n_input+1]=data[n_input+1].substr(0,data[n_input+1].length()-1);
   //cout << data[n_input+1] << '\n';
   n_input++;
 }
 n_input--;
 //printf("n_input: %d\n",n_input);
 // identify each argument
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Machine",dummy0,dummy1,machine,dummy3,2);
   read_input(data[i],"Relativistic Gamma",dummy0,dummy1,dummy2,gamma,3);
   read_input(data[i],"Impedance Length in m",dummy0,L,dummy2,dummy3,1);
   read_input(data[i],"Number of upper layers",N,dummy1,dummy2,dummy3,0);
   //printf("yoyo %d %s %s %s %13.8e %d\n",i,data[i].c_str(),machine.c_str(),gamma.toDec().c_str(),L,N);
   read_input(data[i],"Number of lower layers",M,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Top bottom symmetry (yes or no)",dummy0,dummy1,topbot,dummy3,2);
   read_input(data[i],"exponent (10^) of zmin (in m) of the logarithmic sampling",dummy0,zminlog,dummy2,dummy3,1);
   read_input(data[i],"exponent (10^) of zmax (in m) of the logarithmic sampling",dummy0,zmaxlog,dummy2,dummy3,1);
   read_input(data[i],"linear (1) or logarithmic (0) or both (2) scan in z for the wake",typescan,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Number of points per decade for the logarithmic sampling",nzlog,dummy1,dummy2,dummy3,0);
   read_input(data[i],"zmin in m of the linear sampling",dummy0,zminlin,dummy2,dummy3,1);
   read_input(data[i],"zmax in m of the linear sampling",dummy0,zmaxlin,dummy2,dummy3,1);
   read_input(data[i],"sampling distance in m for the linear sampling",dummy0,dzlin,dummy2,dummy3,1);
   read_input(data[i],"factor weighting the longitudinal impedance error",dummy0,factlong,dummy2,dummy3,1);
   read_input(data[i],"tolerance (in wake units) to achieve",dummy0,tol,dummy2,dummy3,1);
   read_input(data[i],"frequency above which the mesh bisecting is linear",dummy0,freqlin,dummy2,dummy3,1);
   read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
   if (data[i].find("added z") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos){
       n_added=1;zadded[n_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
       while (zadded[n_added] != 0){
         n_added++;zadded[n_added]=strtod(endline,&endline);
       }
       n_added--;
     }
   }
 }
 
 tol/=(2.*(double)pi); // from angular frequancy to frequency units
 
 //printf("%s %s %d %d \n",machine,topbot,N,flag_topbotsym);
 //printf("%13.8e %13.8e %d\n",zadded[1],zadded[2],n_added);
 //printf("%13.8e %13.8e %d %ld\n",zminlog,zmaxlog,typescan,nzlog);
 //printf("%13.8e %13.8e %13.8e\n",freqlin,tol,factlong);

 // flag for top bottom symmetry (1 if there is such a symmetry)
 flag_topbotsym= ((strcmp(topbot.c_str(),"yes")==0 || strcmp(topbot.c_str(),"y")==0) || strcmp(topbot.c_str(),"1")==0);
 
 b.setbounds(1,N+1);thick.setbounds(1,N+1);rho.setbounds(1,N+1);tau.setbounds(1,N+1);
 epsb.setbounds(1,N+1);chi.setbounds(1,N+1);fmu.setbounds(1,N+1);
 if (flag_topbotsym) M=N;
 bm.setbounds(1,M+1);thickm.setbounds(1,M+1);rhom.setbounds(1,M+1);taum.setbounds(1,M+1);
 epsbm.setbounds(1,M+1);chim.setbounds(1,M+1);fmum.setbounds(1,M+1);

 // default values of the layers properties (in case)
 rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;b(2)="Infinity";

 // find inner half gap(s) of the chamber
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner half gap in mm",dummy0,dummy1,dummy2,b(1),3);
   if (flag_topbotsym) bm(1)=b(1);
   else {
     read_input(data[i],"Layer -1 inner half gap in mm",dummy0,dummy1,dummy2,bm(1),3);
   }
 }
 bm(1)=-bm(1);
 //printf("%s %s %s \n",b(1).toDec().c_str(),bm(1).toDec().c_str(),gamma.toDec().c_str());
 
 // find all the layers properties
 for (unsigned int i=1; i<=n_input; i++) {
   // upper layers
   for (unsigned int p=1; p<=N; p++){
     read_input_layer(data[i],"DC resistivity",p,rho(p+1));
     read_input_layer(data[i],"relaxation time for resistivity (ps)",p,tau(p+1));
     read_input_layer(data[i],"real part of dielectric constant",p,epsb(p+1));
     read_input_layer(data[i],"magnetic susceptibility",p,chi(p+1));
     read_input_layer(data[i],"relaxation frequency of permeability",p,fmu(p+1));
     read_input_layer(data[i],"thickness in mm",p,thick(p));
   }
   // lower layers
   if (flag_topbotsym) {
     for (unsigned int p=1; p<=M; p++){
       rhom(p+1)=rho(p+1);
       taum(p+1)=tau(p+1);
       epsbm(p+1)=epsb(p+1);
       chim(p+1)=chi(p+1);
       fmum(p+1)=fmu(p+1);
       thickm(p)=thick(p);
     }
   }
   else {
     for (unsigned int p=1; p<=M; p++){
       read_input_layer(data[i],"DC resistivity",-p,rhom(p+1));
       read_input_layer(data[i],"relaxation time for resistivity (ps)",-p,taum(p+1));
       read_input_layer(data[i],"real part of dielectric constant",-p,epsbm(p+1));
       read_input_layer(data[i],"magnetic susceptibility",-p,chim(p+1));
       read_input_layer(data[i],"relaxation frequency of permeability",-p,fmum(p+1));
       read_input_layer(data[i],"thickness in mm",-p,thickm(p));
     }
   }
 }

 // units conversion to SI (tau was in ps, b in mm and fmu in MHz)
 b(1)=b(1)*1e-3;bm(1)=bm(1)*1e-3;
 for (unsigned int p=2; p<=N+1; p++){
   tau(p)=tau(p)*1e-12;
   b(p)=thick(p-1)*1e-3+b(p-1);
   fmu(p)=fmu(p)*1e6;
   //printf("%s %s %s %s %s %s \n",rho(p).toDec().c_str(),tau(p).toDec().c_str(),epsb(p).toDec().c_str(),
   //	chi(p).toDec().c_str(),fmu(p).toDec().c_str(),b(p).toDec().c_str());
 }
 for (unsigned int p=2; p<=M+1; p++){
   taum(p)=taum(p)*1e-12;
   bm(p)=bm(p-1)-thickm(p-1)*1e-3;
   fmum(p)=fmum(p)*1e6;
   //printf("%s %s %s %s %s %s \n",rhom(p).toDec().c_str(),taum(p).toDec().c_str(),epsbm(p).toDec().c_str(),
   //	chim(p).toDec().c_str(),fmum(p).toDec().c_str(),bm(p).toDec().c_str());
 }
  
 // relativistic velocity factor beta
 beta=amp::sqrt(1-1/amp::sqr(gamma));
 
 // construct the z scan
 // first estimation of the number of z (for memory allocation)
 switch(typescan) {
   case 0:
     nz=(int)ceil((zmaxlog-zminlog)*(double)nzlog)+1+n_added;
     break;
   case 1:
     nz=(int)ceil((zmaxlin-zminlin)/dzlin)+1+n_added;
     break;
   case 2:
     nz=(int)ceil((zmaxlog-zminlog)*(double)nzlog)+1+(int)ceil((zmaxlin-zminlin)/dzlin)+1+n_added;
     break;
   }
 z=new double[nz];
 
 // constructs unsorted version of the array z
 nz=0;z[0]=-1.;

 if (typescan==1) {
   do {
     z[nz]=zminlin+(double)nz * dzlin;
     nz++;
     } while(z[nz-1]<zmaxlin);
   }
 else {
   do {
     z[nz]=pow(10.,zminlog+(double)nz/(double)nzlog);
     nz++;
     } while(z[nz-1]<pow(10.,zmaxlog));
   if (typescan==2) {
     j=0;
     do {
       z[nz]=zminlin+(double)(j) * dzlin;
       nz++;j++;
       } while(z[nz-1]<zmaxlin);
     }
   }
 for (unsigned int i=1;i<=n_added; i++) {
   z[nz]=zadded[i];
   nz++;
   }
 // nz is now the real number of z
 
 // sort the distances z
 gsl_sort(z, 1, nz);
 
 // remove duplicate distances
 nzdup=0;
 for (long i=nz-2;i>=0; i--) {
   if ( (z[i+1]-z[i]) ==0 ) {
       z[i]=DBL_MAX;
       nzdup++;
     }
   }
 gsl_sort(z, 1, nz);
 nz=nz-nzdup;
 // computes times
 t=new double[nz];
 for (unsigned long i=0;i<nz;i++) t[i]=z[i]/double(amp::ampf<Precision>(beta*C).toDouble());
 //printf("%d\n",nz);
 //for (unsigned long i=0;i<nz;i++) printf("%d %13.8e\n",i,z[i]);
 Wakexdip=new double[nz];Wakeydip=new double[nz];
 Wakeyquad=new double[nz];Wakeycst=new double[nz];
 Wakelong=new double[nz];
 Wakexdipold=new double[nz];Wakeydipold=new double[nz];
 Wakeyquadold=new double[nz];Wakeycstold=new double[nz];
 Wakelongold=new double[nz];
 
 
 // set the output files names
 //sprintf(output,"W%s_%dlayersup_%dlayersdown%.2lfmm%s",machine,N,M,
 //	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput);
 sprintf(output,"W%s_%dlayersup_%dlayersdown%.2lfmm%s",machine.c_str(),N,M,
	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput.c_str());
 sprintf(Zxdipoutput,"Zxdip%s.dat",output);
 sprintf(Zydipoutput,"Zydip%s.dat",output);
 sprintf(Zxquadoutput,"Zxquad%s.dat",output);
 sprintf(Zyquadoutput,"Zyquad%s.dat",output);
 sprintf(Zlongoutput,"Zlong%s.dat",output);
 sprintf(Zycstoutput,"Zycst%s.dat",output);
 sprintf(Zxdipoutput2,"Zxdip%s_precise.dat",output);
 sprintf(Zydipoutput2,"Zydip%s_precise.dat",output);
 sprintf(Zxquadoutput2,"Zxquad%s_precise.dat",output);
 sprintf(Zyquadoutput2,"Zyquad%s_precise.dat",output);
 sprintf(Zlongoutput2,"Zlong%s_precise.dat",output);
 sprintf(Zycstoutput2,"Zycst%s_precise.dat",output);
 sprintf(Input,"InputData%s.dat",output);
 sprintf(Wxdipoutput,"Wxdip%s.dat",output);
 sprintf(Wydipoutput,"Wydip%s.dat",output);
 sprintf(Wxquadoutput,"Wxquad%s.dat",output);
 sprintf(Wyquadoutput,"Wyquad%s.dat",output);
 sprintf(Wlongoutput,"Wlong%s.dat",output);
 sprintf(Wycstoutput,"Wycst%s.dat",output);
 sprintf(Wxdipoutput2,"Wxdip%s_precise.dat",output);
 sprintf(Wydipoutput2,"Wydip%s_precise.dat",output);
 sprintf(Wxquadoutput2,"Wxquad%s_precise.dat",output);
 sprintf(Wyquadoutput2,"Wyquad%s_precise.dat",output);
 sprintf(Wlongoutput2,"Wlong%s_precise.dat",output);
 sprintf(Wycstoutput2,"Wycst%s_precise.dat",output);
 
 // writes the InputData file (copy of input file)
 std::ofstream InputFile (Input);
 for (unsigned int i=1; i<=n_input; i++) {
   InputFile << data[i] << '\n';
 }
 InputFile.close();

 // open the output files and write the first line (column description)
 filZxdip=fopen(Zxdipoutput,"w");
 filZydip=fopen(Zydipoutput,"w");
 filZxquad=fopen(Zxquadoutput,"w");
 filZyquad=fopen(Zyquadoutput,"w");
 filZlong=fopen(Zlongoutput,"w");
 filZycst=fopen(Zycstoutput,"w");
 fprintf(filZxdip,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 fprintf(filZycst,"Frequency [Hz]\tRe(Zycst) [Ohm]\tIm(Zycst) [Ohm]\n");
 filZxdip2=fopen(Zxdipoutput2,"w");
 filZydip2=fopen(Zydipoutput2,"w");
 filZxquad2=fopen(Zxquadoutput2,"w");
 filZyquad2=fopen(Zyquadoutput2,"w");
 filZlong2=fopen(Zlongoutput2,"w");
 filZycst2=fopen(Zycstoutput2,"w");
 fprintf(filZxdip2,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip2,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad2,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad2,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong2,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 fprintf(filZycst2,"Frequency [Hz]\tRe(Zycst) [Ohm]\tIm(Zycst) [Ohm]\n");
 filWxdip=fopen(Wxdipoutput,"w");
 filWydip=fopen(Wydipoutput,"w");
 filWxquad=fopen(Wxquadoutput,"w");
 filWyquad=fopen(Wyquadoutput,"w");
 filWlong=fopen(Wlongoutput,"w");
 filWycst=fopen(Wycstoutput,"w");
 fprintf(filWxdip,"Distance [m]\tWake x dip [V/(C.m)]\n");
 fprintf(filWydip,"Distance [m]\tWake y dip [V/(C.m)]\n");
 fprintf(filWxquad,"Distance [m]\tWake x quad [V/(C.m)]\n");
 fprintf(filWyquad,"Distance [m]\tWake y quad [V/(C.m)]\n");
 fprintf(filWlong,"Distance [m]\tWake long [V/C]\n");
 fprintf(filWycst,"Distance [m]\tWake y cst [V/(C.m)]\n");
 filWxdip2=fopen(Wxdipoutput2,"w");
 filWydip2=fopen(Wydipoutput2,"w");
 filWxquad2=fopen(Wxquadoutput2,"w");
 filWyquad2=fopen(Wyquadoutput2,"w");
 filWlong2=fopen(Wlongoutput2,"w");
 filWycst2=fopen(Wycstoutput2,"w");
 fprintf(filWxdip2,"Distance [m]\tWake x dip [V/(C.m)]\n");
 fprintf(filWydip2,"Distance [m]\tWake y dip [V/(C.m)]\n");
 fprintf(filWxquad2,"Distance [m]\tWake x quad [V/(C.m)]\n");
 fprintf(filWyquad2,"Distance [m]\tWake y quad [V/(C.m)]\n");
 fprintf(filWlong2,"Distance [m]\tWake long [V/C]\n");
 fprintf(filWycst2,"Distance [m]\tWake y cst [V/(C.m)]\n");
 
 
 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 // deactivate gsl errors
 gsl_set_error_handler_off();
 
 // some parameters initialization
 param.M=M;
 param.N=N;
 param.rho=rho;
 param.rhom=rhom;
 param.tau=tau;
 param.taum=taum;
 param.epsb=epsb;
 param.epsbm=epsbm;
 param.chi=chi;
 param.chim=chim;
 param.fmu=fmu;
 param.fmum=fmum;
 param.b=b;
 param.bm=bm;
 param.beta=beta;
 param.gamma=gamma;
 param.flag_topbotsym=flag_topbotsym;
 param.L=L;
 param.x=xx;
 
 // first guess for freqmin and freqmax (minimum and maximum frequencies of the interpolation)
 freqmin=freqmin0;
 freqmax0=10.*max(double(amp::ampf<Precision>(beta*gamma*C/(b(1)*amp::twopi<Precision>())).toDouble()),
 	double(amp::ampf<Precision>(beta*gamma*C/(bm(1)*amp::twopi<Precision>())).toDouble()));
 freqmax=freqmax0;
	
 condition_freqmin=true;
 condition_freqmax=true;
 srand(time(NULL));
 
 nf=2; // begins with two frequencies
 freq[0]=freqmin;freq[1]=freqmax;
 condition_int=true;
 
 kmain=0;

 // loop to get accurate enough mesh of the impedances
 while ( condition_int ) {

   printf("Number of frequencies: %ld, iteration nb %ld\n",nf,kmain);
   cout.flush();

   if (kmain==0) {
     for (unsigned long i=0; i<=nf-1; i++) {
       //time(&time1);
       impedance(Zxdipfi[i], Zydipfi[i], Zyquadfi[i], Zlongfi[i], Zycstfi[i], N, M, rho, tau, epsb, chi, fmu, b,
		   rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, freq[i]);
       /*time(&time2);
       dif=difftime(time2,time1);
       printf("Elapsed time during calculation: %.5lf seconds\n",dif);
       cout.flush();*/
     }
     pchipslope(Zxdipdi,freq,Zxdipfi,nf);pchipslope(Zydipdi,freq,Zydipfi,nf);
     pchipslope(Zyquaddi,freq,Zyquadfi,nf);pchipslope(Zycstdi,freq,Zycstfi,nf);
     pchipslope(Zlongdi,freq,Zlongfi,nf);

     // writing the rest of the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.interp_type=0; // pchip interpolation chosen for this step
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zydipfi=Zydipfi;param.Zydipdi=Zydipdi;
     param.Zyquadfi=Zyquadfi;param.Zyquaddi=Zyquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;
     param.Zycstfi=Zycstfi;param.Zycstdi=Zycstdi;

     // parameters for gsl integration
     F.params=&param;

     condition_int=false;
     //newfreq.clear();

     sum=0.;
      
     // make an adaptative integration to get a first trial mesh
     F.function=&integrand_diff; // function to integrate with gsl
     //F.function=&integrand_diff_freq; // function to integrate with gsl
     //time(&time1);
     gsl_integration_qag(&F, log(freq[0]), log(freq[1]), tolintabs, tolintrel, limit,1, w, &sum, &err);
     //gsl_integration_qag(&F, freq[0], freq[1], tolintabs, tolintrel, limit,1, w, &sum, &err);
     // replace freq by all the frequencies in the impedance memory
     //for (unsigned long i=0; i<=impmem-1; i++) newfreq.push_back(freqmem[i]);
     for (unsigned long i=0; i<=impmem-1; i++) {
       if (i<impmem-1) interp_type[i]=0; // pchip is chosen
       freq[i]=freqmem[i];Zxdipfi[i]=Zxdipmem[i];Zydipfi[i]=Zydipmem[i];
       Zyquadfi[i]=Zyquadmem[i];Zycstfi[i]=Zycstmem[i];Zlongfi[i]=Zlongmem[i];
     }
     nf=impmem;
     condition_int=true;
   }
   else if (kmain==1) {
     
     // compute the pchip slopes
     pchipslope(Zxdipdi,freq,Zxdipfi,nf);pchipslope(Zydipdi,freq,Zydipfi,nf);
     pchipslope(Zyquaddi,freq,Zyquadfi,nf);pchipslope(Zycstdi,freq,Zycstfi,nf);
     pchipslope(Zlongdi,freq,Zlongfi,nf);

     // rewrite the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zydipfi=Zydipfi;param.Zydipdi=Zydipdi;
     param.Zyquadfi=Zyquadfi;param.Zyquaddi=Zyquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;
     param.Zycstfi=Zycstfi;param.Zycstdi=Zycstdi;
     
     err=DBL_MIN;sum=0.;
     for (unsigned long i=0; i<nf-1; i++) {
       if (freq[i]>=freqlin) y=(freq[i]+freq[i+1])/2.;	 
       else y=(log(freq[i])+log(freq[i+1]))/2.;
       // first try pchip
       param.interp_type=0;F.params=&param;interp_type[i]=0;
       if (freq[i]>=freqlin) inte[i]=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else inte[i]=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       // then try linear
       param.interp_type=1;F.params=&param;
       if (freq[i]>=freqlin) xlin=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else xlin=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       if (xlin<inte[i]) {
         inte[i]=xlin;interp_type[i]=1;
       }
       if (freq[i]>=freqlin) inte[i]*=(freq[i+1]-freq[i])/6.;
       else inte[i]*=(log(freq[i+1])-log(freq[i]))/6.;
       sum+=inte[i];
       if (inte[i]>err) {
	 err=inte[i];imax=i; // interval with largest error
       }
     }
     cout << "err=" << err << ", imax=" << imax << ", f[imax]=" << freq[imax] << ", f[imax+1]=" << freq[imax+1] << ", sum=" << sum << "\n";
     condition_int=(sum>tol);
     if ( (!condition_int)&&!(sum<=tol) ) condition_int=true; // case when sum=nan -> will still try another loop
     if (condition_int) nf++;
     /*if (condition_int) {
       newfreq.push_back(freq[0]);
       for (unsigned long i=0; i<nf-1; i++) {
         if (i==imax) {
           // refine the worst interval found above with the points in the impedance memory
	   lprov=locate(freqmem,freq[i],impmem-1);lprov2=locate(freqmem,freq[i+1],impmem-1);
	   if ( (lprov>0)&&(freq[i]==freqmem[lprov-1]) ) lprov=lprov-1;
	   if ( (lprov2>0)&&(freq[i+1]==freqmem[lprov2-1]) ) lprov2=lprov2-1;
           for (unsigned long k=lprov+1; k<lprov2; k++) {
	     //cout << freqmem[k] << "\n";
	     newfreq.push_back(freqmem[k]);
	   }
	 }
	 newfreq.push_back(freq[i+1]);
       }
     }*/
   } else {
   
     // first subtract the interval imax from the integral (in 'sum') previously computed
     sum-=inte[imax];

     // bisect the previously found subinterval imax (the one with maximum error)

     // construct the five points array to compute the new pchip slopes on the 3 points freq[imax],
     // freq[imax+1] and the new frequency (half way in log between the two)
     for (unsigned long i=0; i<5; i++) {
       if (i!=2) {
         newfreq[i]=freq[imax-1+i-(i>2)];newZxdipfi[i]=Zxdipfi[imax-1+i-(i>2)];
	 newZydipfi[i]=Zydipfi[imax-1+i-(i>2)];newZyquadfi[i]=Zyquadfi[imax-1+i-(i>2)];
	 newZycstfi[i]=Zycstfi[imax-1+i-(i>2)];newZlongfi[i]=Zlongfi[imax-1+i-(i>2)];
       } else {
         if (freq[imax]>=freqlin) newfreq[i]=(freq[imax]+freq[imax+1])/2.;
         else newfreq[i]=exp((log(freq[imax])+log(freq[imax+1]))/2.);
	 impedance(newZxdipfi[i], newZydipfi[i], newZyquadfi[i], newZlongfi[i], newZycstfi[i], N, M, rho, tau, epsb, chi, fmu, b,
		   rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, newfreq[i]);
       }
     }
     // pchip slopes on those five points
     pchipslope(newZxdipdi,newfreq,newZxdipfi,5);pchipslope(newZydipdi,newfreq,newZydipfi,5);
     pchipslope(newZyquaddi,newfreq,newZyquadfi,5);pchipslope(newZycstdi,newfreq,newZycstfi,5);
     pchipslope(newZlongdi,newfreq,newZlongfi,5);
     // construct the new mesh with impedances and slopes
     for (unsigned long i=nf-1; i>imax+1; i--) {
       if (i<(nf-1)) {
         interp_type[i]=interp_type[i-1];inte[i]=inte[i-1];
       }
       freq[i]=freq[i-1]; 
       Zxdipfi[i]=Zxdipfi[i-1]; Zydipfi[i]=Zydipfi[i-1];
       Zyquadfi[i]=Zyquadfi[i-1]; Zycstfi[i]=Zycstfi[i-1]; Zlongfi[i]=Zlongfi[i-1];
       Zxdipdi[i]=Zxdipdi[i-1]; Zydipdi[i]=Zydipdi[i-1];
       Zyquaddi[i]=Zyquaddi[i-1]; Zycstdi[i]=Zycstdi[i-1]; Zlongdi[i]=Zlongdi[i-1];
     }
     freq[imax+1]=newfreq[2];
     Zxdipfi[imax+1]=newZxdipfi[2]; Zydipfi[imax+1]=newZydipfi[2];
     Zyquadfi[imax+1]=newZyquadfi[2]; Zycstfi[imax+1]=newZycstfi[2]; Zlongfi[imax+1]=newZlongfi[2];
     interp_type[imax]=0; interp_type[imax+1]=0;// choose by default pchip on the two subintervals
     for (unsigned long i=imax; i<=imax+2; i++) {
       Zxdipdi[i]=newZxdipdi[1+i-imax]; Zydipdi[i]=newZydipdi[1+i-imax];
       Zyquaddi[i]=newZyquaddi[1+i-imax]; Zycstdi[i]=newZycstdi[1+i-imax]; Zlongdi[i]=newZlongdi[1+i-imax];
     }
     
     // rewrite the input parameters structure
     param.freqi=freq;param.nf=nf;
     param.Zxdipfi=Zxdipfi;param.Zxdipdi=Zxdipdi;
     param.Zydipfi=Zydipfi;param.Zydipdi=Zydipdi;
     param.Zyquadfi=Zyquadfi;param.Zyquaddi=Zyquaddi;
     param.Zlongfi=Zlongfi;param.Zlongdi=Zlongdi;
     param.Zycstfi=Zycstfi;param.Zycstdi=Zycstdi;
     
     // compute the integral (with Simpson's rule) on the two subintervals created
     for (unsigned long i=imax;i<=imax+1;i++) {
       if (freq[i]>=freqlin) y=(freq[i]+freq[i+1])/2.;	 
       else y=(log(freq[i])+log(freq[i+1]))/2.;
       // first try pchip
       param.interp_type=0;F.params=&param;interp_type[i]=0;
       if (freq[i]>=freqlin) inte[i]=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else inte[i]=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       // then try linear
       param.interp_type=1;F.params=&param;
       if (freq[i]>=freqlin) xlin=integrand_diff_freq(freq[i],F.params)+4.*integrand_diff_freq(y,F.params)+integrand_diff_freq(freq[i+1],F.params);
       else xlin=integrand_diff(log(freq[i]),F.params)+4.*integrand_diff(y,F.params)+integrand_diff(log(freq[i+1]),F.params);
       if (xlin<inte[i]) {
         inte[i]=xlin;interp_type[i]=1;
       }
       if (freq[i]>=freqlin) inte[i]*=(freq[i+1]-freq[i])/6.;
       else inte[i]*=(log(freq[i+1])-log(freq[i]))/6.;
       sum+=inte[i];
     }
     // find largest error 'err' and its index 'imax'
     err=DBL_MIN;
     for (unsigned long i=0;i<nf-1;i++) {
       if (inte[i]>err) {
         err=inte[i];imax=i;
       }
     }
     cout << "err=" << err << ", imax=" << imax << ", f[imax]=" << freq[imax] << ", f[imax+1]=" << freq[imax+1] << ", sum=" << sum << "\n";

     condition_int=(sum>tol); // condition to continue
     if ( (!condition_int)&&!(sum<=tol) ) condition_int=true; // case when sum=nan -> will still try another loop
     if (condition_int) nf++;
   
   }
   
   
   kmain++;

 }
 
 /*cout << "Final mesh and interpolation type: " << "\n";
 for (unsigned long i=0;i<nf-1;i++) {
   cout << freq[i] << " " << interp_type[i] << "\n";
 }*/

 // for the slopes, convert frequencies to angular frequencies
 for (unsigned long i=0; i<=nf-1; i++) {
   Zxdipdi[i]=Zxdipdi[i]/(2.*(double)pi);Zydipdi[i]=Zydipdi[i]/(2.*(double)pi);
   Zyquaddi[i]=Zyquaddi[i]/(2.*(double)pi);Zycstdi[i]=Zycstdi[i]/(2.*(double)pi);Zlongdi[i]=Zlongdi[i]/(2.*(double)pi);
 }
 
 // loop to get low enough minimum frequency
 while (condition_freqmin) {
   
   printf("freq. min= %13.8e\n",freqmin);cout.flush();
   
   // computes omega and delta between successive omegas
   omegai=new long double[nf];delta=new long double[nf-1];
   for (unsigned long i=0; i<=nf-1; i++) omegai[i]=2.L*pi*(long double)freq[i];
   for (unsigned long i=0; i<nf-1; i++) delta[i]=omegai[i+1]-omegai[i];
   
   printf("Wake computation\n");cout.flush();
   // computes the wakes
   if (freqmin==freqmin0) {
     for (unsigned long i=0; i<=nz-1; i++) {
       Wakexdip[i]=std::imag(fourier_integral_inf(Zxdipfi,Zxdipdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakeydip[i]=std::imag(fourier_integral_inf(Zydipfi,Zydipdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakeyquad[i]=std::imag(fourier_integral_inf(Zyquadfi,Zyquaddi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakelong[i]=std::real(fourier_integral_inf(Zlongfi,Zlongdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
       Wakeycst[i]=std::imag(fourier_integral_inf(Zycstfi,Zycstdi,0.,(long double)t[i],omegai,delta,nf,eps,interp_type,1))/pi;
     }
   } else {
     for (unsigned long i=0; i<=nz-1; i++) {
       Wakexdip[i]=Wakexdipold[i]+std::imag(fourier_integral_inf(newZxdipfi,newZxdipdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zxdipfi[1],&Zxdipdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakeydip[i]=Wakeydipold[i]+std::imag(fourier_integral_inf(newZydipfi,newZydipdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zydipfi[1],&Zydipdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakeyquad[i]=Wakeyquadold[i]+std::imag(fourier_integral_inf(newZyquadfi,newZyquaddi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zyquadfi[1],&Zyquaddi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakelong[i]=Wakelongold[i]+std::real(fourier_integral_inf(newZlongfi,newZlongdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::real(fourier_integral_inf(&Zlongfi[1],&Zlongdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
       Wakeycst[i]=Wakeycstold[i]+std::imag(fourier_integral_inf(newZycstfi,newZycstdi,0.,(long double)t[i],omegai,delta,
     		  3,eps,interp_type,0))/pi - std::imag(fourier_integral_inf(&Zycstfi[1],&Zycstdi[1],0.,
		  (long double)t[i],&omegai[1],&delta[1],2,eps,&interp_type[1],0))/pi;
     }
     for (unsigned int i=0;i<=1;i++) {
       Zxdipdi[i]=newZxdipdi[i];Zydipdi[i]=newZydipdi[i];
       Zyquaddi[i]=newZyquaddi[i];Zycstdi[i]=newZycstdi[i];
       Zlongdi[i]=newZlongdi[i];
     }
   }
   
   if (freqmin!=freqmin0) {
     err=DBL_MIN;
     for (unsigned long i=0; i<=nz-1; i++) {
       err=max(max(abs(Wakexdip[i]-Wakexdipold[i]),abs(Wakeydip[i]-Wakeydipold[i])),err);
       err=max(abs(Wakeyquad[i]-Wakeyquadold[i]),err);
       err=max(max(factlong*abs(Wakelong[i]-Wakelongold[i]),abs(Wakeycst[i]-Wakeycstold[i])),err);
     }
     cout << "Max. error between two last minimum frequencies chosen : " << err << "\n";cout.flush();
     condition_freqmin=(err>=tol);
   }
   
   if ( condition_freqmin ) {
     freqmin/=10.;nf++;
     for (unsigned long i=nf-1; i>0; i--) {
       if (i<(nf-1)) interp_type[i]=interp_type[i-1];
       freq[i]=freq[i-1];
       Zxdipfi[i]=Zxdipfi[i-1]; Zydipfi[i]=Zydipfi[i-1];
       Zyquadfi[i]=Zyquadfi[i-1]; Zycstfi[i]=Zycstfi[i-1]; Zlongfi[i]=Zlongfi[i-1];
       Zxdipdi[i]=Zxdipdi[i-1]; Zydipdi[i]=Zydipdi[i-1];
       Zyquaddi[i]=Zyquaddi[i-1]; Zycstdi[i]=Zycstdi[i-1]; Zlongdi[i]=Zlongdi[i-1];
     }
     freq[0]=freqmin;interp_type[0]=0; // use pchip only
     impedance(Zxdipfi[0], Zydipfi[0], Zyquadfi[0], Zlongfi[0], Zycstfi[0], N, M, rho, tau, epsb, chi, fmu, b,
		   rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, freqmin);
     for (unsigned int i=0;i<=2;i++) {
       newfreq[i]=freq[i];newZxdipfi[i]=Zxdipfi[i];
       newZydipfi[i]=Zydipfi[i];newZyquadfi[i]=Zyquadfi[i];
       newZycstfi[i]=Zycstfi[i];newZlongfi[i]=Zlongfi[i];
     }
     pchipslope(newZxdipdi,newfreq,newZxdipfi,3);newZxdipdi[2]=Zxdipdi[2];// last slope remains unchanged
     pchipslope(newZydipdi,newfreq,newZydipfi,3);newZydipdi[2]=Zydipdi[2];// last slope remains unchanged
     pchipslope(newZyquaddi,newfreq,newZyquadfi,3);newZyquaddi[2]=Zyquaddi[2];// last slope remains unchanged
     pchipslope(newZycstdi,newfreq,newZycstfi,3);newZycstdi[2]=Zycstdi[2];// last slope remains unchanged
     pchipslope(newZlongdi,newfreq,newZlongfi,3);newZlongdi[2]=Zlongdi[2];// last slope remains unchanged
     // convert to angular frequencies
     for (unsigned int i=0;i<=1;i++) {
       newZxdipdi[i]=newZxdipdi[i]/(2.*(double)pi);
       newZydipdi[i]=newZydipdi[i]/(2.*(double)pi);
       newZyquaddi[i]=newZyquaddi[i]/(2.*(double)pi);
       newZlongdi[i]=newZlongdi[i]/(2.*(double)pi);
       newZycstdi[i]=newZycstdi[i]/(2.*(double)pi);
     }     
     delete[] omegai;delete[] delta;
   }
   
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdipold[i]=Wakexdip[i];Wakeydipold[i]=Wakeydip[i];
     Wakeyquadold[i]=Wakeyquad[i];Wakeycstold[i]=Wakeycst[i];
     Wakelongold[i]=Wakelong[i];
   }

 }

 // loop to get high enough maximum frequency
 while (condition_freqmax) {
   
   freqmax*=2.;nf++;
   printf("freq. max= %13.8e\n",freqmax);cout.flush();
   freq[nf-1]=freqmax;interp_type[nf-2]=0; // use pchip only
   impedance(Zxdipfi[nf-1], Zydipfi[nf-1], Zyquadfi[nf-1], Zlongfi[nf-1], Zycstfi[nf-1], N, M, rho, tau, epsb, chi, fmu, b,
		 rhom, taum, epsbm, chim, fmum, bm, beta, gamma, flag_topbotsym, L, freqmax);
   for (unsigned long i=0;i<=2;i++) {
     newfreq[i]=freq[nf-3+i];newZxdipfi[i]=Zxdipfi[nf-3+i];
     newZydipfi[i]=Zydipfi[nf-3+i];newZyquadfi[i]=Zyquadfi[nf-3+i];
     newZycstfi[i]=Zycstfi[nf-3+i];newZlongfi[i]=Zlongfi[nf-3+i];
   }
   pchipslope(newZxdipdi,newfreq,newZxdipfi,3);newZxdipdi[0]=Zxdipdi[nf-3]; // first slope remains unchanged
   pchipslope(newZydipdi,newfreq,newZydipfi,3);newZydipdi[0]=Zydipdi[nf-3]; // first slope remains unchanged
   pchipslope(newZyquaddi,newfreq,newZyquadfi,3);newZyquaddi[0]=Zyquaddi[nf-3]; // first slope remains unchanged
   pchipslope(newZycstdi,newfreq,newZycstfi,3);newZycstdi[0]=Zycstdi[nf-3]; // first slope remains unchanged
   pchipslope(newZlongdi,newfreq,newZlongfi,3);newZlongdi[0]=Zlongdi[nf-3]; // first slope remains unchanged
   // convert to angular frequencies
   for (unsigned int i=1;i<=2;i++) {
     newZxdipdi[i]=newZxdipdi[i]/(2.*(double)pi);
     newZydipdi[i]=newZydipdi[i]/(2.*(double)pi);
     newZyquaddi[i]=newZyquaddi[i]/(2.*(double)pi);
     newZycstdi[i]=newZycstdi[i]/(2.*(double)pi);
     newZlongdi[i]=newZlongdi[i]/(2.*(double)pi);
   }     
   
   // computes omega and delta between successive omegas
   delete[] omegai;delete[] delta;
   omegai=new long double[3];delta=new long double[2];
   for (unsigned long i=0; i<=2; i++) omegai[i]=2.L*pi*(long double)freq[nf-3+i];
   for (unsigned long i=0; i<=1; i++) delta[i]=omegai[i+1]-omegai[i];
   
   printf("Wake computation\n");cout.flush();
   // computes the wakes
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdip[i]=Wakexdipold[i]+std::imag(fourier_integral_inf(newZxdipfi,newZxdipdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zxdipfi[nf-3],&Zxdipdi[nf-3],0.,(long double)t[i],omegai,delta,
     		2,eps,&interp_type[nf-3],1))/pi;
     Wakeydip[i]=Wakeydipold[i]+std::imag(fourier_integral_inf(newZydipfi,newZydipdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zydipfi[nf-3],&Zydipdi[nf-3],0.,(long double)t[i],omegai,delta,
     		2,eps,&interp_type[nf-3],1))/pi;
     Wakeyquad[i]=Wakeyquadold[i]+std::imag(fourier_integral_inf(newZyquadfi,newZyquaddi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zyquadfi[nf-3],&Zyquaddi[nf-3],0.,(long double)t[i],omegai,delta,
     		2,eps,&interp_type[nf-3],1))/pi;
     Wakelong[i]=Wakelongold[i]+std::real(fourier_integral_inf(newZlongfi,newZlongdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::real(fourier_integral_inf(&Zlongfi[nf-3],&Zlongdi[nf-3],0.,(long double)t[i],omegai,delta,
     		2,eps,&interp_type[nf-3],1))/pi;
     Wakeycst[i]=Wakeycstold[i]+std::imag(fourier_integral_inf(newZycstfi,newZycstdi,0.,(long double)t[i],omegai,delta,
     		3,eps,&interp_type[nf-3],1))/pi - std::imag(fourier_integral_inf(&Zycstfi[nf-3],&Zycstdi[nf-3],0.,(long double)t[i],omegai,delta,
     		2,eps,&interp_type[nf-3],1))/pi;
   }
   for (unsigned long i=0;i<=1;i++) {
     Zxdipdi[nf-2+i]=newZxdipdi[i+1];Zydipdi[nf-2+i]=newZydipdi[i+1];
     Zyquaddi[nf-2+i]=newZyquaddi[i+1];Zycstdi[nf-2+i]=newZycstdi[i+1];
     Zlongdi[nf-2+i]=newZlongdi[i+1];
   }
   
   err=DBL_MIN;
   for (unsigned long i=0; i<=nz-1; i++) {
     err=max(max(abs(Wakexdip[i]-Wakexdipold[i]),abs(Wakeydip[i]-Wakeydipold[i])),err);
     err=max(abs(Wakeyquad[i]-Wakeyquadold[i]),err);
     err=max(max(factlong*abs(Wakelong[i]-Wakelongold[i]),abs(Wakeycst[i]-Wakeycstold[i])),err);
   }
   cout << "Max. error between two last maximum frequencies chosen : " << err << "\n";cout.flush();
   condition_freqmax=(err>=tol);
   
   for (unsigned long i=0; i<=nz-1; i++) {
     Wakexdipold[i]=Wakexdip[i];Wakeydipold[i]=Wakeydip[i];
     Wakeyquadold[i]=Wakeyquad[i];Wakeycstold[i]=Wakeycst[i];
     Wakelongold[i]=Wakelong[i];
   }

 }
 
 printf("Final wake computation with %ld frequencies\n",impmem);cout.flush();
 // compute the pchip slopes on the full frequency range in memory
 pchipslope(Zxdipdi,freqmem,Zxdipmem,impmem);pchipslope(Zydipdi,freqmem,Zydipmem,impmem);
 pchipslope(Zyquaddi,freqmem,Zyquadmem,impmem);pchipslope(Zycstdi,freqmem,Zycstmem,impmem);
 pchipslope(Zlongdi,freqmem,Zlongmem,impmem);
 // computes omega and delta between successive omegas, and choose interpolation type 
 // on each interval (pchip)
 delete[] omegai;delete[] delta;
 omegai=new long double[impmem];delta=new long double[impmem-1];
 for (unsigned long i=0; i<=impmem-1; i++) {
   omegai[i]=2.L*pi*(long double)freqmem[i];Zxdipdi[i]=Zxdipdi[i]/(2.*(double)pi);Zydipdi[i]=Zydipdi[i]/(2.*(double)pi);
   Zyquaddi[i]=Zyquaddi[i]/(2.*(double)pi);Zycstdi[i]=Zycstdi[i]/(2.*(double)pi);Zlongdi[i]=Zlongdi[i]/(2.*(double)pi);
 }
 for (unsigned long i=0; i<impmem-1; i++) {
   delta[i]=omegai[i+1]-omegai[i];interp_type[i]=0;
 }
 // computes the final wakes using all frequencies in memory
 for (unsigned long i=0; i<=nz-1; i++) {
   Wakexdip[i]=std::imag(fourier_integral_inf(Zxdipmem,Zxdipdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakeydip[i]=std::imag(fourier_integral_inf(Zydipmem,Zydipdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakeyquad[i]=std::imag(fourier_integral_inf(Zyquadmem,Zyquaddi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakelong[i]=std::real(fourier_integral_inf(Zlongmem,Zlongdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
   Wakeycst[i]=std::imag(fourier_integral_inf(Zycstmem,Zycstdi,0.,(long double)t[i],omegai,delta,impmem,eps,interp_type,1))/pi;
 }
 // compare with previous version on the converged mesh
 err=DBL_MIN;
 for (unsigned long i=0; i<=nz-1; i++) {
   err=max(max(abs(Wakexdip[i]-Wakexdipold[i]),abs(Wakeydip[i]-Wakeydipold[i])),err);
   err=max(abs(Wakeyquad[i]-Wakeyquadold[i]),err);
   err=max(max(abs(Wakelong[i]-Wakelongold[i]),abs(Wakeycst[i]-Wakeycstold[i])),err);
 }
 cout << "Max. error between two last frequency meshes : " << err << "\n";
  
   
   
 // writes the final impedances on the final mesh chosen
 for (unsigned long i=0; i<=nf-1; i++) {
   fprintf(filZxdip,"%13.8e %13.8e %13.8e\n",freq[i],Zxdipfi[i].real(),Zxdipfi[i].imag());
   fprintf(filZydip,"%13.8e %13.8e %13.8e\n",freq[i],Zydipfi[i].real(),Zydipfi[i].imag());
   fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],-Zxdipfi[i].real(),-Zxdipfi[i].imag());
   fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],Zyquadfi[i].real(),Zyquadfi[i].imag());
   fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],Zlongfi[i].real(),Zlongfi[i].imag());
   fprintf(filZycst,"%13.8e %13.8e %13.8e\n",freq[i],Zycstfi[i].real(),Zycstfi[i].imag());
 }

 // writes the impedances with the finest possible mesh
 for (unsigned long i=0; i<=impmem-1; i++) {
   fprintf(filZxdip2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zxdipmem[i].real(),Zxdipmem[i].imag());
   fprintf(filZydip2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zydipmem[i].real(),Zydipmem[i].imag());
   fprintf(filZxquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],-Zxdipmem[i].real(),-Zxdipmem[i].imag());
   fprintf(filZyquad2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zyquadmem[i].real(),Zyquadmem[i].imag());
   fprintf(filZlong2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zlongmem[i].real(),Zlongmem[i].imag());
   fprintf(filZycst2,"%13.8e %13.8e %13.8e\n",freqmem[i],Zycstmem[i].real(),Zycstmem[i].imag());
 }

 // writes the final wakes 
 for (unsigned long i=0; i<=nz-1; i++) {
   fprintf(filWxdip,"%13.8e %13.8e\n",z[i],Wakexdipold[i]);
   fprintf(filWydip,"%13.8e %13.8e\n",z[i],Wakeydipold[i]);
   fprintf(filWxquad,"%13.8e %13.8e\n",z[i],-Wakexdipold[i]);
   fprintf(filWyquad,"%13.8e %13.8e\n",z[i],Wakeyquadold[i]);
   fprintf(filWlong,"%13.8e %13.8e\n",z[i],Wakelongold[i]);
   fprintf(filWycst,"%13.8e %13.8e\n",z[i],Wakeycstold[i]);
 } 
 
 // writes the wakes with the finest possible mesh
 for (unsigned long i=0; i<=nz-1; i++) {
   fprintf(filWxdip2,"%13.8e %13.8e\n",z[i],Wakexdip[i]);
   fprintf(filWydip2,"%13.8e %13.8e\n",z[i],Wakeydip[i]);
   fprintf(filWxquad2,"%13.8e %13.8e\n",z[i],-Wakexdip[i]);
   fprintf(filWyquad2,"%13.8e %13.8e\n",z[i],Wakeyquad[i]);
   fprintf(filWlong2,"%13.8e %13.8e\n",z[i],Wakelong[i]);
   fprintf(filWycst2,"%13.8e %13.8e\n",z[i],Wakeycst[i]);
 } 
 
 //finalization
 gsl_integration_workspace_free(w);
 delete[] xx;
 
 delete[] freqmem;delete[] Zxdipmem;delete[] Zydipmem;
 delete[] Zyquadmem;delete[] Zycstmem;delete[] Zlongmem;
 
 delete[] freq;delete[] interp_type; delete[] inte;
 delete[] Zxdipfi;delete[] Zxdipdi;delete[] Zydipfi;delete[] Zydipdi;
 delete[] Zyquadfi;delete[] Zyquaddi;delete[] Zycstfi;delete[] Zycstdi;
 delete[] Zlongfi;delete[] Zlongdi;
 
 delete[] newfreq;
 delete[] newZxdipfi;delete[] newZxdipdi;delete[] newZydipfi;delete[] newZydipdi;
 delete[] newZyquadfi;delete[] newZyquaddi;delete[] newZycstfi;delete[] newZycstdi;
 delete[] newZlongfi;delete[] newZlongdi;
 
 delete[] omegai;delete[] delta;
 delete[] z;delete[] t;
 delete[] Wakexdip;delete[] Wakeydip;
 delete[] Wakeyquad;delete[] Wakeycst;
 delete[] Wakelong;
 delete[] Wakexdipold;delete[] Wakeydipold;
 delete[] Wakeyquadold;delete[] Wakeycstold;
 delete[] Wakelongold;

 fclose(filZxdip);
 fclose(filZydip);
 fclose(filZxquad);
 fclose(filZyquad);
 fclose(filZlong);
 fclose(filZycst);
 fclose(filZxdip2);
 fclose(filZydip2);
 fclose(filZxquad2);
 fclose(filZyquad2);
 fclose(filZlong2);
 fclose(filZycst2);
 fclose(filWxdip);
 fclose(filWydip);
 fclose(filWxquad);
 fclose(filWyquad);
 fclose(filWlong);
 fclose(filWycst);
 fclose(filWxdip2);
 fclose(filWydip2);
 fclose(filWxquad2);
 fclose(filWyquad2);
 fclose(filWlong2);
 fclose(filWycst2);
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
