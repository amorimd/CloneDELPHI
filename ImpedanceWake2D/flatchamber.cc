/*
 *  flatchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the impedance in a flat chamber (see CERN note by N. Mounet and E. Metral, 
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
Layer -2 relaxation time for resistivity (ps)
:	0
Layer -2 real part of dielectric constant:	1.43
Layer -2 magnetic susceptibility:	0.02
Layer -2 relaxation frequency of permeability (MHz):	Infinity
Layer -2 thickness in mm:	Infinity
start frequency exponent (10^) in Hz:	0
stop frequency exponent (10^) in Hz:	14
linear (1) or logarithmic (0) or both (2) frequency scan:	2
sampling frequency exponent (10^) in Hz (for linear):	8
Number of points per decade (for log):	20
when both, fmin of the refinement (in THz):	5.5
when both, fmax of the refinement (in THz):	7
when both, number of points in the refinement:	100
added frequencies [Hz]:	1e-3 1e-1
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If top-bottom symmetry is set (with "yes" or "y" or "1") the lower layers (with a
minus sign) are ignored. Also if there are more layers than indicated by the number of upper (lower) layers,
the additional one are ignored. The last layer is always assumed to go to infinity.

In output one gives six files with the impedances (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar and y constant term). Each have 3 columns :frequency, real part and
imaginary part.

 */

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <sstream>
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
#include <gsl/gsl_errno.h>


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
#define  MAXMEM		50000 // Maximum number of elements if arrays with etas and chis

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
// NB: eps1, eps1m, mu1 and mu1m are actually the inverse of the 'real' eps1, etc.

// global arrays with memory of etas and chis
ap::template_1d_array< amp::ampf<Precision> > kxmem;
ap::template_1d_array< amp::campf<Precision> > eta1mem,eta2mem,chi1mem,chi2mem;
unsigned long mem; // current number of elements of kxmem, eta1mem, etc.


//using std::complex;
  
  
/******************************************************************************
 *** locate: search a table ordered in ascending order		    ***
 ***		      (inspired from Numerical Recipes) 		    ***
 *** Effect         : Function that gives the position lprov (integer) in   ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
 *** Parameters     : table, z, n (table is indexed from 0 to n)            ***
 ******************************************************************************/

unsigned long locate (ap::template_1d_array< amp::ampf<Precision> >& table, 
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
    
    m=mold;

    kx2=amp::sqr(kx);
    kyp=csqrt(kx2+nu2(1));
    fac3=kx/(2*beta);
    
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
    else lprov=locate(kxmem,kx,mem-1);
 
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
      //pcal2=pcal;

      /*for (int k=1; k<=4; k++) {
	for (int j=1; j<=4; j++) {
          printf("%d %d : %s %s\n", k,j,mat(k,j).x.toDec().c_str(),mat(k,j).y.toDec().c_str());
	}
	printf("\n");
      }*/

      //matinv::cmatrixinverse(pcal,4,info,repi);
      matinv4(pcal); // ~50 times quicker but less accurate -> increase Precision
      /*amp::ampf<Precision> sum=0;
      for (unsigned int p=0;p<=1;p++) {
        for (unsigned int q=0;q<=3;q++) sum=amp::maximum(sum,amp::abscomplex((pcal(p,q)-pcal2(p,q))/pcal(p,q)));
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
    //std::cout << des.str() << "\n";
    //std::cout << line << "\n";
    read_input(line,des.str(),dummy0,dummy1,dummy2,param,3);
    
    return;
    
  }
    

/**************************************************
 *** 			main program		***
 ***						***
 **************************************************/

main ()

{
 
 char *endline;
 char output[MAXCHARFILE],Zxdipoutput[MAXCHARFILE+10],Zydipoutput[MAXCHARFILE+10],Zlongoutput[MAXCHARFILE+10];
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Zycstoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,topbot,commentoutput,dummy2;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZycst, *filZlong, *filInput;
 unsigned int N,M,dummy0; // number of upper and lower layers, then dummy parameter
 unsigned int m,n,n_input,n_added,nf,nfdup; /* indices of alphamn (azimuthal mode numbers), number of input lines,
 					number of individually added frequencies, total number of
					frequencies in the scan, number of duplicate frequencies */
 unsigned int flag_topbotsym,typescan,nflog,nflin; /* flag for top-bottom symmetry (1 if such a symmetry), type of frequency scan
 		number of freq. per decade, number of freq. in a lin. scan inside the log scan*/
 ap::template_1d_array< amp::campf<Precision> > eps1,eps1m,mu1,mu1m; /* eps1 and mu1 for upper (without m at
 					the end) and lower layers (with m at the end) */
 ap::template_1d_array< amp::ampf<Precision> > b,bm,thick,thickm; // position of upper and lower boundaries; thickness of the layers 
 ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,chi,fmu,rhom,taum,epsbm,chim,fmum; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<Precision> omega,beta,k,gamma,kovergamma,u,dummy3; // parameters
 amp::campf<Precision> jimagMP; // imaginary constant in multiprecision
 amp::ampf<Precision> mu0,eps0,Z0; // vacuum permeability and permittivity in SI units, and free space impedance
 gsl_integration_workspace *w;
 size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
 double tolintrel=1.e-15; // relative error permitted (w.r.t. to value at the previous frequency) for integration
 double x,y,L,fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15],*freq,dif,dummy1;
 std::complex<double> Zxdip,Zydip,Zxquad,Zyquad,Zlong,Zycst,cst;
 std::complex<double> alpha00,alpha01,alpha02,alpha11,jimag; // alphamn constants and imaginary constant
 std::complex<double> tolintabs00,tolintabs01,tolintabs02,tolintabs11;
 size_t found;
 time_t start,end; // times
					
 // start time
 time(&start);
 
 // memory of etas and chis allocation
 kxmem.setbounds(0,MAXMEM);
 eta1mem.setbounds(0,MAXMEM);eta2mem.setbounds(0,MAXMEM);
 chi1mem.setbounds(0,MAXMEM);chi2mem.setbounds(0,MAXMEM);
 mem=0;

 // constants
 mu0=4e-7*amp::pi<Precision>();
 eps0=1/(mu0*amp::sqr(C));
 Z0=mu0*C;
 jimagMP.x=0;jimagMP.y=1;
 jimag=std::complex<double>(0.,1.);
 
 // default values of the parameters (in case)
 flag_topbotsym=1;
 typescan=0;
 fminlog=2;fmaxlog=13;nflog=10;n_added=0;
 fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;nf=0;
 N=2;M=2;L=1.;gamma="479.6";
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   /* next line is when the input file comes from windows or else and has some ^M characters in the
   end of each line */
   //data[n_input+1]=data[n_input+1].substr(0,data[n_input+1].length()-1);
   //std::cout << data[n_input+1] << '\n';
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
   read_input(data[i],"start frequency exponent (10^)",dummy0,fminlog,dummy2,dummy3,1);
   read_input(data[i],"stop frequency exponent (10^)",dummy0,fmaxlog,dummy2,dummy3,1);
   read_input(data[i],"linear (1) or logarithmic (0) or both (2) frequency scan",typescan,dummy1,dummy2,dummy3,0);
   read_input(data[i],"sampling frequency exponent (10^) in Hz (for linear)",dummy0,fsamplin,dummy2,dummy3,1);
   read_input(data[i],"Number of points per decade (for log)",nflog,dummy1,dummy2,dummy3,0);
   read_input(data[i],"when both, fmin of the refinement",dummy0,fminlin,dummy2,dummy3,1);
   read_input(data[i],"when both, fmax of the refinement",dummy0,fmaxlin,dummy2,dummy3,1);
   read_input(data[i],"when both, number of points in the refinement",nflin,dummy1,dummy2,dummy3,0);
   read_input(data[i],"Comments for the output files names",dummy0,dummy1,commentoutput,dummy3,2);
   if (data[i].find("added frequencies") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos){
       n_added=1;fadded[n_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
       while (fadded[n_added] != 0){
         n_added++;fadded[n_added]=strtod(endline,&endline);
       }
       n_added--;
     }
   }
 }
 //printf("%s %s %d %d \n",machine.c_str(),topbot.c_str(),N,flag_topbotsym);
 //printf("%13.8e %13.8e %d\n",fadded[1],fadded[2],n_added);
 //printf("%13.8e %13.8e %d %d\n",fminlog,fmaxlog,nflog,typescan);

 // flag for top bottom symmetry (1 if there is such a symmetry)
 flag_topbotsym= ((strcmp(topbot.c_str(),"yes")==0 || strcmp(topbot.c_str(),"y")==0) || strcmp(topbot.c_str(),"1")==0);

 
 eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);b.setbounds(1,N+1);thick.setbounds(1,N);
 rho.setbounds(1,N+1);tau.setbounds(1,N+1);epsb.setbounds(1,N+1);chi.setbounds(1,N+1);fmu.setbounds(1,N+1);
 if (flag_topbotsym) M=N;
 eps1m.setbounds(1,M+1);mu1m.setbounds(1,M+1);bm.setbounds(1,M+1);thickm.setbounds(1,M);
 rhom.setbounds(1,M+1);taum.setbounds(1,M+1);epsbm.setbounds(1,M+1);chim.setbounds(1,M+1);fmum.setbounds(1,M+1);

 // default values of the layers properties (in case)
 rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;thick(1)="Infinity";

 // find inner half gap(s) of the chamber
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner half gap in mm",dummy0,dummy1,dummy2,b(1),3);
   if (flag_topbotsym) bm(1)=b(1);
   else {
     read_input(data[i],"Layer -1 inner half gap in mm",dummy0,dummy1,dummy2,bm(1),3);
   }
 }
 bm(1)=-bm(1);
 //printf("%s %s \n",b(1).toDec().c_str(),bm(1).toDec().c_str());
 
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

 // units conversion to SI (fminlin and fmaxlin were in THz, tau was in ps, b in mm and fmu in MHz)
 fminlin*=1.e12;fmaxlin*=1.e12;
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
  
 // first layer (inside the chamber) is always vacuum
 eps1(1)=1;mu1(1)=1;
 eps1m(1)=1;mu1m(1)=1;
 // relativistic velocity factor beta
 beta=amp::sqrt(1-1/amp::sqr(gamma));
 //printf("%s %s\n",eps1(1).x.toDec().c_str(),eps1(1).y.toDec().c_str());
 
 // construct the frequency scan
 // first estimation of the number of frequencies (for memory allocation)
 switch(typescan) {
   case 0:
     nf=(int)ceil((fmaxlog-fminlog)*(double)nflog)+1+n_added;
     break;
   case 1:
     nf=(int)ceil((pow(10.,fmaxlog)-pow(10.,fminlog))/pow(10.,fsamplin))+1+n_added;
     break;
   case 2:
     nf=(int)ceil((fmaxlog-fminlog)*(double)nflog)+1+nflin+1+n_added;
     break;
   }
 freq=new double[nf];
 
 // constructs unsorted version of the array freq
 nf=0;freq[0]=-1.;
 if (typescan==1) {
   do {
     freq[nf]=pow(10.,fminlog)+(double)nf * pow(10.,fsamplin);
     nf++;
     } while(freq[nf-1]<pow(10.,fmaxlog));
   }
 else {
   do {
     freq[nf]=pow(10.,fminlog+(double)nf/(double)nflog);
     nf++;
     } while(freq[nf-1]<pow(10.,fmaxlog));
   if (typescan==2) {
     for (unsigned int i=0; i<=nflin; i++) {
       freq[nf]=fminlin+(double)(i) * (fmaxlin-fminlin)/(double)nflin;
       nf++;
       }
     }
   }
 for (unsigned int i=1;i<=n_added; i++) {
   freq[nf]=fadded[i];
   nf++;
   }
 // nf is now the real number of frequencies
 
 // sort the frequencies
 gsl_sort(freq, 1, nf);
 
 // remove duplicate frequencies
 nfdup=0;
 for (int i=nf-2;i>=0; i--) {
   if ( (freq[i+1]-freq[i]) ==0 ) {
       freq[i]=DBL_MAX;
       nfdup++;
     }
   }
 gsl_sort(freq, 1, nf);
 nf=nf-nfdup;
 //printf("%d\n",nf);
 //for (unsigned int i=0;i<nf;i++) printf("%d %13.8e\n",i,freq[i]);
 
 
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
 sprintf(Input,"InputData%s.dat",output);
 
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
 
 
 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 gsl_set_error_handler_off();
 
 alpha00=std::complex<double>(0.,0.);
 alpha01=std::complex<double>(0.,0.);
 alpha02=std::complex<double>(0.,0.);
 alpha11=std::complex<double>(0.,0.);
 tolintabs00=std::complex<double>(0.,0.);
 tolintabs01=std::complex<double>(0.,0.);
 tolintabs02=std::complex<double>(0.,0.);
 tolintabs11=std::complex<double>(0.,0.);
 
 // impedance computation at each frequency: beginning of the loop
 for (unsigned int i=0; i<nf; i++) {
    
   mem=0;  // initialize memory at eahc frequency
   omega=amp::twopi<Precision>()*freq[i];
   k=omega/(beta*C);
   kovergamma=k/gamma;
   
   // computes the layer properties for the angular freq. omega
   for (unsigned int p=2;p<=N+1; p++) {
     if (rho(p).isFiniteNumber()) {
       eps1(p)=epsb(p)+1/(jimagMP*eps0*rho(p)*omega*(1+jimagMP*omega*tau(p)));
     } else {
       eps1(p)=epsb(p);
     }
     mu1(p)=1+chi(p)/(1+jimagMP*omega/(amp::twopi<Precision>()*fmu(p)));
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
 
   // computes alpha00
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,0,tolintabs00.imag(),limit,w);
   alpha00=std::complex<double>(x,y);
   //printf("alpha00: %13.8e %13.8e\n",alpha00.real(),alpha00.imag());
 
   // computes alpha01
   if (flag_topbotsym==0) {
     // computes alpha01
     x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.real(),limit,w);
     y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,1,tolintabs01.imag(),limit,w);
     alpha01=std::complex<double>(x,y);
     }
   else alpha01=std::complex<double>(0.,0.);
   //printf("alpha01: %13.8e %13.8e\n",alpha01.real(),alpha01.imag());
 
   // computes alpha02
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega, k,kovergamma,0,2,tolintabs02.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,0,2,tolintabs02.imag(),limit,w);
   alpha02=std::complex<double>(x,y);
   //printf("alpha02: %13.8e %13.8e\n",alpha02.real(),alpha02.imag());

   // computes alpha11
   x=integrate(1,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.real(),limit,w);
   y=integrate(0,M,N,b,bm,beta,eps1,eps1m,mu1,mu1m,omega,k,kovergamma,1,1,tolintabs11.imag(),limit,w);
   alpha11=std::complex<double>(x,y);
   //printf("alpha11: %13.8e %13.8e\n",alpha11.real(),alpha11.imag());
   
   // computes and writes the impedances
   cst=jimag*L*double(amp::ampf<Precision>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<Precision>())).toDouble());
   Zlong=cst*alpha00;
   Zycst=cst*alpha01/double(amp::ampf<Precision>(gamma).toDouble());
   cst=cst*double(amp::ampf<Precision>(k/(2*amp::sqr(gamma))).toDouble());
   Zxdip=cst*(alpha02-alpha00);
   Zydip=2.*cst*alpha11;
   Zxquad=-Zxdip;
   Zyquad=cst*(alpha02+alpha00);
   fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],Zlong.real(),Zlong.imag());
   fprintf(filZycst,"%13.8e %13.8e %13.8e\n",freq[i],Zycst.real(),Zycst.imag());
   fprintf(filZxdip,"%13.8e %13.8e %13.8e\n",freq[i],Zxdip.real(),Zxdip.imag());
   fprintf(filZydip,"%13.8e %13.8e %13.8e\n",freq[i],Zydip.real(),Zydip.imag());
   fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],Zxquad.real(),Zxquad.imag());
   fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],Zyquad.real(),Zyquad.imag());
 
  
   tolintabs00=tolintrel*std::complex<double>(std::abs(alpha00.real()),std::abs(alpha00.imag()));  
   tolintabs01=tolintrel*std::complex<double>(std::abs(alpha01.real()),std::abs(alpha01.imag()));  
   tolintabs02=tolintrel*std::complex<double>(std::abs(alpha02.real()),std::abs(alpha02.imag()));  
   tolintabs11=tolintrel*std::complex<double>(std::abs(alpha11.real()),std::abs(alpha11.imag()));
   //std::cout << "tolintabs02: " << tolintabs02.real() << " " << tolintabs02.imag() << "\n";

   }
   // end of loop on frequencies
 
 /*for (int i=1;i<=10000; i++){
   u=i*0.01;
   alpha00=integrand(0, 0, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, u, omega, beta, k, kovergamma);
   //alpha00=integrand(0, 0, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, (1-u)/u, omega, beta, k,kovergamma)/amp::sqr(u);
   alpha02=integrand(0, 2, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, u, omega, beta, k, kovergamma);
   //alpha02=integrand(0, 2, N+1, M+1, eps1, mu1, b, eps1m, mu1m, bm, (1-u)/u, omega, beta, k,kovergamma)/amp::sqr(u);
   printf("%13.8e %13.8e %13.8e\n",double(amp::ampf<Precision>(u).toDouble()),alpha02.real(),alpha02.imag());
 }*/
 
 
 gsl_integration_workspace_free(w);
  
 fclose(filZxdip);
 fclose(filZydip);
 fclose(filZxquad);
 fclose(filZyquad);
 fclose(filZlong);
 fclose(filZycst);
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
