/* 
 *
 * DELPHI.cc
 * Discrete Expansion over Laguerre Polynomials and Headtail modes for Instabilities computation
 *
 * by Nicolas Mounet (nicolas.mounet@cern.ch, nicolas.mounet@m4x.org, followed by imp-inst@cern.ch)
 *
 *
 * For a given impedance and damper gain, solve eigenvalue problem  in transverse, with mode coupling,
 including a bunch-by-bunch damper (see also work by A Burov), in multibunch, for any kind of longitudinal distribution
 given through its expansion on Laguerre polynomial. Uses associated Laguerre polynomial decomposition of the 
 radial function (R_l(r) in Chao's book, chap. 6) as done also by G. Besnier (1982) and Y. Chin (MOSES code, 1985-1988).
 
 
 It gives in output two files (name specified in input):
 - [output_filename].dat: one-line ascii file with mmax and nmax, maximum headtail mode and maximum degree of the 
 associated Laguerre polynomial in the diagonalized matrices finally used
 - [output_filename]_val.dat: eigenvalues Deltaomega (complex betatron 
 angular frequency shift w.r.t. the unperturbed tune), array of size M*(2*mmax+1)*(2*nmax+1)
 (nth radial mode of the mth headtail mode of the Mth multibunch mode). It is an ascii file
 with the array written in "blocks".

(UNUSED: - [output_filename]_coupl.dat: coupling matrix for the impedance (to be multiplied by intensity). This is 
 an array of size M*((2*mmax+1)*(2*nmax+1))^2. It is an ascii file with the matrix written in "blocks". The coupled-bunch
 is written at the beginning of each block of size (2*mmax+1)*(2*nmax+1))^2.
 - [output_filename]_damper.dat: coupling matrix for the damper (to be multiplied by damper gain). This is 
 an array of size ((2*mmax+1)*(2*nmax+1))^2. It is an ascii file with the matrix written in "blocks". )


 
Typical input file:
 

Impedance filename	Z.dat
Longitudinal distribution filename	no
Damper impedance filename	no
Output filename	out
Total bunch length [seconds]	2.5e-9
Number of particles per bunch	5e10
Machine circumference [m]	6911
Relativistic gamma	27.7
Transverse tune	26.129
Synchrotron tune	7.25e-3
Type of particle	proton
Chromaticity (DeltaQ*p/Q*Deltap)	0.1
Momentum compaction factor	1.92e-3
Number of bunches	924
Maximum number of eigenvalues	5
Minimum coupled-bunch mode number to consider	800
Maximum coupled-bunch mode number to consider	10
Coupled-bunch modes added	400 410 420
Maximum damper gain (inverse of number of damping turns)	0.02
Damper phase w.r.t. pi/2 (rad)	0
Parameter a	8
Parameter b	1
Use trapz method	yes
Convergence criterion	5.e-2
Maximum number of radial modes	2
Maximum azimuthal mode number	2

(UNUSED: Use precomputed coupling matrix from impedance	yes
Use precomputed coupling matrix from damper	yes)
 
 
The order of the lines can be whatever, but the exact sentences and the TAB before the parameter indicated, are necessary. 


Some explanatations for the input file:

 - Impedance filename: name of the file containing the dipolar impedance (3 columns, without header: frequency - NOT ANGULAR, real part of the impedance and imaginary part). It should be sorted in ascending order of frequencies. Frequencies should be positive.
 - Initial longitudinal distribution filename: file containing the longitudinal distribution, in terms of coefficients over an expansion 
 on Laguerre polynomials g0(tau)= exp(-b*tau^2/taub^2) sum_(k=0)^N g_k L_k(a*tau^2/taub^2) (for Gaussian, use "no" -> only one term, a can be whatever,
  b=8 - input file value overridden if different, g0= 8/(pi*taub^2) with taub the full bunch length (4 RMS) (Laclare's conventions)
 - Damper impedance filename: either "no" (no filename, bunch-by-bunch damper assumed), or name of the file containing the damper impedance (i.e. frequency dependent gain) (3 columns, without headers: frequency - NOT ANGULAR, real part and imag. part of the impedance). It should be sorted in ascending order of frequencies. Frequencies should span both positive and negative domain.
 - Output filename: eigenvalues are put in [this_filename]_val.dat, eigenvectors in [this_filename]_vec.dat.
 - Type of particle: proton or electron.
 - Maximum number of eigenvalues: kmax: total number of eigenvalues that are kept and accurate within "Convergence criterion" below(choose the ones with lowest imaginary part)
 - Minimum coupled-bunch mode number to consider: nxmin: we compute modes of number between nxmin and M-1 (total number of bunches minus one)
 - Maximum coupled-bunch mode number to consider: nxmax: we also compute modes of number between 0 and nxmax
 - Coupled-bunch modes added: we also compute those modes
 - Maximum damper gain in 1/nturns units: dmax: gain of the damper, expressed as a damping rate=1/(nb damping turns), 
 for mode 0 at 0 chromaticity. Note: this is used to normalize the damper matrix also when using a
 damper impedance function rather than a bunch-by-bunch ideal damper.
 - a and b: parameters for the initial longitudinal distribution (see above) and the Laguerre expansion of the modes. Can play with them 
 to make the code run faster (with a=b formulas are simpler but convergence might be slower).
 - Use trapz method: "no" to perform the impedance sum with "brute force", "yes" to try to speed it up by replacing part of it with an integral. Accuracy in the same, "no" 
 will be better with large number of bunches, "yes" is better in single-bunch.
 - Convergence criterion: crit: maximum error on the imaginary part of the kmax most unstable modes, between the current and the previous matrices diagonalized.
 - (optional) Maximum number of radial modes: nbrad: if used with max_azim, do not take into account convergence 
criterion and compute nbrad radial modes only.
 - (optional) Maximum azimuthal mode number: max_azim: if used with nbrad, do not take into account convergence 
criterion and compute azimuthal modes from -max_azim to +max_azim.

The rest is self-explanatory.

(UNUSED: If using precomputed matrices, they must be located in the same files as the corresponding output files (see above). )


 
*/

#include <iostream>
#include <fstream>
#include <stdio.h>
#include <string.h>
#include <time.h>
#include <vector>

#include <complex>
#include <cmath>
#include <math.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_errno.h>

#include "DELPHI.h"

#define  MAXLINES	30  // Maximum number of lines in the input file
#define  MAXLAG		20  // Maximum nb of Laguerre polynomials in the radial expansion (NOT USED)
#define  MAXHEAD	10  // Maximum headtail mode number (NOT USED)
#define  MAXSIZE	1000 // Maximum size of the matrix to diagonalize

using std::complex;

const double pi = 3.141592653589793238462643383279502884197; // pi
const complex<double> jimag = complex<double>(0.,1.); // imaginary unit

struct params { int l; int lprime; int n; int nprime; double omegaksi; double a; double b; double taub;
		long ng; double *g; long nZ; double *freq; double *Z;};

struct params_modif { int l; int lprime; int n; int nprime; double omegaksi; 
		double a; double b; double taub; long ng; double *g; long nZ; double *freq; double *Z;
		double omega1; double sign;};

extern "C" void zgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, complex<double> *W,
	complex<double> *VL, int *LDVL, complex<double> *VR, int *LDVR, complex<double> *WORK,
	int *LWORK, double *RWORK, int *INFO);


/******************************************************************************
 *** function locate: search a table ordered in ascending order		    ***
 ***		      (inspired from Numerical Recipes) 		    ***
 *** 		      Gives the position lprov (integer)		    ***
 ***                  in table, such that table[lprov-1]<=z<table[lprov]	    ***
 ******************************************************************************/

long locate (double *table, double z, long n)

{ long il, iu, im, lprov;

 il=0;
 iu=n+1;
 while (iu-il >1) {
   im=(iu+il)/2; // midpoint
   if (z >= table[im])
     il=im;
   else
     iu=im;
   }
 if (z==table[0]) lprov=1;
 else if (z==table[n]) lprov=n; 
 else lprov=iu;

 return lprov;

}


/*******************************************************************************
 *** function impedance: computes the impedance at omega		     ***
 ***			input: Z is the table with impedances (size=2*nZ:    ***
 ***			Z[2*i] is real part and Z[2*i+1] the imaginary part),***
 ***			freq is the table with the corresponding frequencies,***
 ***			omega is the angular frequency at which to compute   ***
 ***			the impedance (can be negative) and  nZ is the 	     ***
 ***			number of frequency values in the tables Z and freq. ***
 ***			It uses a linear interpolation of Z.	     ***
 *******************************************************************************/

complex<double> impedance (double *Z, double *freq, double omega, long nZ)

{
 complex<double> Zcomp;
 double reZ,imZ,f,factor;
 long lprov;
 
 f = omega/(2.0*pi);
 if (std::abs(f)>freq[nZ-1]) {
 	//printf ("\n\n    Need higher frequencies in the impedance file (replacing by 0) ! %lf\n",f);
      	return 0.; }
 lprov = locate(freq,std::abs(f),nZ);
 factor = (std::abs(f)-freq[lprov-1])/(freq[lprov]-freq[lprov-1]);
 reZ =  Z[2*(lprov-1)]+(Z[2*lprov]-Z[2*(lprov-1)])*factor;
 imZ =  Z[2*(lprov-1)+1]+(Z[2*lprov+1]-Z[2*(lprov-1)+1])*factor;
 if (f>0) {
 	Zcomp=reZ+jimag*imZ;
 }
 else {
 	if (f<0) {
		Zcomp=-reZ+jimag*imZ;
	}
	else {
		Zcomp=jimag*Z[1];
	}
 }
 
 return Zcomp;

}
 
 
/*******************************************************************************
 *** function damper_gain: computes the "damper impedance" at the frequency f***
 ***			input: d is the "damper impedance" table (size=2*nd: ***
 ***			d[2*i] is real part and d[2*i+1] the imaginary part),***
 ***			freq is the table with the corresponding frequencies,***
 ***			omega is the angular frequency at which to compute   ***
 ***			the impedance (can be negative) and nd is the 	     ***
 ***			number of frequency values in the tables d and freq. ***
 ***			It uses a linear interpolation of d.	     	     ***
 *******************************************************************************/

complex<double> damper_gain (double *d, double *freq, double omega, long nd)

{
 complex<double> dfinal;
 double red,imd,f,factor;
 long lprov;
 
 f = omega/(2.0*pi);
 if ( (f>freq[nd-1])||(f<freq[0]) ) {
 	//printf ("\n\n    Need lower or higher frequencies in the damper file (replacing by 0) ! %lf\n",f);
      	return 0.; }
 lprov = locate(freq,f,nd);
 factor = (f-freq[lprov-1])/(freq[lprov]-freq[lprov-1]);
 red =  d[2*(lprov-1)]+(d[2*lprov]-d[2*(lprov-1)])*factor;
 imd =  d[2*(lprov-1)+1]+(d[2*lprov+1]-d[2*(lprov-1)+1])*factor;
 dfinal=red+jimag*imd;
 
 return dfinal;

}

/*****************************************************************************************
 *** function to compute the factorial of a positive integer			       ***
 *****************************************************************************************/
 
long factorial(int n)

{
 long fact = 1;
 
 while (n>1) {
   fact*=n;
   n--;
 }
 
 return fact;
 
}

/*****************************************************************************************
 *** function to compute factorial(n)/factorial(p) (n must be larger than p)	       ***
 *****************************************************************************************/
 
long factorialpart(int n,int p)

 {
   long fact=1;
   
   for (int i=n; i>=p+1; i--) fact*=i;
   
   return fact;
   
}

/*****************************************************************************************
 *** function to compute (-1)^n							       ***
 *****************************************************************************************/

int minus1pow(int n)

{
 int res;
 
 if (n%2==0) res=1;
 else res= -1;
 
 return res;
 
}

/*****************************************************************************************
 *** function Laguerre: computes the associated Laguerre polynomial L_n^k(x) 	       ***
 *****************************************************************************************/

double Laguerre(int n, int k, double x)
  
{
 /* associated Laguerre polynomial computation */
 
 double coef,res;
 int iini;
 
 if (k>=0) {

     coef=(double)factorial(n+k)/(double)(factorial(k)*factorial(n));
     res=coef;
     iini=1;
 }

 else if (k>=-n) {

     coef=1./(double)factorial(-k)*std::pow(-x,-k);
     res=coef;
     iini=-k+1;
 }
 
 else {

     coef=(double)((double)minus1pow(n)*factorial(-k-1))/(double)(factorial(n)*factorial(-n-k-1));
     res=coef;
     iini=1;
 }

 for(int i=iini; i<=n; i++) {

     coef=coef*(double)(n-i+1)/(double)(i*(k+i))*(-x);
     res=res+coef;
 }
 
 return res;

}


/*****************************************************************************************
 *** function Gln: computes the function G_ln(omega,a) 		       ***
 *****************************************************************************************/

double Gln(int l, int n, double omega, double a, double taub, long ng, double *g)
  
{
 /* computes (2a)^(|l|+1) * taub^(|l|+2) * int_0^inf dtau * tau^(|l|+1) * g0(tau) * exp((b-a)*tau^2) * L_n^|l|(a*tau^2) * J_l(omega*tau)
 where g0(tau) is the initial distribution, L_n^|l| an associated Laguerre polynomial, J_l Bessel function.
 It is independent on b (cancels with exponential in g0(tau) ).
 
 Uses an analytical formula from Erdeliy et al (Table of Integral Transforms, vol. 2, p 43 (8) ).

For the different cases in the Laguerre calculations, uses formula 5.2.1 (p.102) from G. Szego, Orthogonal Polynomials, 1939
*/
 
 double res=0.,lag1,lag2,x=omega*omega/(4.*a);
 long fact;
 int epsl,absl=abs(l);

 if (omega==0.) {
 
   if ( (l!=0)||(n>ng) ) return 0.;
   else return g[n]*taub*taub;
 
 } else {

   if ( (l>=0)||((-l)%2==0) ) epsl=1;
   else epsl=-1;


   for (int k=0; k<ng; k++) {

     if (k<n) {

       fact=factorialpart(n,k);
       // fact is factorial(n)/factorial(k)
       lag1=Laguerre(k,n-k,x)*std::pow(-x,n-k)/(double)fact;

     } else lag1=Laguerre(n,k-n,x);

     if (n+absl<k) {

       fact=factorialpart(k,n+absl);
       // fact is factorial(k)/factorial(n+|l|)
       lag2=Laguerre(n+absl,k-n-absl,x)*std::pow(-x,k-n-absl)/(double)fact;

     } else lag2=Laguerre(k,n+absl-k,x);

     res+=g[k]*(double)minus1pow(n+k)*lag1*lag2;

   } 

   //return res*(double)epsl*std::pow(omega,absl)*std::exp(-x);
   return res*(double)epsl*std::pow(omega*taub,absl)*std::exp(-x)*taub*taub;

 }
 
}

 
/******************************************************************************************
 *** function Iln: computes the function I_ln(omega,a,b) 				***
 ******************************************************************************************/

double Iln(int l, int n, double omega, double a, double b, double taub)
  
{
 /* computes taub^(-|l'|-2) * int_0^inf dtau * tau^(|l|+1) * exp(-b*tau^2) * L_n^|l|(a*tau^2) * J_l(omega*tau)
 where L_n^|l| is an associated Laguerre polynomial, J_l Bessel function.
 
 Uses analytical formulas from Erdeliy et al (Table of Integral Transforms, vol. 2, p 42-43 (2)&(5) ). */
 
 double res=0.,x=omega*omega/(4.*b);
 long fact;
 int epsl,absl=abs(l);

 if (omega==0.) {
 
   if (l!=0) return 0.;
   else return std::pow(1.-a/b,n)/(2.*b*taub*taub);
 
 } else {
 
   if ( (l>=0)||((-l)%2==0) ) epsl=1;
   else epsl=-1;

   if (a==b) {

     //res=std::pow(omega/2.,2*n+absl)/(2.*(double)factorial(n));
     res=std::pow(omega/(2.*b*taub),absl)*std::pow(x,n)/(2.*(double)factorial(n)*b*taub*taub);

   } else {

     //res=std::pow(omega/2.,absl)*std::pow(b-a,n)*Laguerre(n,absl,a*x/(a-b))/2.;
     res=std::pow(omega/(2.*b*taub),absl)*std::pow(1.-a/b,n)*Laguerre(n,absl,a*x/(a-b))/(2.*b*taub*taub);

   }


   //return res*(double)epsl*std::pow(b,-absl-n-1)*std::pow(taub,-absl-2)*std::exp(-x);
   return res*(double)epsl*std::exp(-x);
 
 }
 
}

 
/******************************************************************************************
 *** function integrand: computes the term inside the integral in omega 		***
 ***			(it is also the term in the impedance sum)	     	   	***
 ******************************************************************************************/

complex<double> integrand(double omegap, double omegaksi, int l, int lprime, int n, int nprime, double a,
	double b, double taub, double *g, long ng, double *Z, double *freq, long nZ)
  
{
	
 /* function computing the integrand */

 double I_lprimenprime,G_ln,omega;
 complex<double> Zomegap,inte;

 
 omega=omegap-omegaksi; // for argument of I_lprimenprime and G_ln functions

 // computes impedance, I_lprimenprime and G_ln functions
 Zomegap=impedance(Z,freq,omegap,nZ);
 I_lprimenprime=Iln(lprime,nprime,omega,a,b,taub);
 G_ln=Gln(l,n,omega,a,taub,ng,g);
 
 inte=Zomegap*G_ln*I_lprimenprime;
 
 //printf("Inside integrand: l=%d l'=%d n=%d n'=%d %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",
 //	l,lprime,n,nprime,omega,inte.real(),inte.imag(),G_ln,I_lprimenprime,Zomegap.real(),Zomegap.imag(),omegaksi,a,b);
 
 /*FILE *fil;
 if ((l==1)&&(lprime==1)) {
   fil=fopen("tmp11.dat","a");
   fprintf(fil,"%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",omegap,
 	  Zomegap.real(),Zomegap.imag(),G_ln,I_lprimenprime,inte.real(),inte.imag());
   fclose(fil);
 }
 if ((l==2)&&(lprime==2)) {
   fil=fopen("tmp22.dat","a");
   fprintf(fil,"%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",omegap,
 	  Zomegap.real(),Zomegap.imag(),G_ln,I_lprimenprime,inte.real(),inte.imag());
   fclose(fil);
 }*/

 return inte;
 
}


/******************************************************************************************
 *** function integrand_real: computes the term inside the integral in omega (real part)***
 ***			(replacing part of the sum on the impedance)	     	   	***
 ***			This is in the appropriate format for gsl integration	   	***
 ******************************************************************************************/

double integrand_real(double omegap, void *p)
  
{
	
 /* function computing the real part of the integrand, for gsl integration */

 struct params_modif *param=(struct params_modif *)p;
 int l,lprime,n,nprime;
 long nZ,ng;
 double a,b,taub,omegaksi,*Z,*freq,*g;
 complex<double> inte;

 l=(param->l); // first headtail mode number
 lprime=(param->lprime); // second headtail mode number
 n=(param->n); // first radial mode number
 nprime=(param->nprime); // second radial mode number
 a=(param->a);
 b=(param->b);
 taub=(param->taub);
 omegaksi=(param->omegaksi); // chromatic angular frequency

 ng=(param->ng); // number of elements in Laguerre polynomial decomposition of initial distribution
 g=(param->g); // table of coefficients for Laguerre polynomial decomposition of initial distribution

 nZ=(param->nZ); // number of elements in impedance tables
 freq=(param->freq); // frequency table for impedance
 Z=(param->Z); // impedance table
 
 inte=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
 //printf("Integrand(1): %13.8e %13.8e\n",inte.real(),inte.imag());

 return inte.real();
 
}


/******************************************************************************************
 *** function integrand_imag: computes the term inside the integral in omega (imag part)***
 ***			(replacing part of the sum on the impedance)	     	   	***
 ***			This is in the appropriate format for gsl integration	   	***
 ******************************************************************************************/

double integrand_imag(double omegap, void *p)
  
{
	
 /* function computing the imaginary part of the integrand, for gsl integration */

 struct params_modif *param=(struct params_modif *)p;
 int l,lprime,n,nprime;
 long nZ,ng;
 double a,b,taub,omegaksi,*Z,*freq,*g;
 complex<double> inte;

 l=(param->l); // first headtail mode number
 lprime=(param->lprime); // second headtail mode number
 n=(param->n); // first radial mode number
 nprime=(param->nprime); // second radial mode number
 a=(param->a);
 b=(param->b);
 taub=(param->taub);
 omegaksi=(param->omegaksi); // chromatic angular frequency

 ng=(param->ng); // number of elements in Laguerre polynomial decomposition of initial distribution
 g=(param->g); // table of coefficients for Laguerre polynomial decomposition of initial distribution

 nZ=(param->nZ); // number of elements in impedance tables
 freq=(param->freq); // frequency table for impedance
 Z=(param->Z); // impedance table
 
 inte=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
 //printf("Integrand(2): %13.8e %13.8e\n",inte.real(),inte.imag());
 
 return inte.imag();
 
}


/******************************************************************************************
 *** function integrand_real_modif: computes the term inside the integral, with the 	***
 ***			            change in variable omegap=(omega1+sign*(1-t)/t)	***
 ***			This is in the appropriate format for gsl integration 	   	***
 ******************************************************************************************/

double integrand_real_modif(double t, void *p)
  
{
	
 /* function computing the real part of the integrand, for gsl integration, with change of variable */

 struct params_modif *param=(struct params_modif *)p;
 int l,lprime,n,nprime;
 long nZ,ng;
 double a,b,taub,omegaksi,*Z,*freq,*g,omega1,sign,omegap;
 complex<double> inte;

 l=(param->l); // first headtail mode number
 lprime=(param->lprime); // second headtail mode number
 n=(param->n); // first radial mode number
 nprime=(param->nprime); // second radial mode number
 a=(param->a);
 b=(param->b);
 taub=(param->taub);
 omegaksi=(param->omegaksi); // chromatic angular frequency

 ng=(param->ng); // number of elements in Laguerre polynomial decomposition of initial distribution
 g=(param->g); // table of coefficients for Laguerre polynomial decomposition of initial distribution

 nZ=(param->nZ); // number of elements in impedance tables
 freq=(param->freq); // frequency table for impedance
 Z=(param->Z); // impedance table
 
 omega1=(param->omega1); // bound of integration (in omegap) the closest to zero
 sign=(param->sign); // 1 if integration in omegap is from omega1 to +inf, -1 if integration in omegap is from -inf to omega1
 
 omegap=(omega1+sign*(1.-t)/t);
 
 inte=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
 //printf("Integrand(1): %13.8e %13.8e\n",inte.real(),inte.imag());

 return inte.real()/(t*t);
 
}


/******************************************************************************************
 *** function integrand_imag_modif: computes the term inside the integral, with the 	***
 ***			            change in variable omegap=(omega1+sign*(1-t)/t)	***
 ***			This is in the appropriate format for gsl integration	   	***
 ******************************************************************************************/

double integrand_imag_modif(double t, void *p)
  
{
	
 /* function computing the imaginary part of the integrand, for gsl integration, with change of variable */

 struct params_modif *param=(struct params_modif *)p;
 int l,lprime,n,nprime;
 long nZ,ng;
 double a,b,taub,omegaksi,*Z,*freq,*g,omega1,sign,omegap;
 complex<double> inte;

 l=(param->l); // first headtail mode number
 lprime=(param->lprime); // second headtail mode number
 n=(param->n); // first radial mode number
 nprime=(param->nprime); // second radial mode number
 a=(param->a);
 b=(param->b);
 taub=(param->taub);
 omegaksi=(param->omegaksi); // chromatic angular frequency

 ng=(param->ng); // number of elements in Laguerre polynomial decomposition of initial distribution
 g=(param->g); // table of coefficients for Laguerre polynomial decomposition of initial distribution

 nZ=(param->nZ); // number of elements in impedance tables
 freq=(param->freq); // frequency table for impedance
 Z=(param->Z); // impedance table
 
 omega1=(param->omega1); // bound of integration (in omegap) the closest to zero
 sign=(param->sign); // 1 if integration in omegap is from omega1 to +inf, -1 if integration in omegap is from -inf to omega1
 
 omegap=(omega1+sign*(1.-t)/t);
 
 inte=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
 //printf("Integrand(2): %13.8e %13.8e\n",inte.real(),inte.imag());
 
 return inte.imag()/(t*t);
 
}


/**********************************************************************************
 *** function impedanceint: computes the impedance integral vs. omega, for a	***
 ***			    certain chromatic frequency	omegaksi, headtail mode	***
 ***			    numbers l and lprime, Laguerre polynomial orders n	***
 ***			    and nprime, parameters (for Lag. expansion) a and b.***
 ***			    Integral computed between omega1 and omega2. If 	***
 ***			    omega2=0, it is replaced by	+/-Infinity.		***
 ***			    Laguerre decomposition of initial distribution is   ***
 ***			    in table g, of length ng.				***
 ***			    Impedance tables are in Z and freq, of length nZ.	***
 ***			    Uses "brute-force" gsl integration			***
 **********************************************************************************/

complex<double> impedanceint(double omega1, double omega2, int l, int lprime, int n, int nprime, double a, double b, double taub,
	double omegaksi, double *g, long ng, double *Z, double *freq, long nZ)
  
{

 gsl_integration_workspace *w; // workspace for gsl adaptative integration
 gsl_function F;  // gsl function
 double tolintabs=1.e-30,tolint=1.e-4; // absolute and relative error permitted for gsl adaptative integration
 //struct params param; // input parameters for the integrand functions (gsl integration)
 struct params_modif param; // input parameters for the integrand functions (gsl integration)
 size_t limit=1000,npts; // limit (number of intervals) for gsl integration algorithm
 double integ,err,real_int,imag_int,pts[2];
 int flagreal,status;

 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 // extrema of integration in case of semi-infinite integral and change of variable
 pts[0]=0.;pts[1]=1.;
 
 // parameters for gsl integration
 param.l=l;
 param.lprime=lprime;
 param.n=n;
 param.nprime=nprime;
 param.a=a;
 param.b=b;
 param.taub=taub;
 param.omegaksi=omegaksi;
 param.ng=ng;
 param.g=g;
 param.nZ=nZ;
 param.freq=freq;
 param.Z=Z;
 param.omega1=omega1;
 if (omega1>=0) param.sign = 1.;
 else param.sign = -1.;
 
 F.params=&param;


 for (flagreal=0; flagreal<=1; flagreal++) {
 
   if (flagreal) F.function=&integrand_real; // real part
   else F.function=&integrand_imag; // imaginary part
   
   if (omega2==0.) {

     if (omega1>=0) {

       // gsl adaptative integration on [omega1,+Inf[
       //printf("Integral on [%13.8e, +inf[\n",omega1);
       status=gsl_integration_qagiu(&F, omega1, tolintabs, tolint, limit, w, &integ, &err);

       } else {

       // gsl adaptative integration on ]-Inf,omega1]
       //printf("Integral on ]-inf, %13.8e]\n",omega1);
       status=gsl_integration_qagil(&F, omega1, tolintabs, tolint, limit, w, &integ, &err);

       }
       
       if (status) {

	 if ( ( (status==GSL_EMAXITER)||(status==GSL_EDIVERGE) )||(status==GSL_EROUND) ) {
	   
	   //printf("Warning: status=%d, GSL_EMAXITER=%d, GSL_EDIVERGE=%d, GSL_EROUND=%d\n",status,GSL_EMAXITER,GSL_EDIVERGE,GSL_EROUND);
	   // try on finite interval with change of variable omegap=(omega1+sign(omega1)*(1-t)/t)
	   if (flagreal) F.function=&integrand_real_modif;
	   else F.function=&integrand_imag_modif;
           
	   status=gsl_integration_qagp(&F, pts, 2, tolintabs, tolint, limit, w, &integ, &err);
	   //if ( (status)&&(std::abs(err/integ)>tolint) ) printf("Warning: l=%d, lprime=%d, omega1= %13.8e, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",l,lprime,omega1,flagreal,integ,std::abs(err/integ));
	   
	   } else {
           //printf("Warning: strange error on semi-infinite GSL integration: l=%d, lprime=%d, omega1= %13.8e, flagreal= %d, result= %13.8e, rel. error= %13.8e\n",l,lprime,omega1,flagreal,integ,std::abs(err/integ));
	   }
	 }


     } else {

     // gsl adaptative integration on [omega1,omega2]
     gsl_integration_qag(&F, omega1, omega2, tolintabs, tolint, limit, 1, w, &integ, &err);
     // key=1 -> Gauss-Kronrod with 15 points (Seems optimal)

     }
     
     if (flagreal) real_int=integ;
     else imag_int=integ;

   }
   
   //printf("l=%d, l'=%d, n=%d, n'=%d, omega1= %13.8e, omega2=%13.8e, integral: %13.8e %13.8e\n",l,lprime,n,nprime,omega1,omega2,real_int,imag_int);
   return real_int+jimag*imag_int;

}


/**********************************************************************************
 *** function impedancetrapz: computes the impedance integral vs. omega, for a	***
 ***			    certain chromatic frequency	omegaksi, headtail mode	***
 ***			    numbers l and lprime, Laguerre polynomial orders n	***
 ***			    and nprime, parameters (for Lag. expansion) a and b.***
 ***			    Integral computed between omega1 and omega2. If 	***
 ***			    omega2=0, it is replaced by	+/-Infinity.		***
 ***			    Laguerre decomposition of initial distribution is   ***
 ***			    in table g, of length ng.				***
 ***			    Impedance tables are in Z and freq, of length nZ.	***
 ***			    Uses a simple trapezoid method			***
 **********************************************************************************/

complex<double> impedancetrapz(double omega1, double omega2, int l, int lprime, int n, int nprime, double a, double b, double taub,
	double omegaksi, double *g, long ng, double *Z, double *freq, long nZ)
  
{

 double tolint=1.e-5; // relative error permitted for integration on semi-infinite interval
 complex <double> integ,integold,Z_up,Z_down;
 double G_ln_down,I_lprimenprime_down,G_ln_up,I_lprimenprime_up,omegadown,omegaup,omegatmp,sgn,omegalim;
 long j,ind,npts=100; // npts: we check convergence after every npts intervals
   
 integ=complex<double>(0.,0.);
 omegalim=std::abs(omegaksi)+2.*std::sqrt(std::max(a,b)); // limit in omega below which we should never stop the integration
 
 if (omega1>=0) {
   // integration on [omega1,omega2 or (+Inf)[
   sgn=1.;
   if (omega2==0.) omega2=1.e50;
   //printf("Integral (trapz.) on [%13.8e, %13.8e[\n",omega1,omega2);
 } else {
   // integration on ]omega1 (or -Inf),omega2 (or omega1)]
   sgn=-1.; 
   if (omega2==0.) omega2=-1.e50;
   else {
     omegatmp=omega1;
     omega1=omega2;
     omega2=omegatmp;
   }
   //printf("Integral (trapz.) on ]%13.8e, %13.8e]\n",omega2,omega1);
 }

 // index such that freq[ind-1]<=|omega1|/2pi<freq[ind]
 ind=locate(freq,std::abs(omega1)/(2.*pi),nZ);

 j=ind;
 //omegadown=freq[ind-1]*2.*pi;
 omegadown=std::abs(omega1);Z_down=impedance(Z,freq,omega1,nZ);
 G_ln_down=Gln(l,n,sgn*omegadown-omegaksi,a,taub,ng,g);
 I_lprimenprime_down=Iln(lprime,nprime,sgn*omegadown-omegaksi,a,b,taub);

 do {
   if ((j-ind)%npts==0) integold=integ;
   omegaup=freq[j]*2.*pi;
   G_ln_up=Gln(l,n,sgn*omegaup-omegaksi,a,taub,ng,g);
   I_lprimenprime_up=Iln(lprime,nprime,sgn*omegaup-omegaksi,a,b,taub);
   Z_up=sgn*Z[2*j]+jimag*Z[2*j+1];
   // add trapezoid area
   integ+=( Z_up*G_ln_up*I_lprimenprime_up + 
       	    Z_down*G_ln_down*I_lprimenprime_down )*(omegaup-omegadown)/2.;
   // go to next interval
   omegadown=omegaup;G_ln_down=G_ln_up;I_lprimenprime_down=I_lprimenprime_up;Z_down=Z_up;
   j++;

 } while ( (omegaup<omegalim) || ( ( (omegaup<std::abs(omega2)) && (j<nZ) ) && ( (std::abs((integ-integold)/integ)>tolint)||((j-ind)%npts!=0) ) ) );


 /*printf("ind=%ld, l=%d, l'=%d, n=%d, n'=%d, omega1= %13.8e, omega2=%13.8e, trapz: %13.8e %13.8e, old: %13.8e %13.8e, j=%ld, omegadown=%13.8e, omegaup=%13.8e\n",
 	ind,l,lprime,n,nprime,omega1,omega2,integ.real(),integ.imag(),integold.real(),integold.imag(),j,freq[ind-1]*2.*pi,omegaup);
 printf("Gln_up= %13.8e, Ilprimenprime_up=%13.8e, Gln_down= %13.8e, Ilprimenprime_down=%13.8e, Z_up: %13.8e %13.8e, Z_down: %13.8e %13.8e\n",
 	G_ln_up,I_lprimenprime_up,G_ln_down,I_lprimenprime_down,sgn*Z[2*j],Z[2*j+1],sgn*Z[2*j-2],Z[2*j-1]);*/

 return integ;

}


/**********************************************************************************
 *** function impedancesum: computes the impedance sum, for a certain		***
 ***			    coupled-bunch mode number nx, chromatic frequency	***
 ***			    omegaksi, number of bunches M, angular revolution	***
 ***			    frequency omega0, fractional tune Q, headtail mode	***
 ***			    numbers l and lprime, Laguerre polynomial orders n	***
 ***			    and nprime, parameters (for Lag. expansion) a and b.***
 ***			    Laguerre decomposition of initial distribution is   ***
 ***			    in table g, of length ng.				***
 ***			    abseps is the minimum absolute error tolerated.	***
 ***			    Impedance tables are in Z and freq, of length nZ.	***
 **********************************************************************************/

complex<double> impedancesum(int nx, int M, double omegaksi, double omega0, double Q, int l, int lprime,
	int n, int nprime, double a, double b, double taub, double *g, long ng, double *Z, double *freq, 
	long nZ, int flag_trapz, double abseps)
  
{

 complex<double> sum,sum_terms,oldsum=0.,sum_int=0.;
 double omegap,omegamp,omegap0,omegamp0,Momega0,Momega0over2,eps;
 long p,pmax,p0;
 clock_t tick1,tick2; // number of clock ticks
 

 //pmax=20;
 if (flag_trapz==1) eps=1.e-4;
 else eps=1.e-8;
 Momega0=(double)M*omega0;
 Momega0over2=Momega0/2.;

 //first term (p=0)
 omegap=((double)nx+Q)*omega0;
 omegamp=omegap;//-Momega0;
 //for (int i=0;i<=nZ-1;i++) printf("freq=%13.8e\n",freq[i]);
 //for (int i=0;i<=nZ-1;i++) printf("realZ=%13.8e, imagZ=%13.8e\n",Z[2*i],Z[2*i+1]);
 sum_terms=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
 //sum_terms+=integrand(omegamp, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ); // two terms already computed (one at positive freq, one at negative - to avoid initial imbalance)

 // initialization of correcting term sum_i (with an integral instead of discrete sum)
 //tick1 = clock();
 //sum_int =impedanceint(omegap+Momega0over2, 0.0, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
 //sum_int+=impedanceint(omegap-Momega0over2, 0.0, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
 if (flag_trapz==1) {
   sum_int =impedancetrapz(omegap+Momega0, 0.0, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   sum_int+=impedancetrapz(omegamp-Momega0, 0.0, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
 }
 //tick2 = clock();
 //printf("l= %d, lprime= %d: Elapsed time for calculation of impedance infinite integral: %ld ticks\n",l,lprime,tick2-tick1);

 // first total sum estimate
 sum = sum_int + sum_terms;
 
 p=0;
 //printf("sum_terms: %13.8e %13.8e, sum_int: %13.8e %13.8e\n",sum_terms.real(),sum_terms.imag(),sum_int.real(),sum_int.imag());

 tick1 = clock();
 //while ( ( (std::abs(sum.real()-oldsum.real()))>eps*std::abs(sum.real()) ) 
 //	|| ( (std::abs(sum.imag()-oldsum.imag()))>eps*std::abs(sum.imag()) ) ) {
 while ( ( (std::abs(sum.real()-oldsum.real()))>std::max(abseps,eps*std::abs(sum.real())) ) 
 	|| ( (std::abs(sum.imag()-oldsum.imag()))>std::max(abseps,eps*std::abs(sum.imag())) ) ) {

   oldsum=sum;//printf("p=%ld, p0=%ld, p0+pmax=%ld, omegap0=%13.8e, omegap=%13.8e, omegamp0=%13.8e, omegamp=%13.8e\n",p,p0,p0+pmax,omegap0,omegap,omegamp0,omegamp);
   p0=p;
   pmax=20+(long)floor(omegap/(10.*Momega0+std::abs(omegaksi)));//printf("pmax=%ld\n",pmax);
   omegap0=omegap;
   omegamp0=omegamp;

   //tick1 = clock();
   for(p=p0+1; p<=p0+pmax; p++) {

     omegap  += Momega0;
     omegamp -= Momega0;

     // add terms in the sum for p and -p
     sum_terms+=integrand(omegap, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);
     sum_terms+=integrand(omegamp, omegaksi, l, lprime, n, nprime, a, b, taub, g, ng, Z, freq, nZ);

     }
     
   //printf("sum_terms: %13.8e %13.8e, sum_int: %13.8e %13.8e\n",sum_terms.real(),sum_terms.imag(),sum_int.real(),sum_int.imag());
   
   p--; // because loop always finishes with p=p0+pmax+1
   //tick2 = clock();
   //printf("l= %d, lprime= %d: Elapsed time for calculation of impedance sum: %ld ticks\n",l,lprime,tick2-tick1);
   //printf("omegap0=%13.8e, omegap=%13.8e, omegamp0=%13.8e, omegamp=%13.8e, omegap+Momega0over2=%13.8e, omegamp-Momega0over2=%13.8e\n",
   //	omegap0,omegap,omegamp0,omegamp,omegap+Momega0over2,omegamp-Momega0over2);
   
   //sum_int =impedanceint(omegap+Momega0over2, 0., l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   //sum_int+=impedanceint(omegamp-Momega0over2, 0., l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   if (flag_trapz==1) {
     sum_int =impedancetrapz(omegap+Momega0, 0., l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
     sum_int+=impedancetrapz(omegamp-Momega0, 0., l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   }
   sum = sum_int + sum_terms;
   
   //printf("sum: %13.8e %13.8e, sum_int: %13.8e %13.8e, oldsum: %13.8e %13.8e\n",
   //	sum.real(),sum.imag(),sum_int.real(),sum_int.imag(),oldsum.real(),oldsum.imag());

   // subtract correction (rest of the sum considered as integral -> should suppress redundant terms)
   //tick1 = clock();
   //sum_int =impedanceint(omegap0+Momega0over2, omegap+Momega0over2, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   //sum_int+=impedanceint(omegamp-Momega0over2, omegamp0-Momega0over2, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   //sum_int =impedancetrapz(omegap0+Momega0over2, omegap+Momega0over2, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   //sum_int+=impedancetrapz(omegamp-Momega0over2, omegamp0-Momega0over2, l, lprime, n, nprime, a, b, taub, omegaksi, g, ng, Z, freq, nZ)/Momega0;
   //tick2 = clock();
   //printf("l= %d, lprime= %d: Elapsed time for calculation of impedance integral: %ld ticks\n",l,lprime,tick2-tick1);

   //sum -= sum_int;
   
   }
   
 if ( (omegap>1.e14)||(omegamp<-1.e14) ) {
   printf("l= %d, lprime= %d: very large maximum freq. %13.8e THz, %13.8e THz, sum (real, imag)= %13.8e %13.8e, oldsum (real, imag)= %13.8e %13.8e\n",
      l,lprime,omegap/(1.e12*2.*pi),omegamp/(1.e12*2.*pi),sum.real(),sum.imag(),oldsum.real(),oldsum.imag());
   }
 /*tick2 = clock();
 printf("l=%d, lprime=%d, n=%d, nprime=%d\n",l,lprime,n,nprime);
 printf("p=%ld, omegap final=%13.8e, omegamp final=%13.8e\n",p,omegap,omegamp);	 
 printf("Elapsed time for calculation of impedance sum: %ld ticks, sum (real, imag.)=(%13.8e, %13.8e)\n",
 	tick2-tick1,sum.real(),sum.imag());*/
 return sum;
	
}

/**********************************************************************************
 *** function dampersum: computes the "damper impedance" sum, for a certain	***
 ***			 coupled-bunch mode number nx, number of bunches M,	***
 ***			 rev. frequency omega0, chromatic frequency omegaksi,	***
 ***			 fractional tune Q, headtail mode numbers l and lprime, ***
 ***			 Laguerre polynomial orders n and nprime, parameters	***
 ***			 (for Lag. expansion) a and b.				***
 ***			 Laguerre decomposition of initial distribution is   	***
 ***			 in table g, of length ng.				***
 ***			 abseps is the minimum absolute error tolerated.	***
 ***			 damper impedance tables are in d and freq, of length nd***
 **********************************************************************************/

complex<double> dampersum(int nx, int M, double omegaksi, double omega0, double Q, int l, int lprime,
	int n, int nprime, double a, double b, double taub, double *g, long ng, 
	double *d, double *freq, long nd, double abseps)
  
{

 double eps=1e-8;
 double omegap,omegamp,omegap0,omegamp0,Momega0,I_lprimenprime,G_ln;
 complex<double> dimp,sum,oldsum=0.;
 long p,p0,pmax;
 clock_t tick1,tick2; // number of clock ticks
 
 //tick1 = clock();

 Momega0=(double)M*omega0;
 omegap=((double)nx+Q)*omega0;
 omegamp=omegap;
 I_lprimenprime=Iln(lprime,nprime,omegap-omegaksi,a,b,taub);
 G_ln=Gln(l,n,omegap-omegaksi,a,taub,ng,g);
 dimp=damper_gain(d, freq, omegap, nd);
 sum=dimp*I_lprimenprime*G_ln;
 
 p=0;
 while ( ( (std::abs(sum.real()-oldsum.real()))>std::max(abseps,eps*std::abs(sum.real())) ) 
 	|| ( (std::abs(sum.imag()-oldsum.imag()))>std::max(abseps,eps*std::abs(sum.imag())) ) ) {

   oldsum=sum;//printf("p=%ld, p0=%ld, p0+pmax=%ld, omegap0=%13.8e, omegap=%13.8e, omegamp0=%13.8e, omegamp=%13.8e\n",p,p0,p0+pmax,omegap0,omegap,omegamp0,omegamp);
   p0=p;
   pmax=20+(long)floor(omegap/(10.*Momega0+std::abs(omegaksi)));//printf("pmax=%ld\n",pmax);
   omegap0=omegap;
   omegamp0=omegamp;

   for(p=p0+1; p<=p0+pmax; p++) {

     omegap  += Momega0;
     omegamp -= Momega0;

     // add terms in the sum for p and -p
     I_lprimenprime=Iln(lprime,nprime,omegap-omegaksi,a,b,taub);
     G_ln=Gln(l,n,omegap-omegaksi,a,taub,ng,g);
     dimp=damper_gain (d, freq, omegap, nd);
     sum+=dimp*I_lprimenprime*G_ln;

     I_lprimenprime=Iln(lprime,nprime,omegamp-omegaksi,a,b,taub);
     G_ln=Gln(l,n,omegamp-omegaksi,a,taub,ng,g);
     dimp=damper_gain (d, freq, omegamp, nd);
     sum+=dimp*I_lprimenprime*G_ln;

     }
     
   p--; // because loop always finishes with p=p0+pmax+1
 
 }

 /*tick2 = clock();
 printf("Elapsed time for calculation of impedance sum: %ld ticks, sum (real, imag.)=(%13.8e, %13.8e)\n",
 	tick2-tick1,sum.real(),sum.imag());*/

 return sum;
	
}


/**************************************************
*** read_input: read an input file line		***
***************************************************/

void read_input(std::string line, std::string description, int& param0, double& param1, 
  	char *param2, int type){
	
 /* function reading the number or string at the end of an input line, identifying the description string.
 type=0 if output is an int, 1 if output is a double, 2 if output is a char array.
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

	 //param2=line.substr(found+1,std::string::npos);
	 strcpy(param2,line.substr(found+1,std::string::npos).c_str());
	 break;

     }

   }

 }

 return;

}

/**************************************************
*** read_nZ_impedance: find number of lines in 	***
*** impedance file				***
***************************************************/

long read_nZ_impedance(char *filename){
 
 FILE *filZ;
 double dummy1,dummy2,dummy3;
 int condition;
 long nZ;

 // determines the number of lines nZ
 filZ = fopen(filename,"r");
 if (filZ == NULL)
   { printf ("\n\n    File with impedance does not exist.\n") ;
   exit(2); }
 nZ = 0;
 do {condition = fscanf (filZ, "%lg %lg %lg\n",&dummy1,&dummy2,&dummy3);
   ++nZ;	 
 } while (condition != EOF);
 nZ--;
 fclose(filZ);
 
 return nZ;
 
}

 
/**************************************************
*** read_impedance: read impedance file		***
***************************************************/

void read_impedance(char *filename, long& nZ, double *freq, double *Z){
 
 FILE *filZ;
 int condition;
 long j=0;
  
 // read frequencies and corresponding impedances
 filZ = fopen(filename,"r");
 
 do {condition = fscanf (filZ, "%lg %lg %lg\n",&freq[j],&Z[2*j],&Z[2*j+1]);
 	if ( (j>0) && (freq[j] == freq[j-1]) ) {
	  // remove duplicate frequencies
	  j--;
	  nZ--;
	}
 	j++;
 } while (condition != EOF);
 fclose(filZ);
 
 //for (j=0;j<=nZ-1;j++) printf("%ld %ld %13.8e %13.8e %13.8e\n",nZ,j,freq[j],Z[2*j],Z[2*j+1]);
 
 return;

}


/**************************************************
*** read_nd_damper: find number of lines in 	***
*** damper file					***
***************************************************/

long read_nd_damper(char *filename){
 
 FILE *fild;
 char dummy;
 double dummy1,dummy2;
 int condition;
 long nd;

 // determines the number of lines nd
 fild = fopen(filename,"r");
 if (fild == NULL)
   { printf ("\n\n    File with damper gain data does not exist.\n") ;
   exit(2); }
 nd = 0;
 fscanf (fild,"%s\n",&dummy); // take out first line (headers)
 do {condition = fscanf (fild, "%lg %lg\n",&dummy1,&dummy2);
   ++nd;	 
 } while (condition != EOF);
 nd--;
 fclose(fild);
 
 return nd;
 
}

 
/**************************************************
*** read_damper: read damper file		***
***************************************************/

void read_damper(char *filename, long& nd, double *freqd, double *d){
 
 FILE *fild;
 char dummy;
 int condition;
 long j=0;
  
 // read frequencies and corresponding impedances
 fild = fopen(filename,"r");
 fscanf (fild,"%s\n",&dummy); // take out first line (headers)
 
 do {condition = fscanf (fild, "%lg %lg\n",&freqd[j],&d[j]);
 	if ( (j>0) && (freqd[j] == freqd[j-1]) ) {
	  // remove duplicate frequencies
	  j--;
	  nd--;
	}
 	j++;
 } while (condition != EOF);
 fclose(fild);
 
 for (j=0;j<=nd-1;j++) printf("%ld %ld %13.8e %13.8e\n",nd,j,freqd[j],d[j]);

 return;
 
}

/**************************************************
*** read_ng_longdist: find number of lines in 	***
*** longitudinal distribution file					***
***************************************************/

long read_ng_longdist(char *filename){
 
 FILE *filg;
 double dummy1;
 int condition;
 long ng;

 // determines the number of lines nd
 filg = fopen(filename,"r");
 if (filg == NULL)
   { printf ("\n\n    File with longitudinal distribution does not exist.\n") ;
   exit(2); }
 ng = 0;
 do {condition = fscanf (filg, "%lg\n",&dummy1);
   ++ng;	 
 } while (condition != EOF);
 ng-=2; // note: here ng is the number of elements-1 (we go from 0 to ng)
 fclose(filg);

 return ng;

}


/*********************************************************
*** read_longdist: read longitudinal distribution file ***
**********************************************************/

void read_longdist(char *filename, long& ng, double *g){
 
 FILE *filg;
 int condition;
 long j=0;
 
  // read frequencies and corresponding impedances
 filg = fopen(filename,"r");
 
 do {condition = fscanf (filg, "%lg\n",&g[j]);
 	j++;
 } while (condition != EOF);
 fclose(filg);
 
 for (j=0;j<=ng;j++) printf("%ld %ld %13.8e\n",ng,j,g[j]);

 return;
 
}

/******************************************************************************
 *** 			main program					    ***
 ******************************************************************************/

main ()

{

 const double e=1.60218e-19; // elementary charge in C
 const double clight=299792458; // speed of light in m/s
 double beta,gamma,f0,Ib,Nb,omega0,omegaksi,eta,omegas,tune,tunes,tunefrac,circum,chroma,alphap,taub,m0,dif,tunefracxi,dmax;
 /*double *freq,*omega,*Z,dummy1,dummy2,dummy3,rwork[2*MAXSIZE],integ,omegap,omegal,x,jml,jmlprime,errx,
 	Kxm[2*MAXSIZE*MAXSIZE],Kxmold[2*MAXSIZE*MAXSIZE],eigenvalabs[MAXSIZE];
 complex<double> Zcomp,eigenval[MAXSIZE],*eigenvalold,eigenvect[MAXSIZE*MAXSIZE],work[40*MAXSIZE],vl[1];*/
 double *freq,*Z,*d,*freqd,*g,*rwork,x,*Mat,*eigenvalimag,G_ln,I_lprimenprime,dphase;
 double dummy1,err,diff,maxdiff,fcutoff,a,b,aovertaub2,bovertaub2,fact,crit,abseps=1e-4;
 complex<double> *eigenval,*eigenvalold,*eigenvect,*coupl,*couplold,*damper,*damperold,*work,vl[1],Zsum,coefj,dfactor;
 char Zfilename[300],longfilename[300],dfilename[300],outputfilename[300],outputfilenameval[300],outputfilenamevect[300];//,outputfilenamecoupl[300],outputfilenamedamp[300];
 char particle[10],strcoupl[3],strdamp[3],right,left,*endline,*dummy2,strtrapz[3];
 std::string data[MAXLINES];
 FILE *fileigval, *fileigvect;//, *filcoupl, *fildamper;
 long j,nZ,ng,nd;
 int dummy0,n_input,i_part,flag_coupl,flag_damp,M,nx,nxi,nxmin,nxmax,nxadded[50],*nxiscan,n_added,nmodes,rem,flag_trapz;
 int kmax0,kmax,lmax,lmaxold,nmax,nmaxold,k,kerr,nbrad=0,max_azim=0;
 int p,i,l,lprime,n,nprime,size,sizel,sizeold,sizelold,info,lwork;
 int flagdamperfile,ldvl,flagcrit;
 size_t found,*ind;
 time_t start,end; // times
 clock_t tick1,tick2; // number of clock ticks


 time(&start);
 
 // read input file
 // first read everything to identify the strings in front of each parameters
 n_input=0;
 while (std::cin.eof()==0) {
   std::getline (std::cin,data[n_input+1]);
   /* next line is when the input file comes from windows or else and has some ^M characters in the
   end of each line */
   //data[n_input+1]=data[n_input+1].substr(0,data[n_input+1].length()-1);
   std::cout << data[n_input+1] << '\n';
   n_input++;
 }
 n_input--;
 
 // identify each line of the input file
 for (i=1; i<=n_input; i++) {
   read_input(data[i],"Impedance filename",dummy0,dummy1,Zfilename,2);
   read_input(data[i],"Longitudinal distribution filename",dummy0,dummy1,longfilename,2);
   read_input(data[i],"Damper impedance filename",dummy0,dummy1,dfilename,2);
   read_input(data[i],"Output filename",dummy0,dummy1,outputfilename,2);
   read_input(data[i],"Total bunch length",dummy0,taub,dummy2,1);
   read_input(data[i],"Number of particles per bunch",dummy0,Nb,dummy2,1);
   read_input(data[i],"Machine circumference",dummy0,circum,dummy2,1);
   read_input(data[i],"Relativistic gamma",dummy0,gamma,dummy2,1);
   read_input(data[i],"Transverse tune",dummy0,tune,dummy2,1);
   read_input(data[i],"Synchrotron tune",dummy0,tunes,dummy2,1);
   read_input(data[i],"Type of particle",dummy0,dummy1,particle,2);
   read_input(data[i],"Chromaticity",dummy0,chroma,dummy2,1);
   read_input(data[i],"Momentum compaction",dummy0,alphap,dummy2,1);
   read_input(data[i],"Number of bunches",M,dummy1,dummy2,0);
   read_input(data[i],"Maximum number of eigenvalues",kmax0,dummy1,dummy2,0);
   read_input(data[i],"Minimum coupled-bunch mode number to consider",nxmin,dummy1,dummy2,0);
   read_input(data[i],"Maximum coupled-bunch mode number to consider",nxmax,dummy1,dummy2,0);
   if (data[i].find("Coupled-bunch modes added") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos){
       n_added=0;nxadded[n_added]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
       while (nxadded[n_added] != 0){
         n_added++;nxadded[n_added]=strtod(endline,&endline);
         }
       }
     }
   read_input(data[i],"Maximum damper gain",dummy0,dmax,dummy2,1);
   //read_input(data[i],"Frequency cutoff of damper in MHz",dummy0,fcutoff,dummy2,1);
   //fcutoff*=1.e6;
   read_input(data[i],"Damper phase",dummy0,dphase,dummy2,1);
   read_input(data[i],"Parameter a",dummy0,a,dummy2,1);
   read_input(data[i],"Parameter b",dummy0,b,dummy2,1);
   read_input(data[i],"Use trapz method",dummy0,dummy1,strtrapz,2);
   read_input(data[i],"Convergence criterion",dummy0,crit,dummy2,1);
   read_input(data[i],"Maximum number of radial modes",nbrad,dummy1,dummy2,0);
   read_input(data[i],"Maximum azimuthal mode number",max_azim,dummy1,dummy2,0);

   /*if (data[i].find("Use precomputed coupling matrix from impedance") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(strcoupl,data[i].substr(found+1,std::string::npos).c_str());
     }
   if (data[i].find("Use precomputed coupling matrix from damper") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(strdamp,data[i].substr(found+1,std::string::npos).c_str());
     }*/

 }
 
 /*printf("%s %s %s %d %s\n",distribution,outputfilename,particle,M,dfilename);
 printf("%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",tune,tunes,chroma,taub,gamma,circum,Nb,alphap);
 for (j=0;j<n_added;j++) printf("%ld ",nxadded[j]);
 printf ("%ld\n",n_added);*/

 // deactivate abort on GSL errors
 gsl_set_error_handler_off();
  
 // parameters
 beta=sqrt(1.-1./(gamma*gamma)); // relativistic velocity factor
 f0=beta*clight/circum; // revolution frequency
 omega0=2.*pi*f0; // revolution angular frequency
 Ib=(double)M*e*Nb*f0; // beam current
 eta=alphap-1./(gamma*gamma); // transition parameter (Elias Metral's convention, opposite to Joel Le Duff's one)
 omegaksi=chroma*tune*omega0/eta; // chromatic angular frequency
 omegas=tunes*omega0; // synchrotron angular frequency
 tunefrac=tune-floor(tune); // fractional part of the tune
 
 // coupled-bunch mode numbers considered
 nxiscan=new int[M];
 i=0;
 for(nxi=nxmin;nxi<=M+nxmax;nxi++) {
   if (i<=M) nxiscan[i]=nxi%M;
   i++;
 }
 for(j=0;(int)j<n_added;j++) nxiscan[i+(int)j]=nxadded[j]%M;
 nmodes=i+n_added;
 //printf("%ld initial coupled-bunch modes\n",nmodes);
 gsl_sort_int(nxiscan,1,nmodes);
 // remove duplicates
 rem=0;
 for(i=1;i<nmodes;i++) {
   if (nxiscan[i]==nxiscan[i-1]) {
     nxiscan[i-1]=2*M;
     rem++;
   }
 }
 gsl_sort_int(nxiscan,1,nmodes); // all redundant modes are put at the end of the array
 nmodes -= rem; // now they are virtually suppressed
 printf("%d coupled-bunch modes\n",nmodes);
 for(i=0;i<nmodes;i++) printf("%d ",nxiscan[i]);
 printf("\n");   
 
 

 // default particle: proton
 i_part=0;
 m0=1.6726e-27;

 // beam particles' mass
 if (strcmp(particle,"proton")==0) {
	 i_part=0;
	 m0=1.6726e-27;
 }
 if (strcmp(particle,"electron")==0) {
	 i_part=1;
         m0=9.1094e-31;
 }

 // flag for trapz method
 if ( (strcmp(strtrapz,"yes")==0)||(strcmp(strtrapz,"y")==0) ) flag_trapz=1;
 else flag_trapz=0;
 //printf("%s flag_trapz=%d\n",strtrapz,flag_trapz);

 /*if ( (strcmp(strcoupl,"y")==0)||(strcmp(strcoupl,"1")==0) ) {
	 flag_coupl=1;
 } else {
 	 flag_coupl=0;
 }
 if ( (strcmp(strdamp,"y")==0)||(strcmp(strdamp,"1")==0) ) {
	 flag_damp=1;
 } else {
 	 flag_damp=0;
 }*/


 // read impedance file
 nZ=read_nZ_impedance(Zfilename);
 Z = new double[nZ*2];
 freq = new double[nZ];
 read_impedance(Zfilename, nZ, freq, Z);
 //for (j=0;j<=nZ-1;j++) printf("%ld %ld %13.8e %13.8e %13.8e\n",nZ,j,freq[j],Z[2*j],Z[2*j+1]);
  
 if (strcmp(dfilename,"no")==0) {

   flagdamperfile=0;

 } else {

   flagdamperfile=1;
   // read damper file
   nd=read_nZ_impedance(dfilename);
   d = new double[nd*2];
   freqd = new double[nd];
   read_impedance(dfilename, nd, freqd, d);
   
 }

 if (strcmp(longfilename,"no")==0) {

   // default longitudinal distribution: Gaussian
   b=8.;ng=0;
   g=new double[1];
   g[0]=8./(pi*taub*taub);

 } else {

   // read longitudinal distribution file
   ng=read_ng_longdist(longfilename);
   g = new double[ng+1];
   read_longdist(longfilename, ng, g);

 }

 aovertaub2=a/(taub*taub);
 bovertaub2=b/(taub*taub);
 // those above are the a and b to plug in the analytical formulas for the matrix coefficients

 // output filenames
 strcpy (outputfilenameval,outputfilename) ;
 strcpy (outputfilenamevect,outputfilename) ;
 /*strcpy (outputfilenamecoupl,outputfilename) ;
 strcpy (outputfilenamedamp,outputfilename) ;*/
 strcat (outputfilenameval,"_val.dat") ;
 strcat (outputfilenamevect,"_vec.dat") ;
 /*strcat (outputfilenamecoupl,"_coupl.dat") ;
 strcat (outputfilenamedamp,"_damper.dat") ;*/
 // open files
 fileigval=fopen(outputfilenameval,"w");
 fileigvect=fopen(outputfilenamevect,"w");
 /*filcoupl=fopen(outputfilenamecoupl,"w");
 fildamper=fopen(outputfilenamedamp,"w");*/


 /* note: compiling those arrays in static does not make any difference in term of computation time,
 but reallocating them each time inside the loops makes a big one. */
 eigenvalold=new complex<double>[kmax0];
 ind=new size_t[kmax0];
 eigenval=new complex<double>[MAXSIZE];
 eigenvalimag=new double[MAXSIZE];
 eigenvect=new complex<double>[MAXSIZE*MAXSIZE];
 Mat=new double[2*MAXSIZE*MAXSIZE];
 coupl=new complex<double>[MAXSIZE*MAXSIZE];
 couplold=new complex<double>[MAXSIZE*MAXSIZE];
 damper=new complex<double>[MAXSIZE*MAXSIZE];
 damperold=new complex<double>[MAXSIZE*MAXSIZE];
 work=new complex<double>[40*MAXSIZE];
 rwork=new double[2*MAXSIZE];

 ldvl=1;left='N';right='V';
 
 for(nxi=0; nxi<nmodes; nxi++) {

   nx=nxiscan[nxi];
   // nx: coupled-bunch mode number (each considered individually)
   // initialization of size of Laguerre polynomials expansion (-1: it begins at zero) and maximum headtail mode number
   nmax=0;lmax=max_azim;nmaxold=-1;lmaxold=-1;sizelold=-1;sizeold=-1;
   if (nbrad>1) nmax=nbrad-2;
   if (max_azim>0) lmax=max_azim-1;


   for (k=0; k<=kmax0-1; k++) eigenvalold[k]=complex<double>(1.e50,1.e50);

   // loop for convergence on matrix size
   //tick1 = clock();
   do {

     sizel=2*lmax+1;
     size=sizel*(nmax+1); // matrix size
     kmax=std::min(size,kmax0);//printf("kmax=%d\n",kmax);
     lwork=40*size; // note: not visible impact of lwork on computation time

     // construction of the coupling matrix
     /*tick1 = clock();
     time(&start);*/
     for(l=lmax; l>=-lmax; l--) {

       for(lprime=lmax; lprime>=-lmax; lprime--) {
       
         coefj=std::pow(jimag,lprime-l);

	 for(n=0; n<=nmax; n++) {
	 
	   fact=1./((double)factorialpart(n+abs(l),n)*std::pow(2.,abs(l))); // =factorial(n)/(factorial(n+|l|)*2^|l|)
	   
	   for(nprime=0; nprime<=nmax; nprime++) {

	     if ( (l>=0)&&(lprime>=0) ) {
	     
	       if ( ( (l>lmaxold)||(lprime>lmaxold) ) || ( (n>nmaxold)||(nprime>nmaxold) ) ) {
	       
	         // matrix coefficient was never computed for this l, lprime, n, nprime
		 Zsum=impedancesum(nx, M, omegaksi, omega0, tunefrac, l, lprime, n, nprime, aovertaub2, bovertaub2, taub, g, ng, Z, freq, nZ, flag_trapz, abseps);
		 		 
		 coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coefj*fact*Zsum;
		 //printf("l=%d, l'=%d, n=%d, n'=%d, coupl. matrix: %13.8e %13.8e\n",l,lprime,n,nprime,
		 //	coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].real(),coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].imag());
		 
	       } else {
	       
	         // matrix coefficient was already computed
		 coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=couplold[l+lmaxold+n*sizelold+(lprime+lmaxold+nprime*sizelold)*sizeold];
		 		 
	       }
	       
	     } else if ( (l<0)&&(lprime>=0) ) {
	     
	       coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coupl[-l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];
	     
	     } else if ( (l>=0)&&(lprime<0) ) {
	     
	       coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coupl[l+lmax+n*sizel+(-lprime+lmax+nprime*sizel)*size];
	     
	     } else if ( (l<0)&&(lprime<0) ) {
	     
	       coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coupl[-l+lmax+n*sizel+(-lprime+lmax+nprime*sizel)*size];
	     
	     }
	     
	   }
	 }
       }
     }
     /*time(&end);
     tick2=clock();
     printf("Elapsed time for calculation of coupling matrix: %ld ticks\n",tick2-tick1);	 
     printf("Elapsed time for calculation of coupling matrix: %.2lf seconds\n",difftime(end,start));*/


     // construction of the damper matrix
     //tick1 = clock();
     /*if (flagdamperfile==1) {
       dfactor=dampersum(nx, M, f0, tunefrac, d, freqd, nd, fcutoff);
     } else {
       dfactor=1.;
     }*/
     
     for(l=lmax; l>=-lmax; l--) {

       for(lprime=lmax; lprime>=-lmax; lprime--) {
       
         coefj=std::pow(jimag,lprime-l);

	 for(n=0; n<=nmax; n++) {
	 
	   fact=1./((double)factorialpart(n+abs(l),n)*std::pow(2.,abs(l))); // =factorial(n)/(factorial(n+|l|)*2^|l|)
	   
	   for(nprime=0; nprime<=nmax; nprime++) {

	     if ( (l>=0)&&(lprime>=0) ) {
	     
	       if ( ( (l>lmaxold)||(lprime>lmaxold) ) || ( (n>nmaxold)||(nprime>nmaxold) ) ) {
	       
	         // matrix coefficient was never computed for this l, lprime, n, nprime
		 
		 if (flagdamperfile==1) {
		   // damper from an impedance-like function
		   damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coefj*fact*dampersum(nx, M, omegaksi, omega0, tunefrac, l, lprime, n, nprime, aovertaub2, bovertaub2, taub, g, ng, d, freqd, nd, abseps);

		 } else {
		   // bunch-by-bunch damper
		   I_lprimenprime=Iln(lprime,nprime,-omegaksi,aovertaub2,bovertaub2,taub);
		   G_ln=Gln(l,n,-omegaksi,aovertaub2,taub,ng,g);
		   damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coefj*fact*I_lprimenprime*G_ln;
		 }
		 //printf("l=%d, l'=%d, n=%d, n'=%d, damper matrix: %13.8e %13.8e\n",l,lprime,n,nprime,
		 //	damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].real(),damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].imag());
		 
	       } else {
	       
	         // matrix coefficient was already computed
		 damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size] = 
		 	damperold[l+lmaxold+n*sizelold+(lprime+lmaxold+nprime*sizelold)*sizeold];
	       
	       }
	       
	     } else if ( (l<0)&&(lprime>=0) ) {
	     
	       damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=damper[-l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];
	     
	     } else if ( (l>=0)&&(lprime<0) ) {
	     
	       damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=damper[l+lmax+n*sizel+(-lprime+lmax+nprime*sizel)*size];
	     
	     } else if ( (l<0)&&(lprime<0) ) {
	     
	       damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=damper[-l+lmax+n*sizel+(-lprime+lmax+nprime*sizel)*size];
	     
	     }
	     
	   }
	 }
       }
     }
     //tick2=clock();
     //printf("Elapsed time for calculation of damper matrix: %ld ticks\n",tick2-tick1);	 


     // construction of the final matrix, and old coupling and damper matrices
     
     // normalization factor for damper matrix (at zero chromaticity)
     if (flagdamperfile==1) {
       dfactor=2.*pi*dampersum(nx, M, 0., omega0, tunefrac, 0, 0, 0, 0, aovertaub2, bovertaub2, taub, g, ng, d, freqd, nd, abseps);
     } else {
       I_lprimenprime=Iln(0,0,0.,aovertaub2,bovertaub2,taub);
       G_ln=Gln(0,0,0.,aovertaub2,taub,ng,g);
       dfactor=2.*pi*I_lprimenprime*G_ln;
     }
     // if normalize at the current chromaticity:
     //dfactor=2.*pi*damper[0+lmax+0*sizel+(0+lmax+0*sizel)*size];
     printf("damper matrix norm. factor = %13.8e %13.8e\n",dfactor.real(),dfactor.imag());
     
     for(l=lmax; l>=-lmax; l--) {

       for(lprime=lmax; lprime>=-lmax; lprime--) {
       
	 for(n=0; n<=nmax; n++) {
	 
	   for(nprime=0; nprime<=nmax; nprime++) {

	     Zsum  = (jimag*f0*dmax*2.*bovertaub2/(g[0]*dfactor))*std::exp(jimag*dphase)*damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];
	     Zsum += jimag*Ib*e/(2.*gamma*m0*clight*tune)*coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];
	     if ( (l==lprime)&&(n==nprime) ) Zsum += omegas*(double)l;
	     Mat[(l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size)*2]=Zsum.real();
	     Mat[(l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size)*2+1]=Zsum.imag();

	     damperold[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=damper[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];
	     couplold[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size]=coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size];

	     }
	   }
	 }
       }

       // diagonalization of the matrix Mat using LAPACK routine
       //tick1 = clock();
       zgeev_(&left,&right,&size,Mat,&size,eigenval,vl,&ldvl,eigenvect,&size,work,&lwork,rwork,&info); 
       //printf("%13.8e %13.8e %d %d\n",work[0].real(),work[0].imag(),lwork,info);
       //tick2 = clock();
       //printf("Elapsed time for calculation and diagonalization of one matrix: %ld ticks\n",tick2-tick1);

       // computes the imaginary part of the eigenvalues
       for(p=0; p<size; p++) {
	 eigenvalimag[p]=eigenval[p].imag();
	 //printf("%13.8e %13.8e %13.8e\n",eigenval[p].real(),eigenval[p].imag(),eigenvalimag[p]);
	 }
     
       // take the kmax lowest eigenvalues by imaginary parts
       gsl_sort_smallest_index (ind, kmax, eigenvalimag, 1, size); // ind gives the sorted indices

       // test convergence with respect to matrix size
       err=0.; // maximum relative error for the kmax first eigenvalues
       for(k=0; k<=kmax-1; k++) {
         if (std::abs((eigenval[ind[k]].imag()-eigenvalold[k].imag()))/(std::abs(eigenval[ind[k]].imag())+10.*abseps)>err) {
	   err=std::abs((eigenval[ind[k]].imag()-eigenvalold[k].imag()))/(std::abs(eigenval[ind[k]].imag())+10.*abseps);
	   kerr=k;
	 }
         //err=std::max(err,std::abs(eigenval[ind[k]]-eigenvalold[k])/std::abs(eigenval[ind[k]]));
         eigenvalold[k]=eigenval[ind[k]];
       }
       if ( (nbrad==0)&&(max_azim==0) ) flagcrit=(err >= crit);
       else flagcrit=( (nmax<(nbrad-1))||(lmax<max_azim) );
       
       // extend the matrix
       lmaxold=lmax;sizeold=size;sizelold=sizel;lmax++;nmax+=1;
       
       printf("error=%13.8e, lmaxold=%d, lmaxnew=%d, nmaxnew=%d, oldsizel=%d, oldsize=%d\n",err,lmaxold,lmax,nmax,sizel,size);
     } while (flagcrit);
     //tick2 = clock();
     //printf("Elapsed time for convergence of eigenvalues: %ld ticks\n",tick2-tick1);
     
     lmax=lmaxold;sizel=2*lmax+1;nmax-=1;size=sizel*nmax;

     // print the sorted eigenvalues and eigenvectors
     fprintf(fileigval,"Coupled-bunch mode %d, max. headtail mode %d, nb radial modes %d, max relative error on imag. parts %13.8e (note: absolute precision cannot be less than %lf)\n",nx,lmax,nmax+1,err,10.*abseps);
     fprintf(fileigvect,"Coupled-bunch mode %d, max. headtail mode %d, nb radial modes %d, max relative error on imag. parts %13.8e (note: absolute precision cannot be less than %lf)\n",nx,lmax,nmax+1,err,10.*abseps);
     for(k=0; k<=kmax-1; k++) {
       fprintf(fileigval,"%20.15e %20.15e\n",eigenval[ind[k]].real(),eigenval[ind[k]].imag());
       for(l=-lmax; l<=lmax; l++) {
         for(n=0; n<=nmax; n++) {
           fprintf(fileigvect,"%20.15e %20.15e\n",
	 		eigenvect[l+lmax+n*sizel+ind[k]*size].real(),
	 		eigenvect[l+lmax+n*sizel+ind[k]*size].imag());
	 }
         fprintf(fileigvect,"\n");
       }
       fprintf(fileigvect,"\n\n");
     }
     fprintf(fileigval,"\n");
     fprintf(fileigvect,"\n");
     
     // print the coupling and damper matrices
     /*fprintf(filcoupl,"Coupled-bunch mode %ld, max. headtail mode %d, nb radial modes %d\n",nx,lmax,nmax);
     fprintf(fildamp,"Coupled-bunch mode %ld, max. headtail mode %d, nb radial modes %d\n",nx,lmax,nmax);
     for(l=-lmax; l<=lmax; l++) {
       for(n=0; n<=nmax; n++) {
         for(lprime=-lmax; lprime<=lmax; lprime++) {
	   for(nprime=0; nprime<=nmax; nprime++) {
	     fprintf(filcoupl,"%20.15e %20.15e\n",coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].real(),
	     		coupl[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].imag());
	     fprintf(fildamp,"%20.15e %20.15e\n",damp[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].real(),
	     		damp[l+lmax+n*sizel+(lprime+lmax+nprime*sizel)*size].imag());
	   }
	 }
         fprintf(filcoupl,"\n");
         fprintf(fildamp,"\n");
       }
     }
     fprintf(filcoupl,"\n");
     fprintf(fildamp,"\n");*/
     
     fprintf(fileigval,"\n");
     fprintf(fileigvect,"\n");
   }
 // end of coupled-bunch mode loop


 delete [] Z;
 delete [] freq;
 if (flagdamperfile==1) {
   delete [] d;
   delete [] freqd;
 }
 delete [] g;
 delete [] eigenvalold;
 delete [] ind;
 delete [] eigenval;
 delete [] eigenvalimag;
 delete [] eigenvect;
 delete [] Mat;
 delete [] coupl;
 delete [] couplold;
 delete [] damper;
 delete [] damperold;
 delete [] work;
 delete [] rwork;
 delete [] nxiscan;

 fclose(fileigval);
 fclose(fileigvect);
 /*fclose(filcoupl);
 fclose(fildamper);*/

 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);
 


}
