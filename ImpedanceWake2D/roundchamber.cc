/*
 *  roundchamber.cc 
 *  
 *  by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *  computes the impedance in a round chamber (see CERN note by N. Mounet and E. Metral, 
 "Electromagnetic fields created by a macroparticle in an infinitely long
and axisymmetric multilayer beam pipe", 2009, CERN-BE-2009-039)
 
In input : typical input file is

Machine:	LHC
Relativistic Gamma:	479.6
Impedance Length in m:	2.8
Number of layers:	3
Layer 1 inner radius in mm:	4
Layer 1 DC resistivity (Ohm.m):	4.3e-7
Layer 1 relaxation time for resistivity (ps):	0
Layer 1 real part of dielectric constant:	1
Layer 1 magnetic susceptibility:	0
Layer 1 relaxation frequency of permeability (MHz):	Infinity
Layer 1 thickness in mm:	0.003
Layer 2 DC resistivity (Ohm.m):	4e+12
Layer 2 relaxation time for resistivity (ps):	0
Layer 2 real part of dielectric constant:	4
Layer 2 magnetic susceptibility:	0
Layer 2 relaxation frequency of permeability (MHz):	Infinity
Layer 2 thickness in mm:	54
Layer 3 DC resistivity (Ohm.m):	7.2e-7
Layer 3 relaxation time for resistivity (ps):	0
Layer 3 real part of dielectric constant:	1
Layer 3 magnetic susceptibility:	0
Layer 3 relaxation frequency of permeability (MHz):	Infinity
Layer 3 thickness in mm:	Infinity
start frequency exponent (10^) in Hz:	-1
stop frequency exponent (10^) in Hz:	15
linear (1) or logarithmic (0) or both (2) frequency scan:	2
sampling frequency exponent (10^) in Hz (for linear):	8
Number of points per decade (for log):	100
when both, fmin of the refinement (in THz):	0.0008
when both, fmax of the refinement (in THz):	0.1
when both, number of points in the refinement:	40000
added frequencies [Hz]:	1e-9 1e-8 1e-7 1e-6 1e-5 0.0001 0.001 0.01 1e+16 
Yokoya factors long, xdip, ydip, xquad, yquad:	1 1 1 0 0
Comments for the output files names:	_some_element

The order of the lines can be whatever, but the exact sentences and the TAB before the parameter
indicated, are necessary. If there are more layers than indicated by the number of layers,
the additional one are ignored. The last layer is always assumed to go to infinity.

In output one gives five files with the impedances (longitudinal, x dipolar, y dipolar,
x quadrupolar, y quadrupolar). Each have 3 columns :frequency, real part and
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
#include <amp.h>
#include <ablas.h>
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


#define  Precision	160  // Number of bits for the mantissa of all numbers (classic doubles have 53 bits)
#define  MAXLINES	100  // Maximum number of lines in the input file (correspond to 5 layers in top and
			    // bottom parts)
#define  MAXCHAR	200  // Maximum number of characters in one line in the input file
#define  MAXCHARFILE	200  // Maximum number of characters in the output files name extension

const  amp::ampf<Precision> C=299792458;    // velocity of light [m/s]
const  amp::ampf<Precision> euler="0.577215664901532860606512090082402431042159335939923598805767234884867726777664670936947063"; // Euler's gamma constant


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
*** clog: Complex natural logarithm in multiprecision
***
***************************************************/

  amp::campf<Precision> clog(amp::campf<Precision> z){
  
    /* Complex natural logarithm of z in multiprecision. */
  
    amp::campf<Precision> result;
    amp::ampf<Precision> rho,phi;
    
    rho=amp::log(amp::abscomplex(z));
    phi=amp::atan2(z.y,z.x);
    
    result.x=rho;
    result.y=phi;
    
    return result;
  
  }
  

/**************************************************
*** bessel: compute modified Bessel functions I and K,
*** of order m and complex argument z, in multiprecision
*** 
*** This is better optimized for low values of the order m
***
***************************************************/

  void bessel(unsigned int m,amp::campf<Precision> z, 
    	amp::campf<Precision>& besi, amp::campf<Precision>& besk,
	amp::campf<Precision>& besinorm, amp::campf<Precision>& besknorm){
  
    // order m, argument z
    // outputs are besi=Im(z) and besk=Km(z), and normalized ones 
    // besinorm=Im(z)/exp(z) and besknorm=Km(z)*exp(z)

    // use Taylor series and asymptotic expansions from Abramowitz-Stegun. 
    
    amp::campf<Precision> power,sumi,sumk,zover2,z2,sumf,zover2m;
    unsigned int condition,factm=1;
    amp::ampf<Precision> fact,abspower,k,expoabs,absz2,psi,m1powm,m1powk,mu,mMP;
    //timeval c1,c2; // clock ticks
    
    //gettimeofday(&c1,0);
    
    // m1powm = (-1)^m
    if (m % 2 ==0) m1powm=1;
    else m1powm=-1;
    mMP=amp::ampf<Precision>(m); // m in multiprecision
      
    zover2=z/2; // z/2
    z2=amp::csqr(zover2); // z^2/4
    absz2=amp::abscomplex(zover2); // |z/2|
    
    if (absz2<30) {
      // use Taylor series for both I and K

      expoabs=amp::exp(amp::abscomplex(z2)); // exp(|z^2/4|)
      for (unsigned int i=1; i<=m; i++){
        factm *= i;
      }
      power.x=1/(amp::ampf<Precision>(factm));power.y=0; // power=1/m!
      psi=-2*euler; // -2*Eulergamma
      for (unsigned int i=1; i<=m; i++){
        psi += 1/amp::ampf<Precision>(i);
      }
      // psi=psi(1)+psi(m+1)
      sumi.x=0;sumi.y=0;sumk.x=0;sumk.y=0;
      k=0;
      
      // compute infinite sums in I and K (simultaneously)
      do {  
      	    // power=(z/2)^(2k)/(k! (m+k)!), psi=psi(k+1)+psi(m+k+1) (digamma function)
	    sumi += power;
	    sumk += psi*power; 
 	    k += 1;
	    fact= 1/(k*(mMP+k));
	    power *= z2*fact;
	    psi += (mMP+2*k)*fact;
	    abspower = amp::abscomplex(power)*expoabs;
	    condition = ( (abspower + sumi)!=sumi )&&( (abspower*(mMP+k) + sumk)!=sumk );
      	    } while ( condition ) ;
	    
      // compute (z/2)^m
      zover2m.x=1;zover2m.y=0;
      for (unsigned int i=1; i<=m; i++){
        zover2m *= zover2;
      }
      
      // another finite sum for Km(z)
      sumf.x=0;sumf.y=0;
      if (m!=0) {
        power.x=amp::ampf<Precision>(factm)/mMP;power.y=0; // =(m-1)!
        for (unsigned int i=0; i<=m-1; i++){
	  sumf += power;
          if (i==(m-2)) power *= -z2/amp::ampf<Precision>(i+1);
	  else power *= -z2/amp::ampf<Precision>((i+1)*(m-i-2));
        }
	sumf /= 2*zover2m;
      }
      
      // final results (unormalized)
      besi = zover2m*sumi;
      besk = sumf + m1powm*(-clog(zover2)*besi + zover2m*sumk/2);
      // normalized ones
      power=cexp(z);
      besinorm = besi/power;
      besknorm = besk*power;
      
    } else {
      // use asymptotic formulae for large |z| (note: accuracy control is less good)
      // Note: here z comes from the square root of a complex number so is of |argument| less
      // or equal to pi/2 (case when exactly equal to pi/2 is pathological)

      sumi.x=0;sumi.y=0;sumk.x=0;sumk.y=0;z2=8*z;
      k=0;m1powk=1;mu=amp::ampf<Precision>(4*m*m);
      power.x=1;power.y=0;
      
      // compute infinite sums in I and K (simultaneously)
      do {  
      	    // power=(prod_{n=1}^{n=k}(4m^2 - (2n-1)^2))/(k! (8z)^k), m1powk=(-1)^k
	    sumi += m1powk*power;
	    sumk += power; 
 	    k += 1;
	    power *= (mu-amp::sqr(2*k-1))/(k*z2);
	    m1powk = -m1powk;
	    abspower = 5*amp::abscomplex(power); // 5 is a security factor (no clue how much we should put exactly)
	    condition = ( (abspower + sumi)!=sumi )&&( (abspower + sumk)!=sumk );
      	    } while ( condition ) ;
      
      if (amp::atan2(z.y,z.x)>=amp::halfpi<Precision>()) std::cout << "Pb in bessel: argument of z=" <<
      		amp::atan2(z.y,z.x).toDec().c_str() << "\n";
      
      // final results (normalized)
      power=1/csqrt(amp::twopi<Precision>()*z);
      besinorm = sumi*power;
      besknorm = sumk*power*amp::pi<Precision>();
      // unormalized ones
      power=cexp(z);
      besi = besinorm*power;
      besk = besknorm/power;
      
    }
   
    //gettimeofday(&c2,0);
    //std::cout << "bessel: time= " << (c2.tv_sec-c1.tv_sec)*1.e6+(c2.tv_usec-c1.tv_usec) << " usec\n";
    //std::cout << m << " " << z.x.toDec().c_str() << " " << z.y.toDec().c_str() << "\n";
    /*std::cout << besi.x.toDec().c_str() << " " << besi.y.toDec().c_str() << "\n";
    std::cout << besk.x.toDec().c_str() << " " << besk.y.toDec().c_str() << "\n";
    std::cout << besinorm.x.toDec().c_str() << " " << besinorm.y.toDec().c_str() << "\n";
    std::cout << besknorm.x.toDec().c_str() << " " << besknorm.y.toDec().c_str() << "\n";*/
    
    return;
  
  }

  
/**************************************************
*** multilayerm: computes matrix for field matching (multiprecision)
***
***************************************************/

  void multilayerm(ap::template_2d_array< amp::campf<Precision> >& mat,
  	unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<Precision> >& eps1, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	amp::ampf<Precision> k, amp::ampf<Precision> beta){
	
    /* computes the matrix M for the multilayer field matching, from mu1, eps1, b of each of the N layers
    and the azimuthal mode number m, relativistic velocity factor beta
    and wave number k. */

    amp::campf<Precision> nup,nu2p,nup1,nu2p1,epsp1overnup1,epspovernup,mup1overnup1,mupovernup,xpp,xp1p;
    amp::campf<Precision> Impp,Imp1p,Kmpp,Kmp1p,Imppnorm,Imp1pnorm,Kmppnorm,Kmp1pnorm; // bessel functions of order m
    amp::campf<Precision> Im1pp,Im1p1p,Km1pp,Km1p1p,Im1ppnorm,Im1p1pnorm,Km1ppnorm,Km1p1pnorm; // bessel functions of order m-1
    amp::campf<Precision> Imppratio,Kmppratio,Imp1pratio,Kmp1pratio; // ratio of bessel functions
    amp::campf<Precision> ImIm,ImKm,KmIm,KmKm; // products of bessel functions
    amp::campf<Precision> tmp1;
    amp::ampf<Precision> mMP,beta2,k2;
    ap::template_2d_array< amp::campf<Precision> > mp1p;
    ap::template_2d_array< amp::campf<Precision> > matold;
   
    mp1p.setbounds(1,4,1,4);
    matold.setbounds(1,4,1,4);
    
    mMP=amp::ampf<Precision>(m);
  
    for (int i=1; i<=4; i++) {
      for (int j=1; j<=4; j++) {
        if (i==j) matold(i,j)=1;
	else matold(i,j)=0;
      }
    }
    
    mat=matold;

    beta2=amp::sqr(beta);
    k2=amp::sqr(k);
    // first layer
    nu2p=k2*(1-beta2*eps1(1)*mu1(1));
    nup=csqrt(nu2p);
    epspovernup=eps1(1)/nup;
    mupovernup=mu1(1)/nup;
    
    for (unsigned int p=1; p<=N-1; p++) {

      nu2p1=k2*(1-beta2*eps1(p+1)*mu1(p+1));
      nup1=csqrt(nu2p1);
      epsp1overnup1=eps1(p+1)/nup1;
      mup1overnup1=mu1(p+1)/nup1;
      
      xpp=nup*b(p);
      xp1p=nup1*b(p);
      
      bessel(m,xpp,Impp,Kmpp,Imppnorm,Kmppnorm);
      bessel(m,xp1p,Imp1p,Kmp1p,Imp1pnorm,Kmp1pnorm);
      if (m==0) {
        bessel(1,xpp,Im1pp,Km1pp,Im1ppnorm,Km1ppnorm);
        bessel(1,xp1p,Im1p1p,Km1p1p,Im1p1pnorm,Km1p1pnorm);
	Imppratio = Im1pp/Impp;
	Kmppratio =-Km1pp/Kmpp;
	Imp1pratio= Im1p1p/Imp1p;
	Kmp1pratio=-Km1p1p/Kmp1p;
      } else {
        bessel(m-1,xpp,Im1pp,Km1pp,Im1ppnorm,Km1ppnorm);
        bessel(m-1,xp1p,Im1p1p,Km1p1p,Im1p1pnorm,Km1p1pnorm);
	Imppratio = Im1pp/Impp - mMP/xpp;
	Kmppratio =-Km1pp/Kmpp - mMP/xpp;
	Imp1pratio= Im1p1p/Imp1p - mMP/xp1p;
	Kmp1pratio=-Km1p1p/Kmp1p - mMP/xp1p;
      }
      
      ImIm=Impp*Imp1p;
      ImKm=Impp*Kmp1p;
      KmIm=Kmpp*Imp1p;
      KmKm=Kmpp*Kmp1p;
      
      // submatrix P
      tmp1=-xp1p/epsp1overnup1;
      mp1p(1,1)=tmp1*ImKm*(epsp1overnup1*Kmp1pratio-epspovernup*Imppratio);
      mp1p(1,2)=tmp1*KmKm*(epsp1overnup1*Kmp1pratio-epspovernup*Kmppratio);
      mp1p(2,1)=tmp1*ImIm*(-epsp1overnup1*Imp1pratio+epspovernup*Imppratio);
      mp1p(2,2)=tmp1*KmIm*(-epsp1overnup1*Imp1pratio+epspovernup*Kmppratio);
      
      // submatrix Q
      if (m==0) {
	mp1p(1,3)=0;
	mp1p(1,4)=0;
	mp1p(2,3)=0;
	mp1p(2,4)=0;
      } else {
        tmp1=-(nu2p1/nu2p-1)*mMP/(beta*eps1(p+1));
	mp1p(1,3)=-tmp1*ImKm;
	mp1p(1,4)=-tmp1*KmKm;
	mp1p(2,3)=tmp1*ImIm;
	mp1p(2,4)=tmp1*KmIm;
      }
      
      // submatrix R
      tmp1=-xp1p/mup1overnup1;
      mp1p(3,3)=tmp1*ImKm*(mup1overnup1*Kmp1pratio-mupovernup*Imppratio);
      mp1p(3,4)=tmp1*KmKm*(mup1overnup1*Kmp1pratio-mupovernup*Kmppratio);
      mp1p(4,3)=tmp1*ImIm*(-mup1overnup1*Imp1pratio+mupovernup*Imppratio);
      mp1p(4,4)=tmp1*KmIm*(-mup1overnup1*Imp1pratio+mupovernup*Kmppratio);
      
      // submatrix Q
      if (m==0) {
	mp1p(3,1)=0;
	mp1p(3,2)=0;
	mp1p(4,1)=0;
	mp1p(4,2)=0;
      } else {
	tmp1=eps1(p+1)/mu1(p+1);
	mp1p(3,1)=tmp1*mp1p(1,3);
	mp1p(3,2)=tmp1*mp1p(1,4);
	mp1p(4,1)=tmp1*mp1p(2,3);
	mp1p(4,2)=tmp1*mp1p(2,4);
      }
      
      // matrix multiplication
      ablas::cmatrixgemm<Precision>(4,4,4,1,mp1p,1,1,0,matold,1,1,0,0,mat,1,1);
      
      nu2p=nu2p1;
      nup=nup1;
      epspovernup=epsp1overnup1;
      mupovernup=mup1overnup1;

      matold=mat;
    }
    
  
    return;
  }
  
  
/**************************************************
*** alphaTM: computes alphaTM (multiprecision)
***
***************************************************/

  std::complex<double> alphaTM(unsigned int N, unsigned int m,
	ap::template_1d_array< amp::campf<Precision> >& eps1, 
  	ap::template_1d_array< amp::campf<Precision> >& mu1, 
	ap::template_1d_array< amp::ampf<Precision> >& b,
	amp::ampf<Precision> k, amp::ampf<Precision> beta){
	
    /* function that computes alphaTM for a given azimuthal mode number m, from mu1, eps1, b 
    of each of the N layers, and from the relativistic velocity factor beta and wave number k. */
    
    amp::campf<Precision> alphaTM;
    std::complex <double> result;
    ap::template_2d_array< amp::campf<Precision> > mat; // 4*4 field matching matrix
    //timeval c1,c2; // clock ticks


    /* setting bounds for matrix mat */
    mat.setbounds(1,4,1,4);

    // compute the field matching 4*4 matrices (upper and lower layers)
    //gettimeofday(&c1,0);
    multilayerm(mat, N, m, eps1, mu1, b, k, beta);
    //gettimeofday(&c2,0);
    //std::cout << "multilayer: time= " << (c2.tv_sec-c1.tv_sec)*1.e6+(c2.tv_usec-c1.tv_usec) << " usec\n";

    // compute alphaTM
    alphaTM=(mat(1,2)*mat(3,3)-mat(3,2)*mat(1,3))/(mat(1,1)*mat(3,3)-mat(1,3)*mat(3,1));
    
    // conversion to double complex
    result=std::complex<double>(amp::ampf<Precision>(alphaTM.x).toDouble(),amp::ampf<Precision>(alphaTM.y).toDouble());

    return result;
    
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
 char Zxquadoutput[MAXCHARFILE+10],Zyquadoutput[MAXCHARFILE+10],Input[MAXCHARFILE+10];
 std::string data[MAXLINES],machine,commentoutput,dummy2;
 FILE *filZxdip, *filZydip, *filZxquad, *filZyquad, *filZlong, *filInput;
 unsigned int N,dummy0; // number of layers, then dummy parameter
 unsigned int m,n_input,n_added,nf,nfdup; /* azimuthal mode number, number of input lines,
 					number of individually added frequencies, total number of
					frequencies in the scan, number of duplicate frequencies */
 unsigned int typescan,nflog,nflin; /* type of frequency scan, number of freq. per decade, number of freq. in a lin. scan inside the log scan*/
 ap::template_1d_array< amp::campf<Precision> > eps1,mu1; /* eps1 and mu1 for the layers */
 ap::template_1d_array< amp::ampf<Precision> > b,thick; // position of boundaries; thickness of the layers 
 ap::template_1d_array< amp::ampf<Precision> > rho,tau,epsb,chi,fmu; /* layers
 						properties (DC resistivity, resistivity relaxation time,
						dielectric constant, magnetic susceptibility=mu_r-1,
						relaxation frequency of permeability)*/
 amp::ampf<Precision> omega,beta,k,gamma,dummy3; // parameters
 amp::ampf<Precision> mu0,eps0,Z0; // vacuum permeability and permittivity in SI units, and free space impedance
 amp::campf<Precision> jimagMP; // imaginary constant in multiprecision
 size_t limit=1000; // limit (number of intervals) for gsl integration algorithm
 double x,y,L,fminlog,fmaxlog,fminlin,fmaxlin,fsamplin,fadded[15],*freq,dif,dummy1,yokoya[5];
 std::complex<double> Zxdip,Zydip,Zxquad,Zyquad,Zlong,cst;
 std::complex<double> alphaTM0,alphaTM1,jimag; // alphaTM constants and imaginary constant
 size_t found;
 time_t start,end; // times
					
 // start time
 time(&start);
 
 // to test bessel functions
 /*amp::campf<Precision> zz,bes1,bes2,bes3,bes4;
 zz.x=1;zz.y=1;
 bessel(0,zz,bes1,bes2,bes3,bes4);bessel(0,zz*100000,bes1,bes2,bes3,bes4);
 bessel(2,zz,bes1,bes2,bes3,bes4);bessel(2,zz*100000,bes1,bes2,bes3,bes4);
 bessel(1,zz,bes1,bes2,bes3,bes4);bessel(1,zz*10,bes1,bes2,bes3,bes4);
 bessel(1,zz*100,bes1,bes2,bes3,bes4);
 zz.x="2.3";zz.y=0;
 bessel(0,zz,bes1,bes2,bes3,bes4);bessel(0,zz*100000,bes1,bes2,bes3,bes4);
 bessel(2,zz,bes1,bes2,bes3,bes4);bessel(2,zz*100000,bes1,bes2,bes3,bes4);
 bessel(1,zz,bes1,bes2,bes3,bes4);bessel(1,zz*10,bes1,bes2,bes3,bes4);
 bessel(1,zz*100,bes1,bes2,bes3,bes4);*/
 
 // constants
 mu0=4e-7*amp::pi<Precision>();
 eps0=1/(mu0*amp::sqr(C));
 Z0=mu0*C;
 jimagMP.x=0;jimagMP.y=1;
 jimag=std::complex<double>(0.,1.);
 
 // default values of the parameters (in case)
 typescan=0;
 fminlog=2;fmaxlog=13;nflog=10;n_added=0;
 fsamplin=8;fminlin=1;fmaxlin=2;nflin=100;nf=0;
 N=2;L=1.;gamma="479.6";
 
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
   read_input(data[i],"Number of layers",N,dummy1,dummy2,dummy3,0);
   //printf("yoyo %d %s %s %s %13.8e %d\n",i,data[i].c_str(),machine.c_str(),gamma.toDec().c_str(),L,N);
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
   if (data[i].find("Yokoya factors long, xdip, ydip, xquad, yquad") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos){
       yokoya[0]=strtod(data[i].substr(found+1,std::string::npos).c_str(),&endline);
       for (int k=1; k<=4; k++) {
         yokoya[k]=strtod(endline,&endline);
       }
     }
   }
 }
 //printf("%s %d \n",machine.c_str(),N);
 //printf("%13.8e %13.8e %d\n",fadded[1],fadded[2],n_added);
 //printf("%13.8e %13.8e %13.8e %13.8e %13.8e\n",yokoya[0],yokoya[1],yokoya[2],yokoya[3],yokoya[4]);
 //printf("%13.8e %13.8e %d %d\n",fminlog,fmaxlog,nflog,typescan);

 
 eps1.setbounds(1,N+1);mu1.setbounds(1,N+1);b.setbounds(1,N+1);thick.setbounds(1,N);
 rho.setbounds(1,N+1);tau.setbounds(1,N+1);epsb.setbounds(1,N+1);chi.setbounds(1,N+1);fmu.setbounds(1,N+1);

 // default values of the layers properties (in case)
 rho(2)="1e5";tau(2)=0;epsb(2)=1;chi(2)=0;fmu(2)="Infinity";b(1)=2;thick(1)="Infinity";

 // find inner radius(ii) of the chamber layers
 for (unsigned int i=1; i<=n_input; i++) {
   read_input(data[i],"Layer 1 inner radius in mm",dummy0,dummy1,dummy2,b(1),3);
 }
 //printf("%s \n",b(1).toDec().c_str());
 
 // find all the layers properties
 for (unsigned int i=1; i<=n_input; i++) {
   for (unsigned int p=1; p<=N; p++){
     read_input_layer(data[i],"DC resistivity",p,rho(p+1));
     read_input_layer(data[i],"relaxation time for resistivity (ps)",p,tau(p+1));
     read_input_layer(data[i],"real part of dielectric constant",p,epsb(p+1));
     read_input_layer(data[i],"magnetic susceptibility",p,chi(p+1));
     read_input_layer(data[i],"relaxation frequency of permeability",p,fmu(p+1));
     read_input_layer(data[i],"thickness in mm",p,thick(p));
   }
 }

 // units conversion to SI (fminlin and fmaxlin were in THz, tau was in ps, b in mm and fmu in MHz)
 fminlin*=1.e12;fmaxlin*=1.e12;
 b(1)=b(1)*1e-3;
 for (unsigned int p=2; p<=N+1; p++){
   tau(p)=tau(p)*1e-12;
   b(p)=thick(p-1)*1e-3+b(p-1);
   fmu(p)=fmu(p)*1e6;
   //printf("%s %s %s %s %s %s \n",rho(p).toDec().c_str(),tau(p).toDec().c_str(),epsb(p).toDec().c_str(),
   //	chi(p).toDec().c_str(),fmu(p).toDec().c_str(),b(p).toDec().c_str());
 }
  
 // first layer (inside the chamber) is always vacuum
 eps1(1)=1;mu1(1)=1;
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
 sprintf(output,"W%s_%dlayers%.2lfmm%s",machine.c_str(),N,
	1e3*double(amp::ampf<Precision>(b(1)).toDouble()),commentoutput.c_str());
 sprintf(Zxdipoutput,"Zxdip%s.dat",output);
 sprintf(Zydipoutput,"Zydip%s.dat",output);
 sprintf(Zxquadoutput,"Zxquad%s.dat",output);
 sprintf(Zyquadoutput,"Zyquad%s.dat",output);
 sprintf(Zlongoutput,"Zlong%s.dat",output);
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
 fprintf(filZxdip,"Frequency [Hz]\tRe(Zxdip) [Ohm/m]\tIm(Zxdip) [Ohm/m]\n");
 fprintf(filZydip,"Frequency [Hz]\tRe(Zydip) [Ohm/m]\tIm(Zydip) [Ohm/m]\n");
 fprintf(filZxquad,"Frequency [Hz]\tRe(Zxquad) [Ohm/m]\tIm(Zxquad) [Ohm/m]\n");
 fprintf(filZyquad,"Frequency [Hz]\tRe(Zyquad) [Ohm/m]\tIm(Zyquad) [Ohm/m]\n");
 fprintf(filZlong,"Frequency [Hz]\tRe(Zlong) [Ohm]\tIm(Zlong) [Ohm]\n");
 
 
 // impedance computation at each frequency: beginning of the loop
 for (unsigned int i=0; i<nf; i++) {
    
   omega=amp::twopi<Precision>()*freq[i];
   k=omega/(beta*C);
   
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
   //printf("%s \n",b(1).toDec().c_str());
   //printf("%s %s \n",omega.toDec().c_str(),k.toDec().c_str());
 
   // computes alphaTM0
   alphaTM0=alphaTM(N+1,0,eps1,mu1,b,k,beta);
   //printf("alphaTM0: %13.8e %13.8e\n",alphaTM0.real(),alphaTM0.imag());
 
   // computes alphaTM1
   alphaTM1=alphaTM(N+1,1,eps1,mu1,b,k,beta);
   //printf("alphaTM1: %13.8e %13.8e\n",alphaTM1.real(),alphaTM1.imag());
 
   
   // computes and writes the impedances (applying yokoya factors)
   cst=jimag*L*double(amp::ampf<Precision>(k*Z0/(beta*amp::sqr(gamma)*amp::twopi<Precision>())).toDouble());
   Zlong=cst*alphaTM0*yokoya[0];
   cst=cst*double(amp::ampf<Precision>(k/(2*amp::sqr(gamma))).toDouble());
   Zxdip=cst*alphaTM1*yokoya[1];
   Zydip=cst*alphaTM1*yokoya[2];
   if ( (yokoya[0]==1)&& ( (yokoya[1]==1)&&(yokoya[2]==1) ) ) {
     // case of axisymmetric geometry -> small quadrupolar impedance (see ref. cited at the beginning)
     Zxquad=cst*alphaTM0;
     Zyquad=Zxquad;
   } else {     
     Zxquad=cst*alphaTM1*yokoya[3];
     Zyquad=cst*alphaTM1*yokoya[4];
   }
   fprintf(filZlong,"%13.8e %13.8e %13.8e\n",freq[i],Zlong.real(),Zlong.imag());
   fprintf(filZxdip,"%13.8e %13.8e %13.8e\n",freq[i],Zxdip.real(),Zxdip.imag());
   fprintf(filZydip,"%13.8e %13.8e %13.8e\n",freq[i],Zydip.real(),Zydip.imag());
   fprintf(filZxquad,"%13.8e %13.8e %13.8e\n",freq[i],Zxquad.real(),Zxquad.imag());
   fprintf(filZyquad,"%13.8e %13.8e %13.8e\n",freq[i],Zyquad.real(),Zyquad.imag());
  
   }
   // end of loop on frequencies
 
  
 fclose(filZxdip);
 fclose(filZydip);
 fclose(filZxquad);
 fclose(filZyquad);
 fclose(filZlong);
 
 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);

}
