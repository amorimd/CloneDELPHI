/* 
 *
 * laclare.cc
 *
 * by Nicolas Mounet (Nicolas.Mounet@cern.ch)
 *
 *
 * For a given impedance, solve Laclare's eigenvalue problem  in transverse, in the case of low intensity perturbations (no mode coupling), for different kinds of 
longitudinal distributions.
 
 
 It gives in output two files (name specified in input):
 - [output_filename].dat: one-line ascii file with plmax, the number of lines in the diagonalized matrices finally used
 - [output_filename]_val.dat: eigenvalues Deltaomega_cm (complex betatron 
 angular frequency shift of the modes m), array of size M*(2*mmax+1)*(2*plmax+1)
 (pth eigenvalue of the mth mode of the Mth multibunch mode). It is an ascii file
 with the array written in "blocks".
 - [output_filename]_vec.dat: the corresponding eigenvectors sigma_(x,m). This is 
 an array of size M*(2*mmax+1)*(2*plmax+1)*(2*plmax+1) (lth component of the pth 
 eigenvector of the mth mode of the Mth multibunch mode). It is an ascii file with the array written in "blocks".

 
Typical input file:
 

Impedance filename	Z.dat
Output filename	out
Total bunch length [seconds]	2.5e-9
Number of particles per bunch	5e10
Machine circumference [m]	6911
Relativistic gamma	27.7
Transverse tune	26.129
Synchrotron tune	7.25e-3
Longitudinal distribution	parabolicline
Type of particle	proton
Chromaticity (DeltaQ*p/Q*Deltap)	0.1
Momentum compaction factor	1.92e-3
Number of bunches	924
Maximum headtail mode considered	1
Maximum number of eigenvalues	5
Minimum coupled-bunch mode number to consider	800
Maximum coupled-bunch mode number to consider	10
 
 
The order of the lines can be whatever, but the exact sentences and the TAB before the parameter indicated, are necessary. 


Some explanatations for the input file:

 - Impedance filename: name of the file containing the dipolar impedance (3 columns, without header: frequency - NOT ANGULAR, real part of the impedance and imaginary part). It should be sorted in ascending order of frequencies. Frequencies should be positive.
 - Output filename: eigenvalues are put in [this_filename]_val.dat, eigenvectors in [this_filename]_vec.dat.
 - Longitudinal distribution: type of longitudinal distribution:
 waterbag (water bag bunch),parabolicline (parabolic line density - usually good for protons), parabolicamp (parabolic amplitude density), gaussiancut (gaussian with a cut at the bunch length - assumed to be at 4 sigmas, useful for comparison with Headtail, gaussian (gaussian with infinite tail - usually good for electrons) (see Laclare).
 - Type of particle: proton or electron.
 - Maximum headtail mode considered: mmax: modes considered are from -mmax to mmax.
 - Maximum number of eigenvalues: kmax: number of eigenvalues that are kept and accurate within 0.1%
 - Minimum coupled-bunch mode number to consider: nxmin: we compute modes of number between nxmin and M-1 (total number of bunches minus one)
 - Maximum coupled-bunch mode number to consider: nxmax: we also compute modes of number between 0 and nxmax

The rest is self-explanatory.

See Elias Metral's USPAS 2009 course : Bunched beams transverse coherent
 instabilities, and J.L. Laclare lectures at CAS (1987, p.264)
 
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
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_errno.h>

#define  MAXLINES	20  // Maximum number of lines in the input file
#define  MAXMAT		5000 // Maximum plmax 
#define  MAXSIZE	10001 // Maximum size of the matrix to diagonalize = 2*MAXMAT+1

using std::complex;

const double pi = 3.141592653589793238462643383279502884197; // pi
const complex<double> jimag = complex<double>(0.,1.); // imaginary unit

struct params {int m; int i; double taub; double omega1; double omega2;};

extern "C" void zgeev_(char *JOBVL, char *JOBVR, int *N, double *A, int *LDA, complex<double> *W,
	complex<double> *VL, int *LDVL, complex<double> *VR, int *LDVR, complex<double> *WORK,
	int *LWORK, double *RWORK, int *INFO);

/******************************************************************************
 *** function locate: search a table ordered in ascending order		    ***
 ***		      (inspired from Numerical Recipes) 		    ***
 *** 		      Gives the position lprov (integer)		    ***
 ***                  in table, such that table[lprov-1]<z<table[lprov]	    ***
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
 	printf ("\n\n    Need higher frequencies in the impedance file !\n");
	printf ("%lf",freq[nZ-1]);
      	exit(2); }
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
 *** function integrand: computes the term inside the integral in tauhat,   ***
 ***			for the	matrix elements	of the eigensystem.	     ***
 ***			This is in the appropriate format for gsl integration***
 *******************************************************************************/

double integrand(double tauhat, void *p)
  
{
	
 /* function computing the integrand, for gsl integration */

 struct params *param=(struct params *)p;
 int m,i;
 double inte,omega1,omega2,taub,J1,J2,g0,x;

 m=(param->m); // headtail mode number
 i=(param->i); // 0<=i<=4: indicates which longitudinal distribution
 omega1=(param->omega1); // angular frequency in front of tauhat in the first Bessel function
 omega2=(param->omega2); // angular frequency in front of tauhat in the second Bessel function
 taub=(param->taub); // bunch length in seconds
 
 J1=gsl_sf_bessel_Jn(m,omega1*tauhat);
 J2=gsl_sf_bessel_Jn(m,omega2*tauhat);
 
 switch (i) {

   case 0:
     g0=4.;
     break;

   case 1:
     x=2.*tauhat/taub;
     g0=6.*sqrt(1.-x*x);
     break;

   case 2:
     x=2.*tauhat/taub;
     g0=8.*(1.-x*x);
     break;

   case 3:
   case 4:
     x=2.*tauhat/taub;
     g0=8.*std::exp(-2.*x*x);
     break;

 }
 
 inte=g0*J1*J2*tauhat;
 //printf("%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",inte,x,g0,J1,J2,tauhat);
 return inte;
 
}
 
/*******************************************************************************
 *** function integrand_modif: computes the term inside the integral in	     ***
 ***			       tauhat, with the change of variable 	     ***
 ***			       tauhat=(1-t)/t), for the	matrix elements	of   ***
 ***			       the eigensystem.	Case of gaussian distribution***
 ***			      only (with infinite tails - i_dist=4).	     ***
 ***			This is in the appropriate format for gsl integration***
 *******************************************************************************/

double integrand_modif(double t, void *p)
  
{
	
 /* function computing the integrand, for gsl integration */

 struct params *param=(struct params *)p;
 int m,i;
 double inte,omega1,omega2,taub,J1,J2,g0,x,tauhat;

 m=(param->m); // headtail mode number
 i=(param->i); // 0<=i<=4: indicates which longitudinal distribution
 omega1=(param->omega1); // angular frequency in front of tauhat in the first Bessel function
 omega2=(param->omega2); // angular frequency in front of tauhat in the second Bessel function
 taub=(param->taub); // bunch length in seconds
 
 tauhat=(1.-t)/t;
 
 J1=gsl_sf_bessel_Jn(m,omega1*tauhat);
 J2=gsl_sf_bessel_Jn(m,omega2*tauhat);
 
 x=2.*tauhat/taub;
 g0=8.*std::exp(-2.*x*x);
 
 inte=g0*J1*J2*tauhat/(t*t);
 //printf("%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",inte,x,g0,J1,J2,tauhat);
 return inte;
 
}
 
/******************************************************************************
 *** 			main program					    ***
 ******************************************************************************/

main ()

{

 const double e=1.60218e-19; // elementary charge in C
 const double clight=299792458; // speed of light in m/s
 double beta,gamma,f0,Ib,Nb,omega0,omegaksi,eta,omegas,tune,tunes,tunefrac,circum,chroma,alphap,taub,m0,dif,err,tunefracxi;
 /*double *freq,*omega,*Z,dummy1,dummy2,dummy3,rwork[2*MAXSIZE],integ,omegap,omegal,x,jml,jmlprime,errx,
 	Kxm[2*MAXSIZE*MAXSIZE],Kxmold[2*MAXSIZE*MAXSIZE],eigenvalabs[MAXSIZE];
 complex<double> Zcomp,eigenval[MAXSIZE],*eigenvalold,eigenvect[MAXSIZE*MAXSIZE],work[40*MAXSIZE],vl[1];*/
 double *freq,*omega,*Z,dummy1,dummy2,dummy3,*rwork,integ,omegap,omegal,x,jml,jmlprime,errx,
 	*Kxm,*Kxmold,*eigenvalabs;
 complex<double> Zcomp,*eigenval,*eigenvalold,*eigenvect,*work,vl[1];
 char Zfilename[300],outputfilename[300],outputfilenamevect[300],outputfilenameval[300];
 char distribution[15],particle[10],right,left;
 std::string data[MAXLINES];
 FILE *filZ, *fileigval, *fileigvect;
 long i,j,n_input,nZ,nx,nxi,nxmin,nxmax;
 int i_part,i_dist,M,m,kmax,mmax,plmax,plmax0,plmaxold,k,l,p,info,size,sizeold,lwork,ldvl,condition,status;
 gsl_integration_workspace *w; // workspace for gsl adaptative integration
 gsl_function F;  // gsl function
 double tolintabs=1.e-15,tolint=1.e-6,pts[2]; // absolute and relative error permitted for gsl adaptative integration, and extrema of integration
 struct params param; // input parameters for the integrand functions (gsl integration)
 size_t limit=1000,found,*ind; // limit (number of intervals) for gsl integration algorithm
 time_t start,end; // times
 clock_t tick1,tick2; // number of clock ticks


 time(&start);

// extrema of integration in case of semi-infinite integral and change of variable
 pts[0]=0.;pts[1]=1.; 

 // deactivate abort on GSL errors
 gsl_set_error_handler_off();  
 
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
   if (data[i].find("Impedance filename") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(Zfilename,data[i].substr(found+1,std::string::npos).c_str());
     }            
   if (data[i].find("Output filename") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(outputfilename,data[i].substr(found+1,std::string::npos).c_str());
     }           
   if (data[i].find("Total bunch length") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) taub=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Number of particles per bunch") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) Nb=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Machine circumference") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) circum=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Relativistic gamma") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) gamma=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Transverse tune") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) tune=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Synchrotron tune") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) tunes=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Longitudinal distribution") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(distribution,data[i].substr(found+1,std::string::npos).c_str());
     }           
   if (data[i].find("Type of particle") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) strcpy(particle,data[i].substr(found+1,std::string::npos).c_str());
     }           
   if (data[i].find("Chromaticity") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) chroma=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Momentum compaction") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) alphap=strtod(data[i].substr(found+1,std::string::npos).c_str(),NULL);
     }
   if (data[i].find("Number of bunches") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) M=atoi(data[i].substr(found+1,std::string::npos).c_str());
     }
   if (data[i].find("Maximum headtail mode") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) mmax=atoi(data[i].substr(found+1,std::string::npos).c_str());
     }
   if (data[i].find("Maximum number of eigenvalues") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) kmax=atoi(data[i].substr(found+1,std::string::npos).c_str());
     }
   if (data[i].find("Minimum coupled-bunch mode number to consider") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) nxmin=atoi(data[i].substr(found+1,std::string::npos).c_str());
     }
   if (data[i].find("Maximum coupled-bunch mode number to consider") != std::string::npos){
     found=data[i].find("\t");
     if (found != std::string::npos) nxmax=atoi(data[i].substr(found+1,std::string::npos).c_str());
     }
 }
 
 //printf("%s %s %s %d %d %d\n",distribution,outputfilename,particle,M,mmax,plmax);
 //printf("%13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e %13.8e\n",tune,tunes,chroma,taub,gamma,circum,Nb,alphap);

 beta=sqrt(1.-1./(gamma*gamma)); // relativistic velocity factor
 f0=beta*clight/circum; // revolution frequency
 omega0=2.*pi*f0; // revolution angular frequency
 Ib=(double)M*e*Nb*f0; // bunch intensity
 eta=alphap-1./(gamma*gamma); // transition parameter (Elias Metral's convention, opposite to Joel Le Duff's one)
 omegaksi=chroma*tune*omega0/eta; // chromatic angular frequency
 omegas=tunes*omega0; // synchrotron angular frequency
 tunefrac=tune-floor(tune); // fractional part of the tune
 /* 10/06/2011: modification to include integer part of chromatic tune into the mode number 
 (matrix will be smaller -> faster computation): */
 // bug fix 28/10/2011: only a multiple of M should be included in to the mode number + rename to avoid using it in the impedance
 //tunefracxi=tunefrac-(chroma*tune/eta-(double)M*floor(chroma*tune/eta/(double)M));
 // new bug fix (Nov. 2011) : this apparently does not work (comparing with HEADTAIL) -> back to previous version before 10/06/2011
 tunefracxi=tunefrac-chroma*tune/eta;

 // default: proton and water bag distribution
 i_part=0;
 m0=1.6726e-27;
 i_dist=0;

 // beam particles' mass
 if (strcmp(particle,"proton")==0) {
	 i_part=0;
	 m0=1.6726e-27;
 }
 if (strcmp(particle,"electron")==0) {
	 i_part=1;
         m0=9.1094e-31;
 }

 // longitudinal distribution switch
 if (strcmp(distribution,"waterbag")==0) i_dist=0;
 if (strcmp(distribution,"parabolicline")==0) i_dist=1;
 if (strcmp(distribution,"parabolicamp")==0) i_dist=2;
 if (strcmp(distribution,"gaussiancut")==0) i_dist=3;
 if (strcmp(distribution,"gaussian")==0) i_dist=4;
 

 // read impedance file
 // first determines the number of lines nZ
 filZ = fopen(Zfilename,"r");
 if (filZ == NULL)
   { printf ("\n\n    File with impedance does not exist.\n") ;
   exit(2); }
 nZ = 0;
 do {condition = fscanf (filZ, "%lg %lg %lg\n",&dummy1,&dummy2,&dummy3);
   ++nZ;	 
 } while (condition != EOF);
 nZ--;
 fclose(filZ);
 // read frequencies and corresponding impedances
 filZ = fopen(Zfilename,"r");
 Z = new double[nZ*2];
 freq = new double[nZ];
 i=0;
 do {condition = fscanf (filZ, "%lg %lg %lg\n",&freq[i],&Z[2*i],&Z[2*i+1]);
 	if ( (i>0) && (freq[i] == freq[i-1]) ) {
	  // remove duplicate frequencies
	  i--;
	  nZ--;
	}
 	i++;
 } while (condition != EOF);
 fclose(filZ);
 //for (i=0;i<=nZ-1;i++) printf("%ld %ld %13.8e %13.8e %13.8e\n",nZ,i,freq[i],Z[2*i],Z[2*i+1]);

 // output filenames
 strcpy (outputfilenameval,outputfilename) ;
 strcpy (outputfilenamevect,outputfilename) ;
 strcat (outputfilenameval,"_val.dat") ;
 strcat (outputfilenamevect,"_vec.dat") ;
 // open files
 fileigval=fopen(outputfilenameval,"w");
 fileigvect=fopen(outputfilenamevect,"w");

 // workspace allocation for gsl adaptative integration
 w=gsl_integration_workspace_alloc(limit);
 
 /* note: compiling those arrays in static does not make any difference in term of computation time,
 but redefining them each time inside the loops makes a big one. */
 eigenvalold=new complex<double>[kmax];
 ind=new size_t[kmax];
 eigenval=new complex<double>[MAXSIZE];
 eigenvalabs=new double[MAXSIZE];
 eigenvect=new complex<double>[MAXSIZE*MAXSIZE];
 Kxm=new double[2*MAXSIZE*MAXSIZE];
 Kxmold=new double[2*MAXSIZE*MAXSIZE];
 work=new complex<double>[40*MAXSIZE];
 rwork=new double[2*MAXSIZE];

 // first estimation of plmax needed
 //plmax0=std::min(std::max(0,(int)ceil((20./taub+omegaksi)/((double)M*omega0)))+kmax,MAXMAT);
 //10/06/2011: modification
 plmax0=std::min(std::max(0,(int)ceil((20./taub)/((double)M*omega0)))+kmax,MAXMAT);
 printf("%d\n",plmax0);

 ldvl=1;left='N';right='V';
 
 for(nxi=nxmin; nxi<=M+nxmax; nxi++) {

   nx=nxi%M;
   // nx: coupled-bunch mode number   

   for(m=-mmax; m<=mmax; m++) {
     // consider each headtail mode m individually (low intensity, no coupling)
     // construction of the matrix K^(x,m)
     plmax=plmax0;plmaxold=-1;sizeold=-1;

     for (k=0; k<=kmax-1; k++) eigenvalold[k]=complex<double>(1.e50,1.e50);
     // loop for convergence on matrix size
     //tick1 = clock();
     do {
       //tick1 = clock();
       size=2*plmax+1;lwork=40*size; // note: not visible impact of lwork on computation time
       for(p=-plmax; p<=plmax; p++) {
	 omegap=((double)(nx+p*M)+tunefrac)*omega0+(double)m*omegas;
	 Zcomp=impedance(Z,freq,omegap,nZ);
	 //printf("%13.8e\n",omegap);

	 for(l=-plmax; l<=plmax; l++) {
	 
	   if ( (abs(p)<=plmaxold) && (abs(l)<=plmaxold) ) {
	     // if it was already computed before, useless to recompute it
	     Kxm[(l+plmax+(p+plmax)*size)*2]=Kxmold[(l+plmaxold+(p+plmaxold)*sizeold)*2];
	     Kxm[(l+plmax+(p+plmax)*size)*2+1]=Kxmold[(l+plmaxold+(p+plmaxold)*sizeold)*2+1];
	   }
	   else {
	     /*omegap=((double)(nx+p*M)+tunefrac)*omega0-omegaksi;
	     omegal=((double)(nx+l*M)+tunefrac)*omega0-omegaksi;*/
	     // 10/06/2011: modification
	     // 28/10/2011: bug fix (chromatic frequency)
	     omegap=((double)(nx+p*M)+tunefracxi)*omega0;
	     omegal=((double)(nx+l*M)+tunefracxi)*omega0;
	     // parameters for gsl integration
	     param.m=m;
	     param.i=i_dist;
	     param.taub=taub;
	     param.omega1=omegap;
	     param.omega2=omegal;
	     F.params=&param;
	     F.function=&integrand;

	     switch (i_dist) {

               case 0:
	       // case of water-bag distribution: integral can be computed analytically
               x=taub/2.;
               if (p==l) {
		 jml=gsl_sf_bessel_Jn(m,omegal*x);
		 jmlprime=-gsl_sf_bessel_Jn(m+1,omegal*x) + (double)m*jml/(omegal*x);
        	 integ=(x*x/2.)*jmlprime*jmlprime + 
	       		  (x*x-(double)(m*m)/(omegal*omegal))*jml*jml/2.;
               } else {
        	 integ=(x/(omegal*omegal-omegap*omegap))*(omegal*gsl_sf_bessel_Jn(m,omegap*x)* 
	       		  gsl_sf_bessel_Jn(m+1,omegal*x) - 
	       		  omegap*gsl_sf_bessel_Jn(m,omegal*x)*gsl_sf_bessel_Jn(m+1,omegap*x));
               }
	       integ=integ*4.;
	       break;

               case 1:
	       case 2:
	       case 3:
		 // case of distributions of finite extension (-taub/2 -> taub/2)
		 // gsl adaptative integration on [0,taub/2]
		 gsl_integration_qag(&F, 0., taub/2., tolintabs, tolint, limit, 
	     		    2, w, &integ, &errx); // key=2 -> Gauss-Kronrod with 21 points -> seems optimum
		 break;

	       case 4:
		 // case of infinite gaussian distribution
		 // gsl adaptative integration on [0,+Inf[
		 // Note: change of variable x=2.*tauhat/taub;
	     	 param.taub=2.;
		 param.omega1=omegap*taub/2.;
		 param.omega2=omegal*taub/2.;
	     	 F.params=&param;
		 status=gsl_integration_qagiu(&F, 0., tolintabs, tolint, limit, w, &integ, &errx);

		 if ( ( (status==GSL_EMAXITER)||(status==GSL_EDIVERGE) )||(status==GSL_EROUND) ) {
		   // try on finite interval with change of variable tauhat=(1-t)/t
		   printf("Semi-infinite integration failure - so now we try with change of variable instead");
		   F.function=&integrand_modif;
		   status=gsl_integration_qagp(&F, pts, 2, tolintabs, tolint, limit, w, &integ, &errx);
		   if ( (status)&&(std::abs(err/integ)>tolint) ) printf("Warning: integration: result= %13.8e, rel. error= %13.8e\n",integ,std::abs(errx/integ));
		 }
		 integ*=taub*taub/4.;

		 break;

       	     }

	     //printf("%13.8e %13.8e %13.8e %13.8e %d\n",integ,Zcomp.real(),Zcomp.imag(),errx,w->size);
	     // Matrix elements (indices order chosen to be compatible with Fortran LAPACK function ZGEEV)
	     Kxm[(l+plmax+(p+plmax)*size)*2]=Zcomp.real()*integ; //real part
	     Kxm[(l+plmax+(p+plmax)*size)*2+1]=Zcomp.imag()*integ; //imaginary part
	   }
	 }
       }
       for(p=-plmax; p<=plmax; p++) {
	 for(l=-plmax; l<=plmax; l++) {
	   Kxmold[(l+plmax+(p+plmax)*size)*2]=Kxm[(l+plmax+(p+plmax)*size)*2];
	   Kxmold[(l+plmax+(p+plmax)*size)*2+1]=Kxm[(l+plmax+(p+plmax)*size)*2+1];
	 }
       }

       // diagonalization of Kxm using LAPACK routine
       //tick1 = clock();
       zgeev_(&left,&right,&size,Kxm,&size,eigenval,vl,&ldvl,eigenvect,&size,work,&lwork,rwork,&info); 
       //printf("%13.8e %13.8e %d %d\n",work[0].real(),work[0].imag(),lwork,info);
       //tick2 = clock();
       //printf("Elapsed time for calculation and diagonalization of one matrix: %ld ticks\n",tick2-tick1);
       // apply constant multiplication factor and computes the absolute values of eigenvalues
       for(p=-plmax; p<=plmax; p++) {
	 eigenval[p+plmax]=eigenval[p+plmax]*jimag*Ib*e/(2.*gamma*m0*clight*tune*pi*taub*taub);
	 eigenvalabs[p+plmax]=std::abs(eigenval[p+plmax]);
	 //printf("%13.8e %13.8e %13.8e\n",eigenval[p+plmax].real(),eigenval[p+plmax].imag(),eigenvalabs[p+plmax]);
	 }
     
       // take the kmax highest eigenvalues by absolute values
       gsl_sort_largest_index (ind, kmax, eigenvalabs, 1, size); // ind gives the sorted indices
       // test convergence with respect to matrix size
       err=0.; // maximum relative error
       for(k=0; k<=kmax-1; k++) {
         err=std::max(err,std::abs(eigenval[ind[k]]-eigenvalold[k])/eigenvalabs[ind[k]]);
         eigenvalold[k]=eigenval[ind[k]];
       }
       plmaxold=plmax;sizeold=size;plmax=plmax+std::max(plmax0/10,10);
       printf("%13.8e %d\n",err,plmax);
     } while (err >= 5e-3);
     //tick2 = clock();
     //printf("Elapsed time for convergence of eigenvalues: %ld ticks\n",tick2-tick1);
     
     plmax=plmaxold;size=2*plmax+1;
     // print the sorted eigenvalues and eigenvectors
     for(k=0; k<=kmax-1; k++) {
       fprintf(fileigval,"%20.15e %20.15e\n",eigenval[ind[k]].real(),eigenval[ind[k]].imag());
       fprintf(fileigvect,"%d\n",plmax);
       for(l=-plmax; l<=plmax; l++) {
         fprintf(fileigvect,"%20.15e %20.15e\n",
	 		eigenvect[(l+plmax+ind[k]*size)].real(),
	 		eigenvect[(l+plmax+ind[k]*size)].imag());
       }
       fprintf(fileigvect,"\n");
     }

     fprintf(fileigval,"\n");
     fprintf(fileigvect,"\n");
   }
   fprintf(fileigval,"\n");
   fprintf(fileigvect,"\n");
 }


 gsl_integration_workspace_free(w);
 delete[] Z;
 delete[] freq;
 delete [] eigenvalold;
 delete [] ind;
 delete [] eigenval;
 delete [] eigenvalabs;
 delete [] eigenvect;
 delete [] Kxm;
 delete [] Kxmold;
 delete [] work;
 delete [] rwork;

 fclose(fileigval);
 fclose(fileigvect);

 time(&end);
 dif=difftime(end,start);
 printf("Elapsed time during calculation: %.2lf seconds\n",dif);
 

}
