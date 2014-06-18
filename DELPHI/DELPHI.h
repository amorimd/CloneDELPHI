#ifndef DELPHI_H
#define DELPHI_H

#include <complex>

using std::complex;

extern "C"
{

  long locate (double *table, double z, long n);

  complex<double> impedance (double *Z, double *freq, double omega, long nZ);

  double damper_gain (double *d, double *freq, double f, long nd, double fcutoff);

  long factorial(int n);

  long factorialpart(int n,int p);

  int minus1pow(int n);

  double Laguerre(int n, int k, double x);

  double Gln(int l, int n, double omega, double a, double taub, long ng, double *g);

  double Iln(int l, int n, double omega, double a, double b, double taub);

  complex<double> integrand(double omegap, double omegaksi, int l, int lprime, int n, int nprime, double a,
	  double b, double taub, double *g, long ng, double *Z, double *freq, long nZ);

  double integrand_real(double omegap, void *p);

  double integrand_imag(double omegap, void *p);

  double integrand_real_modif(double t, void *p);

  double integrand_imag_modif(double t, void *p);

  complex<double> impedanceint(double omega1, double omega2, int l, int lprime, int n, int nprime, double a, double b, double taub,
	  double omegaksi, double *g, long ng, double *Z, double *freq, long nZ);

  complex<double> impedancetrapz(double omega1, double omega2, int l, int lprime, int n, int nprime, double a, double b, double taub,
	  double omegaksi, double *g, long ng, double *Z, double *freq, long nZ);

  complex<double> impedancesum(int nx, int M, double omegaksi, double omega0, double Q, int l, int lprime,
	  int n, int nprime, double a, double b, double taub, double *g, long ng, double *Z, double *freq,
	  long nZ, int flag_trapz, double abseps);

  complex<double> dampersum(int nx, int M, double omegaksi, double omega0, double Q, int l, int lprime,
	  int n, int nprime, double a, double b, double taub, double *g, long ng, 
	  double *d, double *freq, long nd, double abseps);

  void read_input(std::string line, std::string description, int& param0, double& param1, 
  	  char *param2, int type);

  long read_nZ_impedance(char *filename);

  void read_impedance(char *filename, long& nZ, double *freq, double *Z);

  long read_nd_damper(char *filename);

  void read_damper(char *filename, long& nd, double *freqd, double *d);

  long read_ng_longdist(char *filename);

  void read_longdist(char *filename, long& ng, double *g);
 
}

#endif
