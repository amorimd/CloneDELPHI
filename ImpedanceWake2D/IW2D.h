#ifndef IW2D_H
#define IW2D_H

#include <complex>

using std::complex;

extern "C"
{

  unsigned long locate (double *table, double z, unsigned long n);

  complex<double> pchip(double z, double *zi, complex<double> *fi, complex<double> *di, unsigned long nz);
	
  complex<double> interp(double z, double *zi, complex<double> *fi, complex<double> *di, unsigned long nz, unsigned int interp_type);
	
  void pchipslope(complex<double>* d, double* x, complex<double> * y, unsigned long length);

  complex<long double> Phi(complex<long double> x,double eps);
    
  complex<long double> Psi(complex<long double> x,double eps);
    
  complex<long double> Lambda(complex<long double> x,double eps);
    
  complex<double> fourier_integral_inf(complex<double>* fi,complex<double>* d,complex<double> df, long double t, long double* omegai, long double* delta, unsigned long length, double eps, unsigned int* interp_type, int flaginf);

}

#endif
