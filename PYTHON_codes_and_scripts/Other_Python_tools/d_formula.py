#!/usr/bin/python2.6

import numpy as np
import sys

i=sys.argv[1]; # number of particles in 10^11 units
Qx=0.62; # fractional part of tune

lambd=np.exp(1j*2*np.pi*Qx-float(i)*0.265435);
f=(lambd**2-np.cos(2*np.pi*Qx)*lambd)/(np.cos(2*np.pi*Qx)*lambd-1);
eps=(1-f)/2; # damper gain for HEADTAIL

print round(eps.real*10000)/10000.;
