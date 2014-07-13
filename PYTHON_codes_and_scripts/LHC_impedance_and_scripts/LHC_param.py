#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);

from string import *
import numpy as np
from DELPHI import Qs_from_RF_param

def LHC_param(E0,E=7e12):

    # E is the energy in eV, V the RF voltage in V
    if (E==7e12): V=16e6;Estr='7TeV';taub=1.e-9; # full length in s
    elif (E==6.5e12): V=16e6;Estr='6p5TeV';taub=1.e-9; # full length in s
    elif (E==4e12): V=12e6;Estr='4TeV';taub=1.25e-9; # full length in s
    elif (E==3.5e12): V=12e6;Estr='3p5TeV';taub=1.2e-9; # full length in s
    elif (E==450e9): V=6e6;Estr='450GeV';taub=1.3e-9; # full length in s
    else:
    	print "LHC energy not recognized; take taub=1.25ns, V=12 MV";
	V=12e6;Estr=float_to_str(E/1e12)+'TeV';taub=1.25e-9; # full length in s

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='LHC';
    h=35640; # harmonic number
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    sigmaz=taub*beta*c/4.; # RMS bunch length (m)
    circ=26658.883; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    Qx=64.31;Qxfrac=Qx-np.floor(Qx);
    Qy=59.32;Qyfrac=Qy-np.floor(Qy);
    alphap=3.225e-4; # momentum compaction factor
    eta=alphap-1./(gamma*gamma); # slip factor
    
    # compute Qs from RF parameters
    Qs=Qs_from_RF_param(V,h,gamma,eta,phis=0.,particle='proton');
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    omegas=Qs*omega0;
    dphase=0.; # additional damper phase

    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,V,h;

