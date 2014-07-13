#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);

import numpy as np
from string import *
from string_lib import *


def SPS_param(E0,E=26e9,optics='Q26'):

    # E is the energy in eV

    e=1.602176487e-19; # elementary charge
    c=299792458;
    # fixed parameters
    machine='SPS';
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=6911; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    
    if optics=='Q26':
	Qx=26.13;
	Qy=26.18;
	alphap=1.9181e-3; # momentum compaction factor
	if (E==26e9): Qs=0.00725;Estr='26GeV';taub=2.6e-9; # full length in s. This corresponds to 3MV
	elif (E==450e9): Qs=0.00467;Estr='450GeV';taub=1.5e-9; # full length in s
	else:
    	    print "SPS energy not recognized; Q26 injection parameters taken";
	    Qs=7.25e-3;Estr=float_to_str(E/1e9)+'GeV';taub=2.6e-9; # full length in s

    elif optics=='Q20':
	Qx=20.13;
	Qy=20.18;
	alphap=3.1e-3; # momentum compaction factor
	if (E==26e9): Qs=0.01513;Estr='26GeV';taub=2.9e-9; # full length in s. This corresponds to 3MV
	elif (E==450e9): Qs=0.00595;Estr='450GeV';taub=1.6e-9; # full length in s
	else:
    	    print "SPS energy not recognized; Q20 injection parameters taken";
	    Qs=0.01513;Estr=float_to_str(E/1e9)+'GeV';taub=2.9e-9; # full length in s

    else: print 'SPS_param: only Q20 & Q26 optics implemented';sys.exit();
    
    sigmaz=taub*beta*c/4.; # RMS bunch length (m)
    Qxfrac=Qx-np.floor(Qx);
    Qyfrac=Qy-np.floor(Qy);
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    dphase=0.; # additional damper phase

    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr;
