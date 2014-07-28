#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,re,dateutil,random,pytz
from numpy import *
from string import split, replace
from optparse import OptionParser
import matplotlib
from coll_settings_find_sigma import read_settings


def parsse():
    parser = OptionParser()
    parser.add_option("-o", "--output",help="Specify output filename",
                      default="coll.txt",dest="OUT")
    parser.add_option("-f", "--file",
                      help="Specify the settings file name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--bunchlength",type=float,
                      help="Specify the total bunch length (4 RMS) in ns",
                      metavar="BUNLEN", default=1.3,dest="BUNLEN")
    parser.add_option("-v", "--rfvoltage",type=float,
                      help="Specify the RF voltage in MV",
                      metavar="VRF", default=6.,dest="VRF")
    parser.add_option("-n", "--nbunch", type=int,
    		      help="Specify the number of bunches (default=1)",
		      metavar="NB",default=1,dest="NB")
    parser.add_option("-i", "--intensity",type=float,
                      help="Specify the intensity per bunch in number of particles",
                      metavar="INT", default=1.e11,dest="INT")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    return opt, args

    

if __name__ == "__main__":

    opt,args=parsse();
        
    e=1.602176487e-19; # elementary charge
    c=299792458; # speed of light
    sigz=opt.BUNLEN*c/4.e9; # bunch length RMS in m 
    vrf=opt.VRF*1.e6; # RF voltage in V
    gam=1.2254; # gamma(3/4)
    mu0=4e-7*pi;
    Z0=mu0*c;
	
    
    # read settings file
    name,mat,angle,length,halfgap,betax,betay=read_settings(opt.FILE);
 
    # open output file
    fileout=open(opt.OUT,'w');
    print >> fileout, "name\thalfgap[m]\tphase_shift[deg]"

    totalphase=0.;
    
    for k,hgap in enumerate(halfgap):
    
	# compute phase shift assuming cos(Phis)=1 and classic thick wall formula
	# with Gaussian bunch
	
	if mat[k]=='C': cond=1/1.5e-5;
	elif mat[k]=='CFC': cond=1/0.5e-5;
	elif mat[k]=='HBN': cond=1/2.5e-6;
	elif mat[k]=='CU': cond=1/17e-9;
	elif mat[k]=='W': cond=1/54e-9;
	
	phase=180.*opt.NB*opt.INT*e*c*gam*length[k]*sqrt(Z0/(2.*cond*sigz**3))/(hgap*vrf*4.*pi**3);
		
	print >> fileout, name[k], "\t%5.2e"% hgap, "\t%5.2e"% phase;
	
	totalphase+=phase;
	
	
    print "Total phase shift in deg. from this configuration (classic thick wall formula, Gaussian bunch): ", totalphase;
    
    fileout.close()
	
	
    sys.exit()

