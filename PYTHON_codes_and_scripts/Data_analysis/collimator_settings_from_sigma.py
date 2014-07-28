#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from read_ADT_fit import build_colors;
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,re,dateutil,random,pytz
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from optparse import OptionParser
from Timber import parseout
from plot_lib import init_figure,end_figure
from io_lib import read_ncol_file_identify_header
from string_lib import find_ind_names
from collimator_settings import compute_onesigma_halfgap,read_ref_coll_file
import matplotlib
import matplotlib.dates


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--emittance",type=float,
                      help="Specify the normalized emittance (in mm.mrad) in x and y. Default=3.5 mm.mrad.",
                      metavar="EMIT", default=3.5,dest="EMIT")
    parser.add_option("-b", "--betafile",
                      help="Specify the beta functions file name. Default=None (beta functions provided in ref. file from -r option).",
                      metavar="BETAFILE", default=None,dest="BETAFILE")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify the energy in GeV. Default= 450 GeV",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-i", "--inj",action="store_false",
                      help="Specify if half-gaps of injection protection collimators are computed from the sigmas (default=take the same half-gaps as in the sigma file)",
                      metavar="INJ", default=True,dest="INJ")
    parser.add_option("-o", "--output",help="Specify output filename",
                      default="coll.txt",dest="OUT")
    parser.add_option("-r", "--reffile",
                      help="Specify the reference file name (where coll. names, beta functions and angles can be found)",
                      metavar="RFILE", default=None,dest="RFILE")
    parser.add_option("-s", "--sigmafile",
                      help="Specify the file with the halfgaps and sigmas",
                      metavar="SFILE", default=None,dest="SFILE")
    (opt, args) = parser.parse_args()
    print "Selected Sigma File:", opt.SFILE
    print "Selected Reference file:", opt.RFILE
    return opt, args

def read_sigma_coll_file(filename):

    # read collimator settings file with number of sigmas indicated
    nsigma=[];halfgap=[];collname=[];
    fid=open(filename,"r")
    
    for j,l in enumerate(fid.readlines()):
	ll=l.strip().split();
	#sys.stdout.write(l);sys.stdout.flush();
    	if (j==0):
	    # find column name "name"
	    for k,col in enumerate(ll):
	    	if col.startswith('name'): namecol=k;
	    	if col.startswith('halfgap'): hgcol=k;
	    	if col.startswith('nsigma'): nsigcol=k;
    	else:
	    collname.append(ll[namecol]);
	    halfgap.append(float(ll[hgcol]));
	    nsigma.append(float(ll[nsigcol]));
    fid.close()
    
    return collname,halfgap,nsigma;


if __name__ == "__main__":
    opt,args=parsse();

    
    if (opt.RFILE==None): print "Specify reference file for the collimator settings!"; sys.exit();

    if (opt.SFILE==None): print "Specify nsigma file for the collimator settings!"; sys.exit();
    
    # read reference file
    collname,angle,betax,betay=read_ref_coll_file(opt.RFILE);
    
    if opt.BETAFILE!=None:
	# different file with beta functions
	names=read_ncol_file_identify_header(opt.BETAFILE,'[nN]ame');
	betax=read_ncol_file_identify_header(opt.BETAFILE,'[bB]etax');
	betay=read_ncol_file_identify_header(opt.BETAFILE,'[bB]etay');
	# reorder such that the coll. names match with namesref
	ind=find_ind_names(collname,names);
	betax=betax[ind];betay=betay[ind];
    
    # read nsigma file
    collnamesig,halfgap,nsigma=read_sigma_coll_file(opt.SFILE);

    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    
    # open output settings file
    file=open(opt.OUT,'w');
    print >> file, "name\thalfgap[m]\tnsigma"

    for k,name in enumerate(collname):
    
	# compute 1 sigma half-gap (collimator characteristics, energy and emittance)
	gap1sig=compute_onesigma_halfgap(betax[k],betay[k],angle[k],gamma,opt.EMIT*1.e-6)
	    
	for knam,nam in enumerate(collnamesig):
	    if nam.startswith(name):
	    	nsig=nsigma[knam];hg2=halfgap[knam];
	
	if ((not(name.startswith('TDI')))and(not(name.startswith('TCLI'))))or(not(opt.INJ)):
	    hg=nsig*gap1sig;
	else:
	    hg=hg2;

	print >> file, name, "\t", hg, "\t", nsig;
	#print name,dt,hg,halfgap,gu[0]
	
	
    file.close()
	
    sys.exit()

