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
from io_lib import *


def parsse():
    parser = OptionParser()
    parser.add_option("-o", "--output",help="Specify output filename (default=coll.txt)",
                      default="coll.txt",dest="OUT")
    parser.add_option("-f", "--file",
                      help="Specify the settings file name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify the momentum in GeV/c (default=450)",
                      metavar="EN", default=450,dest="EN")
    parser.add_option("-x", "--emitx",type=float,
                      help="Specify the horizontal emittance in mm.mrad (default=2)",
                      metavar="EPSX", default=2,dest="EPSX")
    parser.add_option("-y", "--emity",type=float,
                      help="Specify the vertical emittance in mm.mrad (default=2)",
                      metavar="EPSY", default=2,dest="EPSY")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    return opt, args


def read_settings(filename):

    # read a file with collimator settings
    # OBSOLETE
    
    name=[];mat=[];angle=[];halfgap=[];
    length=[];betax=[];betay=[];
    fid=open(filename,"r")
    for j,l in enumerate(fid.readlines()):
	ll=l.strip().split();
	#sys.stdout.write(l);sys.stdout.flush();
    	if (j==0):
	    # find column headers
	    for k,col in enumerate(ll):
	    	if col.startswith('name'): namecol=k;
		if col.startswith('Material'): matcol=k;
		if (col.startswith('Angle'))or(col.startswith('angle')): angcol=k;
		if col.startswith('Length'): lencol=k;
		if col.startswith('halfgap'): gapcol=k;
		if col.startswith('betax'): betaxcol=k;
		if col.startswith('betay'): betaycol=k;
    	else:
	    name.append(ll[namecol]);
	    mat.append(ll[matcol]);
	    angle.append(float(ll[angcol]));
	    length.append(float(ll[lencol]));
	    halfgap.append(float(ll[gapcol]));
	    betax.append(float(ll[betaxcol]));
	    betay.append(float(ll[betaycol]));
	    
	    
    fid.close();
    
    return name,mat,angle,length,halfgap,betax,betay;
    

if __name__ == "__main__":

    opt,args=parsse();
    
    # compute gamma
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=sqrt(1.-1./(gamma*gamma));
    print 'gamma: ',gamma;
    
    
    # read settings file
    #name,mat,angle,length,halfgap,betax,betay=read_settings(opt.FILE);
    name=read_ncol_file_identify_header(opt.FILE,'[nN]ame');
    angle=read_ncol_file_identify_header(opt.FILE,'[aA]ngle');
    halfgap=read_ncol_file_identify_header(opt.FILE,'[hH]alfgap');
    betax=read_ncol_file_identify_header(opt.FILE,'[bB]etax');
    betay=read_ncol_file_identify_header(opt.FILE,'[bB]etay');
    
 
    # open output file
    fileout=open(opt.OUT,'w');
    print >> fileout, "name\thalfgap[m]\tnsigma"

    
    for k,hgap in enumerate(halfgap):
    
	onesig=sqrt(1.e-6*(betax[k]*opt.EPSX*cos(angle[k])**2+betay[k]*opt.EPSY*sin(angle[k])**2)/beta/gamma)
	nsig=hgap/onesig;
	
	print >> fileout, name[k], "\t%5.2e"% hgap, "\t%5.2e"% nsig;
	
	
    fileout.close()
	
	
    sys.exit()

