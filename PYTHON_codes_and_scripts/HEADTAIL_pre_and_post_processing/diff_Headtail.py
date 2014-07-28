#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);

import numpy as np
from string import split, replace
from optparse import OptionParser
from tables_lib import diffshape,diff
import math
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail files root name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-t", "--threshold",type=float,
                      help="Specify above what threshold value we print the differences (default=0.)",
                      metavar="THR", default=0.,dest="THR")
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    #print "Selected bunches:", opt.BNUM
    return opt, args


def read(filename):
    # function to read data
    
    data=[];
    for line in open(filename):
    	d=[float(k) for k in split(line)];
    	if (len(d)!=0): data.append(d)
    # convert data to array
    return np.array(data);
    

if __name__ == "__main__":
    opt,args=parsse(); 

    # revolution frequency and period (assumes protons)
    #print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the LHC circumference
    Trev=1./frev;
    
    if (len(opt.FILE)!=2): print "Two files should be specified !";sys.exit();
    
    listoutput=['prt','hdtl','prb','bunchds','sample','trk'];
    
    # check if bpm files also present
    listbpmname=glob.glob(opt.FILE[0]+'_bpm*.dat');
    listbpmname.sort();
    for l in listbpmname:
    	j=l.find('bpm');
    	#print l[j:-4]
    	listoutput.append(l[j:-4]);
	

    # main loop: check the differences
    for i,datatype in enumerate(listoutput):
       
    	# read files
	data1=read(opt.FILE[0]+'_'+datatype+'.dat')
	data2=read(opt.FILE[1]+'_'+datatype+'.dat')
	
	if (diffshape(data1,data2)):
	    d1,d2,d3=diff(data1,data2);
	    if (d1>=opt.THR)or(d2>=opt.THR)or(d3>=opt.THR):
	        print datatype+':'
		print "    rel (where data1 non zero):\t",d1;
		print "    abs (where data1 is zero):\t",d2;
		print "    rel over mean:\t\t",d3;
	else:
	    print datatype+':'
	    print "    shapes are not identical:\t",data1.shape,' ',data2.shape
    	    


    sys.exit()

