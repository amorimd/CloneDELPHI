#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab
import numpy as np
from string import split, replace
from optparse import OptionParser
from construct_bunchtable_from_scheme import check_nbunch


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",
                      help="Specify the name of the FBCT file with slots and bunch intensities.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-o", "--output",
                      help="Specify the output file name (without extension).",
                      metavar="OUT", default="out",dest="OUT")
    parser.add_option("-t", "--threshold",type=float,
                      help="Specify intensity threshold (in 10^11 p+/b) below which bunches are not taken into account (default=0).",
                      metavar="THRES", default=0.,dest="THRES")
    (opt, args) = parser.parse_args()
    return opt, args


def read_file(filename,nheader=0):

    # read a file organized in columns
    # there can be a header (number of lines = nheader)
    data=[];
    
    for il,line in enumerate(open(filename)):
    	if (len(split(line))>0)and(il>=nheader):
	    data.append(np.array(split(line)))
	
    return data;
    

                  
if __name__ == "__main__":
    opt,args=parsse();
    

    file=open(opt.OUT+'.bunch','w');

    # read FBCT file
    data=read_file(opt.FILE,nheader=1);

    # first column is the 25 ns slot number
    slotFBCT=[int(data[i][0]) for i in range(len(data))];
    # second column is intensity in 10^11 p+/b
    intFBCT=[float(data[i][1]) for i in range(len(data))];
    print slotFBCT,intFBCT;
    
    # construct bunch table
    for slot in range(max(slotFBCT)+1):
    	
	if (slot in slotFBCT)and(intFBCT[slotFBCT.index(slot)]>opt.THRES): print >> file, intFBCT[slotFBCT.index(slot)];
	else: print >> file, 0;
	

    file.close()
    
    file=open(opt.OUT+'.bunch','r');    
    n,m=check_nbunch(file);
    file.close()
    print "Number of non empty bunches = %d" % n
    print "Total number of 25ns buckets = %d" % m
    

    sys.exit()

