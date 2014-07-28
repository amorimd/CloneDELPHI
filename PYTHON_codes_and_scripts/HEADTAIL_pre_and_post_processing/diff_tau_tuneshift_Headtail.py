#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);

import numpy as np
from string import split, replace
from parser_lib import *
from tables_lib import diffshape
import math

def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) to compare on (several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail files root name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-m", "--mminmmax",type=int,nargs=2,
                      help="Specify minimum and maximum headtail mode to compare",
                      metavar="MMAX", default=[0,0],dest="MMAX")
    parser.add_option("-s", "--suffix",
    		      help="Specify suffix of the files",
                      default="tau.txt",dest="SUF")
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
    for j,line in enumerate(open(filename)):
        # skip first line
	if (j>0):
    	    d=[float(k) for k in split(line)];
    	    if (len(d)!=0): data.append(d)

    # convert data to array
    return np.array(data);
    

def diff(data1,data2):
    # compare two sets of data up to the line 'line' (not included)
    if (len(data1)!=0):
    	d=np.max(np.abs((data1[:,:5]-data2[:,:5])/data1[:,:5]),axis=0);
    	#d3=np.max(np.max(np.abs(data1-data2),axis=0)/np.abs(np.average(data1,axis=0))); # does not work (nan when one mean is zero))
    else: d=np.zeros(5);
    return d;


if __name__ == "__main__":

    print "\n"
    
    opt,args=parsse(); 

    mmax=opt.MMAX[1]
    mmin=opt.MMAX[0]
    
    if (opt.BNUM!=None):
	bunches=np.array(opt.BNUM);
	bunches.sort();
	print "Selected bunches:", bunches

    if (len(opt.FILE)!=2): print "Two files should be specified !";sys.exit();
    
   
    # main loop: check the differences
    for j,m in enumerate(range(mmin,mmax+1)):
       
    	# read files
	data1=read(opt.FILE[0]+'m'+str(m)+opt.SUF)
	data2=read(opt.FILE[1]+'m'+str(m)+opt.SUF)
	
	#if (diffshape(data1,data2)):
	
	if (opt.BNUM==None):
	    d=diff(data1,data2);
	    if (d[0]!=0): print "Not the same bunch number !!",d[0]
	else:
	    d=np.zeros(5)
	    for bnum in bunches:
		ind1=np.where(data1[:,0]==bnum);ind1=ind1[0];
		ind2=np.where(data2[:,0]==bnum);ind2=ind2[0];
		d1=diff(data1[ind1,:],data2[ind2,:])
		d=np.max([d,d1],axis=0)

	if (np.max(d)>=opt.THR):
	    print 'mode m='+str(m)+':'
	    print "    relative differences (taux, tauy, tunex, tuney): %.2g" % d[1],"%.2g" % d[2],"%.2g" % d[3],"%.2g" % d[4];
	#else:
	#    print "    shapes are not identical:\t",data1.shape,' ',data2.shape
    	    


    sys.exit()

