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



def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--batch",type=int,
                      help="Specify the number of bunches in each elementary batch", 
                      metavar="BATCH", default=72,dest="BATCH")
    parser.add_option("-e", "--empty",type=int,action="append",
                      help="Specify number of empty buckets after each train (25ns buckets). Use either -f or -e options.", 
                      metavar="EMPTY", default=None,dest="EMPTY")
    parser.add_option("-f", "--first",type=int,action="append",
                      help="Specify first bucket number of each train (2.5ns buckets, as in injection scheme viewer). Use either -f or -e options.", 
                      metavar="FNUM", default=None,dest="FNUM")
    parser.add_option("-i", "--firstbatch",type=int,
                      help="Specify number of bunches in the first batch (considered separately)", 
                      metavar="FBATCH", default=None,dest="FBATCH")
    parser.add_option("-m", "--morempty",action="store_true",
                      help="Specify if there are more empty 25ns slots between batches: instead of 8 (default), it becomes 8 + nb of empty slots from the spacing (i.e. 1 if 50ns, 2 if 75ns, etc.).",
                      metavar="MORE", default=False,dest="MORE")
    parser.add_option("-o", "--output",
                      help="Specify the output file name (without extension)",
                      metavar="OUT", default="out",dest="OUT")
    parser.add_option("-s", "--spacing",type=int,
                      help="Specify the spacing in ns (25, 50, 75, etc.)",
                      metavar="SPAC", default=25,dest="SPAC")
    parser.add_option("-t", "--train",type=int,action="append",
                      help="Specify the number of elementary batches in each train (i.e. at each injection), spaced by 8 empty 25ns buckets. You can specify it for each -f or -e option, except for the first one when -i option used", 
                      metavar="TRAIN", default=None,dest="TRAIN")
    (opt, args) = parser.parse_args()
    return opt, args


def wrt_batch(file,n,sp,count):
    # writes n consecutive bunches spaced by sp ns to the file in input, in 
    # terms of 25ns buckets, and increments the 25ns bucket counter count
    for i in range(0,n):
    	print >> file, 1;
	for l in range(0,sp/25-1): print >> file, 0;
    count=count+n*sp/25;
    return count;
    
    
def check_nbunch(file):
    # check the total number of bunches
    nbunch=0;nline=0;
    for l in file.readlines(): 
    	nbunch=nbunch+(l.strip()!='0');
	nline=nline+1;
    return nbunch,nline

                  
if __name__ == "__main__":
    opt,args=parsse();
    
    if (len(opt.TRAIN)==1):
        m=opt.TRAIN[0];
	if (opt.FNUM!=None):
    	    for j in range(1,len(opt.FNUM)): opt.TRAIN.append(m);
	elif (opt.EMPTY!=None):
    	    for j in range(1,len(opt.EMPTY)): opt.TRAIN.append(m);

    if (opt.FBATCH!=None):
        n=len(opt.TRAIN)
        for j,t in enumerate(opt.TRAIN[::-1]):
	    if (j==0): opt.TRAIN.append(t);
	    else: opt.TRAIN[n-j]=t;
        opt.TRAIN[0]=1;


    # parameter to subtract to the number of empty slots between batches
    if opt.MORE: parempty=0;
    else: parempty=(opt.SPAC/25-1);
    
    #print opt.TRAIN
    
    file=open(opt.OUT+'.bunch','w');
    fnumold=0;
    if (opt.FNUM!=None):
	for j,fnum in enumerate(opt.FNUM):
            if (j==0):
	        if (opt.FBATCH!=None): nbatch=opt.FBATCH;
		else: nbatch=opt.BATCH;
		for k in range(0,opt.TRAIN[0]-1):
		    fnumold=wrt_batch(file,nbatch,opt.SPAC,fnumold);
		    # write 8 empty buckets between elementary batches
		    for i in range(0,8-parempty): print >> file, 0;
		    fnumold=fnumold+8-parempty;
		fnumold=wrt_batch(file,nbatch,opt.SPAC,fnumold);
	    else:
		# fills with 0 between trains
		for i in range(fnumold,(fnum-opt.FNUM[0])/10):
	    	    print >> file, 0;
		    fnumold=fnumold+1;
		# now the train
		for k in range(0,opt.TRAIN[j]-1):
		    fnumold=wrt_batch(file,opt.BATCH,opt.SPAC,fnumold);
		    # writes 8 empty buckets between elementary batches
		    for i in range(0,8-parempty): print >> file, 0;
		    fnumold=fnumold+8-parempty;
		fnumold=wrt_batch(file,opt.BATCH,opt.SPAC,fnumold);
		
    else:
	for j,nemp in enumerate(opt.EMPTY):
            if (j==0):
	        if (opt.FBATCH!=None): nbatch=opt.FBATCH;
		else: nbatch=opt.BATCH;
		for k in range(0,opt.TRAIN[0]):
		    #print k,opt.TRAIN[0],j,nemp
		    fnumold=wrt_batch(file,nbatch,opt.SPAC,fnumold);
		    # write 8 empty buckets between elementary batches
		    for i in range(0,8-parempty): print >> file, 0;
		    fnumold=fnumold+8-parempty;
		# fills with 0 between trains (empty 25ns buckets)
		for i in range(0,nemp):
	    	    print >> file, 0;
		    fnumold=fnumold+1;
	    else:
		# now the train
		for k in range(0,opt.TRAIN[j]):
		    #print k,opt.TRAIN[j],j,nemp
		    fnumold=wrt_batch(file,opt.BATCH,opt.SPAC,fnumold);
		    # writes 8 empty buckets between elementary batches
		    for i in range(0,8-parempty): print >> file, 0;
		    fnumold=fnumold+8-parempty;
		# fills with 0 between trains (empty 25ns buckets)
		for i in range(0,nemp):
	    	    print >> file, 0;
		    fnumold=fnumold+1;
    	
	if (j+1<len(opt.TRAIN)):
	    for k in range(0,opt.TRAIN[j+1]):
		#print k,opt.TRAIN[j+1],j+1
		fnumold=wrt_batch(file,opt.BATCH,opt.SPAC,fnumold);
		# writes 8 empty buckets between elementary batches
		for i in range(0,8-parempty): print >> file, 0;
		fnumold=fnumold+8-parempty;
	    
    
    file.close()
    
    file=open(opt.OUT+'.bunch','r');    
    n,m=check_nbunch(file);
    file.close()
    print "Number of non empty bunches = %d" % n
    print "Total number of 25ns buckets = %d" % m
    

    sys.exit()

