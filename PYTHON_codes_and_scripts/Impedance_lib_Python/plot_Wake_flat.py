#!/usr/bin/python2.6

import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab
import numpy as np
import math,cmath
from string import split, replace
from optparse import OptionParser
from plot_lib import set_fontsize,init_figure,end_figure
from io_lib import list_files
from string_lib import takeout_common


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the file name (longitudinal wake, flat case). Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-m", "--meter",action="store_true",
                      help="Specify if the wakes are per meter length",
                      metavar="METER", default=False,dest="METER")
    parser.add_option("-c", "--constant",action="store_true",
                      help="Specify if we should plot the constant term Wycst",
                      metavar="CST", default=False,dest="CST")
    parser.add_option("-l", "--loglog",action="store_true",
                      help="Specify if loglog plot (semilogy by default)",
                      metavar="LOG", default=False,dest="LOG")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args

def logplot(z,W,lab,leg1,leg2,col):
    # plot wake and put legend and label
    pylab.loglog(z,np.abs(W),col+'-',label="|"+leg1+"|"+leg2,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Distance behind the source [m]")
    pylab.ylabel("W "+lab);

def semilogplot(z,W,lab,leg,col):
    # plot wake with labels, legend and color in semilogy
    pylab.semilogy(z,W,col+'-',label=leg,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Distance behind the source [m]")
    pylab.ylabel("W "+lab);

def plot(z,W,lab,leg,col):
    # plot wake with labels, legend and color
    pylab.plot(z,W,col+'-',label=leg,lw=2.5,ms=10.,mew=2.5)
    pylab.xlabel("Distance behind the source [m]")
    pylab.ylabel("W "+lab);


def read(filename):
    # read wake file (2 columns: z, wake)
    z=[];W=[];
    fh=open(filename,"r")
    for l in fh.readlines():
    	if not l.startswith('Dist'):
		ll=l.strip().split();
		z.append(ll[0]);
		W.append(ll[1]);
    fh.close()
    z=np.array([float(j) for j in z])
    W=np.array([float(j) for j in W])
    return z,W
    
                  
if __name__ == "__main__":
    opt,args=parsse();
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);

    # create list of associated legends (either with opt.LEG or 
    # the list of file names taking out all the common parameters
    # in the names)
    if (opt.LEG!=None):
        listleg=opt.LEG;
    else:
        listleg=takeout_common(listname);
    print listname,listleg


    if (opt.CST): nfig=6;
    else: nfig=5;
 
    for i in range(nfig): init_figure();
    
    color=['b','r','g','m','c','y','k'];

    for j,fil in enumerate(listname):

	# string for legend
	leg=listleg[j].replace(".dat","");
	#if opt.LEG==None: leg=fil.replace("WlongW","",1).replace("_"," ");
	#else: leg=opt.LEG[j];
	
	# constructs all the files names (flat chamber impedances)
	filefxdip=fil.replace("Wlong","Wxdip",1);
	filefxquad=fil.replace("Wlong","Wxquad",1);
	filefydip=fil.replace("Wlong","Wydip",1);
	filefyquad=fil.replace("Wlong","Wyquad",1);
	filefycst=fil.replace("Wlong","Wycst",1);

	# read longitudinal impedance and plot
	z,W=read(fil)
	if not opt.METER: lab=" [V / C]"
	else: lab=" [V / (C.m)]"
	pylab.figure(1);
	if opt.LOG: logplot(z,W,lab,"$W_{\|\|}$",", "+leg,color[j]);
	else: semilogplot(z,W,lab,"$W_{\|\|}$, "+leg,color[j]);
	#else: plot(z,W,lab,"$W_{\|\|}$, "+leg,color[j]);pylab.grid(pylab.gca());

	if opt.CST:
    	    # read constant vertical impedance and plot
    	    pylab.figure(6);
	    z,W=read(filefycst)
    	    if opt.LOG: logplot(z,W,lab,"$W_y^{cst}$",", "+leg,color[j]);
    	    else: semilogplot(z,W,lab,"$W_y^{cst}$, "+leg,color[j]);


	if not opt.METER: lab=" [V / (C.m)]"
	else: lab=" [V / (C.m$^2$)]"
	# read transverse dipolar impedances and plot them    
	z1,W1=read(filefxdip)
	z2,W2=read(filefydip)
	pylab.figure(2);
	if opt.LOG:
    	    logplot(z1,W1,lab,"$W_x^{dip}$",", "+leg,color[j]);
    	    pylab.figure(3);logplot(z2,W2,lab,"$W_y^{dip}$",", "+leg,color[j]);
	else:
    	    semilogplot(z1,W1,lab,"$W_x^{dip}$, "+leg,color[j]);
    	    pylab.figure(3);semilogplot(z2,W2,lab,"$W_y^{dip}$, "+leg,color[j]);

	# read transverse quadrupolar impedances and plot them
	z1,W1=read(filefxquad)
	z2,W2=read(filefyquad)
	pylab.figure(4);
	if opt.LOG:
    	    logplot(z1,W1,lab,"$W_x^{quad}$",", "+leg,color[j]);
    	    pylab.figure(5);logplot(z2,W2,lab,"$W_y^{quad}$",", "+leg,color[j]);
	else:
    	    semilogplot(z1,W1,lab,"$W_x^{quad}$, "+leg,color[j]);
    	    pylab.figure(5);semilogplot(z2,W2,lab,"$W_y^{quad}$, "+leg,color[j]);
    
    for i in range(5):
    	pylab.figure(i+1);
    	pylab.legend(loc=0);
	set_fontsize(pylab.figure(i+1),'xx-large');
	
    pylab.show();sys.exit()

