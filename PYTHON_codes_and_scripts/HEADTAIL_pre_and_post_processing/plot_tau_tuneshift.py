#!/usr/bin/python2.6

import sys
#sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
#import pylab,re,random
import numpy as np
from string import split, replace
from parser import *
import math
import pylab
from diff_tau_tuneshift_Headtail import read
from plot_lib import set_fontsize,init_figure

def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) to compare on (several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--bounds",type=float,nargs=2,
                      help="Specify the bounds for the y axis of the rise time plot (default=autoscale)",
                      metavar="BOUND", default=None,dest="BOUND")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the files names (several files possible)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--legend",action="append",
                      help="Specify the legend for each file to be plotted",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-o", "--offset",type=float,nargs=2,action="append",
                      help="Specify the offset (in number of bunches), and multiplication factor, to apply to the first column of each file (to deal with artificial variation of bunch numbers)",
                      metavar="OFF", default=None,dest="OFF")
    parser.add_option("-p", "--plane",
                      help="Specify the plane 'x' or 'y' (default= both x and y)",
                      metavar="PLANE", default=None,dest="PLANE")
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    #print "Selected bunches:", opt.BNUM
    return opt, args


if __name__ == "__main__":

    print "\n"
    
    opt,args=parsse(); 

    print 'Energy: ',opt.EN,' GeV, circumference: ',opt.CIRC;
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the circumference
    Trev=1./frev;
   
    
    if (opt.BNUM!=None):
	bunches=np.array(opt.BNUM);
	bunches.sort();
	print "Selected bunches:", bunches


    if (opt.OFF==None): opt.OFF=np.zeros((len(opt.FILE),2))
    else:
        if (len(opt.OFF)<len(opt.FILE)): 
	    print "Not enough number of bunch offsets defined !"
	    sys.exit()
	
	
    if (opt.LEG!=None):
        if (len(opt.LEG)<len(opt.FILE)):
	    print "Not enough number of legends defined !"
	    sys.exit()


    figtau,axtau=init_figure()
    figrealtune,axrealtune=init_figure()
    figcomptune,axcomptune=init_figure()
	
    col=['b','r','g','m','k','c','y']
   
    # main loop: plot for each file
    for i,filename in enumerate(opt.FILE):
    
    	print filename;
	fil=filename.replace("_"," ").replace(".txt","").replace("/"," ").replace("-"," ")
	
    	# string for the legend
	if (opt.LEG==None):
            leg=fil
	else:
	    leg=opt.LEG[i]
    	
	# read files
	data=read(filename)

	if (opt.BNUM==None)and(i==0): bunches=(np.arange(len(data[:,0]))+opt.OFF[0][0])*opt.OFF[0][1]
	
	# select the right lines
	ind=[];
	for bnum in bunches:
	    a=np.int_(data[:,0]*opt.OFF[i][1]+opt.OFF[i][0])
	    #print a,bnum
	    j=np.where(a==bnum);#print j
	    ind.append(j[0][0])
	
	bunchi=data[ind,0]*opt.OFF[i][1]+opt.OFF[i][0]
	print bunchi
	
	if (data.shape[1]>=5)and(opt.PLANE==None):
	    # plot rise times in seconds
	    axtau.plot(bunchi,data[ind,1],col[i%7]+'-x',label='Horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    axtau.plot(bunchi,data[ind,2],col[i%7]+'--x',label='Vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	    # plot real and complex tune shift (if present in the data)
	    axrealtune.plot(bunchi,data[ind,3],col[i%7]+'-x',label='Horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    axrealtune.plot(bunchi,data[ind,4],col[i%7]+'--x',label='Vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,data[ind,3],col[i%7]+'.-x',label='Real, horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,data[ind,4],col[i%7]+':x',label='Real, vertical '+leg,lw=2.5,ms=10.,mew=2.5);
    	    axcomptune.plot(bunchi,Trev/(data[ind,1]*2.*np.pi),col[i%7]+'-x',label='Imag. , horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,Trev/(data[ind,2]*2.*np.pi),col[i%7]+'--x',label='Imag. , vertical '+leg,lw=2.5,ms=10.,mew=2.5);

	elif (data.shape[1]==3):
	    # if only 3 columns -> ADT data, plot only rise times in seconds
	    if (fil.find('Hor')!=-1):
	    	axtau.plot(bunchi,data[ind,1],col[i%7]+'-x',label='Horizontal '+leg,lw=2.5,ms=10.,mew=2.5);    
	    else:
	        axtau.plot(bunchi,data[ind,1],col[i%7]+'-x',label='Vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	
	elif (data.shape[1]>=5)and(opt.PLANE=='x'):   
	    # plot rise times in seconds
	    axtau.plot(bunchi,data[ind,1],col[i%7]+'-x',label='Horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    # plot real and complex tune shift (if present in the data)
	    axrealtune.plot(bunchi,data[ind,3],col[i%7]+'-x',label='Horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,data[ind,3],col[i%7]+'-x',label='Real, horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
    	    axcomptune.plot(bunchi,Trev/(data[ind,1]*2.*np.pi),col[i%7]+'--x',label='Imag. , horizontal '+leg,lw=2.5,ms=10.,mew=2.5);
	
	elif (data.shape[1]>=5)and(opt.PLANE=='y'):
	    # plot rise times in seconds
	    axtau.plot(bunchi,data[ind,2],col[i%7]+'-x',label='Vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	    # plot real and complex tune shift (if present in the data)
	    axrealtune.plot(bunchi,data[ind,4],col[i%7]+'-x',label='Vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,data[ind,4],col[i%7]+'-x',label='Real, vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	    axcomptune.plot(bunchi,Trev/(data[ind,2]*2.*np.pi),col[i%7]+'--x',label='Imag. , vertical '+leg,lw=2.5,ms=10.,mew=2.5);
	
	
	else: print "Not the right shape for the data!", filename;sys.exit()  


    axtau.set_xlabel("Bunch number");
    axtau.set_ylabel("Rise time [s]");
    if (opt.BOUND!=None): axtau.set_ylim(opt.BOUND)
    axtau.legend(loc=0);
    set_fontsize(figtau,'xx-large');

    axrealtune.set_xlabel("Bunch number");
    axrealtune.set_ylabel("Real tune shift");
    axrealtune.legend(loc=0);
    set_fontsize(figrealtune,'xx-large');

    axcomptune.set_xlabel("Bunch number");
    axcomptune.set_ylabel("Complex tune shift");
    axcomptune.legend(loc=0);
    set_fontsize(figcomptune,'xx-large');
    
    pylab.show()

