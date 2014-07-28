#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from parser_lib import *
import math
import matplotlib
import matplotlib.dates
from plot_lib import plot,init_figure,end_figure
from io_lib import read_ncol_file
from tables_lib import intersect,count_between

# plot two variables against each other, taken at the same 25ns-slot number

def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bar",action="callback",callback=multifloatint_parse,
                      help="Specify if, instead of a plot a y vs x, we want to plot bars (i.e. an histogram) of the number of y data points above a threshold (given by first argument) vs x, with a discretization in x given by the second argument. If a third argument is provided to this option, it is a threshold such that we count the number of x data points above this threshold in each range, and divide the initial bars by this. Default=no bars, simply plot y vs x.",
                      metavar="BAR", default=None,dest="BAR")
    parser.add_option("-l", "--log",type=int,
                      help="Specify if we want to plot in linear scale (0 - by default), in semilogx (1), in semilogy (2) or in loglog (3)",
                      metavar="LOG", default=0,dest="LOG")
    parser.add_option("-x", "--xaxis",action="callback",callback=multistring_parse,
                      help="Specify the file name for the x-axis. It should be a file with at least two columns and a single-line header. First column: 25-ns slot number or equivalent (see later), second column: variable for x-axis. If there 2 strings or more following -x, the second one is the label of the x-axis (by default it is simply the filename). If there are 3 strings, the last one indicates how to interpret the first column: 0-> directly as 25ns-slots numbers (by default), filename-> interpret it as the place of the bunch in the bunch train that is given by an FBCT-like file given here (this should be a two column - at least - file, with a single line header, as obtained by read_FBCT.py - average intensities vs non-empty 25ns-slot numbers).",
                      metavar="XAXIS", default=None,dest="XAXIS")
    parser.add_option("-y", "--yaxis",action="callback",callback=multistring_parse,
                      help="Specify the file name for the y-axis. It should be a file with at least two columns and a single-line header. First column: 25-ns slot number or equivalent (see later), second column: variable for y-axis. If there 2 strings or more following -y, the second one is the label of the x-axis (by default it is simply the filename). If there are 3 strings, the last one indicates how to interpret the first column: 0-> directly as 25ns-slots numbers (by default), filename-> interpret it as the place of the bunch in the bunch train that is given by an FBCT-like file given here (this should be a two column - at least - file, with a single line header, as obtained by read_FBCT.py - average intensities vs non-empty 25ns-slot numbers).",
                      metavar="YAXIS", default=None,dest="YAXIS")
    parser.add_option("-o", "--output",help="Specify output filename prefix. Default: ''correlation''",
                      default="correlation",dest="OUT")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")

    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args


def transform_bunches_nb_into_slots_from_scheme(bunches,scheme):
    # transform an array 'bunches' giving the position of bunches in the bunch train
    # (ordered from 0 to nbunch-1), into an array giving instead the slots numbers
    # of all these bunches
    # 'scheme' contains the filling scheme, i.e. the slot number of all bunches in the 
    # train.
    # Note that there can be less bunches in 'bunches' than in the filling scheme.
    
    slots=[];
    for bunch in bunches:
        slots.append(scheme[bunch]);
	
    return np.array(slots);


def transform_bunches_nb_into_slots_from_schemefilename(bunches,scheme_filename):
    # Wrapper of the previous function with name of the file where to find
    # the filling scheme.
    # This should be a two column (at least) file, with a single line header, 
    # as obtained by read_FBCT.py (average intensities vs non-empty 25ns-slot numbers).
    
    s=read_ncol_file(scheme_filename,ignored_rows=1);
    scheme=s[:,0];
    return transform_bunches_nb_into_slots_from_scheme(bunches,scheme);
    

if __name__ == "__main__":
    opt,args=parsse(); 

    # read x-axis data
    datax=read_ncol_file(opt.XAXIS[0],ignored_rows=1);
    
    if (len(opt.XAXIS)>1): xlab=opt.XAXIS[1];
    else: xlab=opt.XAXIS[0];

    if (len(opt.XAXIS)>2)and(opt.XAXIS[2]!=0):
    	slotsx=transform_bunches_nb_into_slots_from_schemefilename(datax[:,0],opt.XAXIS[2]); 
    else:
    	slotsx=datax[:,0];
    
    # read y-axis data
    datay=read_ncol_file(opt.YAXIS[0],ignored_rows=1);
    
    if (len(opt.YAXIS)>1): ylab=opt.YAXIS[1];
    else: ylab=opt.YAXIS[0];
   
    if (len(opt.YAXIS)>2)and(opt.YAXIS[2]!=0):
    	slotsy=transform_bunches_nb_into_slots_from_schemefilename(datay[:,0],opt.YAXIS[2]); 
    else:
    	slotsy=datay[:,0];
	
    # find the corresponding slots
    slots,indx,indy=intersect(slotsx,slotsy);
	    
    fig,ax=init_figure();
    if (opt.BAR==None):
	# plot y vs x
	plot(datax[indx,1],datay[indy,1],'','.',ylab,ax,opt.LOG,xlab=xlab);

    else:
        # plot bars
	# x discretization
	xbar=np.arange(np.min(datax[indx,1])-opt.BAR[1]/2.,np.max(datax[indx,1])+opt.BAR[1]/2.,opt.BAR[1]);
	ybar=np.zeros(len(xbar));
	for ix,x in enumerate(xbar[:-1]):
	    ybar[ix]=count_between(datax[indx,1],datay[indy,1],opt.BAR[0],x,xbar[ix+1]);
	
	if (len(opt.BAR)>2):
	    ynorm=np.zeros(len(xbar));
	    for ix,x in enumerate(xbar[:-1]):
		ynorm[ix]=count_between(datax[indx,1],datay[indy,1],opt.BAR[2],x,xbar[ix+1]);
		if ynorm[ix]!=0: ybar[ix]=ybar[ix]/float(ynorm[ix]);
	
	#outfile=open('toto.dat','w');
	#for ix,x in enumerate(xbar): print >> outfile, x,ybar[ix];
	#outfile.close();
	ax.bar(xbar+opt.BAR[1]/2.,ybar,facecolor='k',edgecolor='k',width=opt.BAR[1]*0.9)
	ax.set_xlabel(xlab);ax.set_ylabel(ylab);
	
    
    end_figure(fig,ax,save=opt.SAVE*(opt.OUT));
    
    if not(opt.SAVE): pylab.show();

    sys.exit()
