#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,random,pytz
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from parser_lib import *
import math
import matplotlib
from SussixNM import *
from multiturndata_riccardodemaria_modifNico import *
from plot_lib import init_figure,cmap,plot,end_figure,make_movie,build_colors,plot2D,plot_save_hist
from io_lib import list_files
from datetime_lib import plotdate,tt,tb
from read_Headtail_prt_fit import slideaver,envelop,fit
from read_Headtail_prt_sussix import extract
from read_Headtail_prt_coupledbunch import gettune


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunchtable",action="append",
                      help="Specify the bunch table (HEADTAIL like) filename for the filling scheme. There can be two -b options: one file for each beam.",
                      metavar="BUNCH", default=None,dest="BUNCH")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the ADT (.data) name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--makemovie",action="store_true",
                      help="Specify if we make a movie with the bunch positions histogram",
                      metavar="MAKE", default=False,dest="MAKE")
    parser.add_option("-n", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify 25ns-slot number (beginning at 0) or list of slots (i.e. 0:6:2 = [0 2 4 6]) for the raw data plot (several -n options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-o", "--output",help="Specify some suffix for the output filenames (empty string by default)",
                      default="",dest="OUT")
    parser.add_option("-q", "--tunes",type=float,nargs=2,
                      help="Specify the factional part of the tunes (x and y) (default=[0.28,0.31])",
                      metavar="TUNES", default=[0.28,0.31],dest="TUNES")
    parser.add_option("-s", "--timeoffset",type=float,
                      help="Specify time offset for conversion to gmt time, in seconds (default=3600 seconds). This might vary with summer/winter time...",
                      metavar="OFFSET", default=3600.,dest="OFFSET")
    parser.add_option("-t", "--threshold",type=float,
                      help="Specify threshold (in ADT units): if value in 2D plot is above this, the slot number is dumped in a file (default: we do not dump anything).",
                      metavar="THRES", default=None,dest="THRES")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    (opt, args) = parser.parse_args()
    #print "Selected files:", opt.FILE
    return opt, args


def readLHCDamper(filename):
    	
    # extract various data from the ADT file
    
    # read file
    if filename.startswith('LHCDamperBPM'):
    	# Riccardo's script kind of data
	adt=LHCDamperBPMData(filename);
    else:
    	print 'Problem in the data filename: not recognized';
	sys.exit();
	
    # number of bunches, of turns and bunch bucket numbers (25ns buckets) 
    # for the selected bunches
    nbunch=adt.bunches;
    nturns=adt.turns;
    turns=np.arange(float(nturns));

    return turns,nbunch,nturns,adt
	

def plot2Dvsbunchestime(time2D,data2D,dat,optsave=None,tit=''):

    # complete procedure to do a 2D plot vs bunches and time

    nbunch=len(data2D[0,:]);
    
    # sort time2D, and data2D accordingly
    ind=np.argsort(time2D)
    time2D.sort();
    data2Dbis=data2D[ind,:];

    timeinterval=np.ceil((time2D[-1]-time2D[0])*86400/10.);
    print "Time interval for plot: ", timeinterval

    # re-interpolate data on a regular mesh (for the y-axis i.e. time)
    timeintervalinterp=np.average(np.diff(time2D))/10;
    newtime2D=np.arange(time2D[0],time2D[-1],timeintervalinterp);
    newdata2D=np.zeros((len(newtime2D),nbunch),dtype=float);
    for bnum in range(nbunch):
	newdata2D[:,bnum]=np.interp(newtime2D,time2D,data2Dbis[:,bnum]);
    # 2D plot
    fig2D,ax2D=init_figure();
    plot2D(newdata2D,0,nbunch-1,newtime2D[0],newtime2D[-1],'25ns-slot number',
	'Local time on '+dat.strftime("%Y-%m-%d"),tit,ax2D);
    ax2D.yaxis_date(tz=None);
    ax2D.yaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
    ax2D.yaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))

    if (optsave!=None): end_figure(fig2D,ax2D,save=optsave);
    else: end_figure(fig2D,ax2D);
	

def alignADTdata(pos,nturns,nbunch,bunchtable):

    # re-align position data such that it fits with the filling scheme (provided as a bunch table,
    # like for HEADTAIL: each slot is at 1 or 0)
    # pos is a 2D nturns*nbunch array
    
    # read bunch table (note: slots between end of table and nbunch are supposed to be empty)
    scheme=np.zeros(nbunch);
    for l,line in enumerate(open(bunchtable)):
	#print line[0]
	scheme[l]=int(line[0]);

    # try different bunch offset until we find a correct matching
    minsum=nbunch+1;
    for boff in range(nbunch):
    	newpos=concatenate((pos[0,boff:],pos[0,:boff]));
    	#a=newpos*scheme;
	#b=a[scheme==1];
    	a=np.abs(newpos)+scheme;
	b=a[scheme==0];
	if (np.all(b==0)):
	    boffsol=boff; # offset solution
	    if (minsum==0): print "WARNING: data re-alignment has more than one possible solution !";
	    minsum=0;
	elif (np.sum(b!=0)<minsum):
	    boffsol=boff; # best solution so far
	    minsum=np.sum(b!=0);
    
    if (minsum>0): print "WARNING: data re-alignment: still ",str(minsum)," errors."
    possol=concatenate((pos[:,boffsol:],pos[:,:boffsol]),axis=1);
    #print scheme[:70];
    #print possol[0,:70];
    #print possol[0,:70]*scheme[:70]

    return possol;
    

if __name__ == "__main__":
    opt,args=parsse(); 

    gmt=pytz.timezone('Europe/Amsterdam');
    utc=pytz.timezone('UTC')
    
    # number of seconds to convert to gmt time (for plot...). This might vary.
    toff=opt.OFFSET;

    # tunes
    Qx=opt.TUNES[0]
    Qy=opt.TUNES[1]
    Qx=Qx-floor(Qx)
    Qy=Qy-floor(Qy)
    
    # revolution frequency and period (assumes protons)
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the circumference
    Trev=1./frev;
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    #print listname
    
    # read first file to get device names (assume they are always the same)
    turns,nbunch,nturns,adt=readLHCDamper(listname[0])
    
    # color pattern
    col=build_colors(nbunch);

    for j,dev in enumerate(adt.deviceNames):
    
	# initialization of 2D data
	data2D=np.zeros((len(listname),nbunch),dtype=float);
	tunem=np.zeros((len(listname),nbunch),dtype=float);
	ampm=np.zeros((len(listname),nbunch),dtype=float);
	time2D=np.zeros(len(listname),dtype=float);
	posaveturns=np.zeros(nbunch,dtype=float);
	
	if (opt.MAKE):
	    # initialization of histogram figure
	    fighist,axhist=init_figure();

	if (dev.find('Hor')!=-1):
	    Q=Qx;plane='x';
	else: 
	    Q=Qy;plane='y';

	print "Selected device: ", dev;
	
	if (opt.THRES!=None): 
	    fileunst=open('unstable_bunches_'+dev+opt.OUT+'.txt','w');
	    print >> fileunst, "slot (25ns)\ttime\taverage ADT abs. value"
	
	for i,filename in enumerate(listname):

    	    # string for the title
            #fil=filename.replace("_"," ").replace(".data","").replace("/"," ").replace("-"," ");
    	    print "Selected file:", filename;

    	    # read file
	    turns,nbunch,nturns,adt=readLHCDamper(filename);
	    
	    timestamp=getattr(adt,dev).acqStamp;
	
            if (i==0):
		# initial date and time
		t0=timestamp;
		dat=datetime.fromtimestamp(t0,tz=gmt).date();
		print dat,datetime.fromtimestamp(t0,tz=gmt).ctime();
		    

    	    # extract position data: 2D array of shape=(nturns,nbunch)
    	    pos=getattr(getattr(adt,dev),'pos');
	    
	    # re-align the slots with the filling scheme
	    if (len(opt.BUNCH)<=1): pos=alignADTdata(pos,nturns,nbunch,opt.BUNCH[0])
	    else:
	    	beam=int(dev[-1]); # beam number
		pos=alignADTdata(pos,nturns,nbunch,opt.BUNCH[beam-1]);
	    
	    pos=pos-np.average(pos,0);
	    
	    # bunch-by-bunch Sussix (DOES NOT WORK)
	    #for bnum in range(nbunch):
	    #    x=pos[:,bnum];
	    #	if (np.max(np.abs(x))>0):
	#	    suss=gettune(x,np.zeros(len(x)),Q,0.,0.07)
	#	    # sort and extract only tunes at 0.5 from nominal tune Q
	#	    tunes,amps=extract(suss.ox,suss.ax,Q,0.5)
	#	    #print suss.ox,suss.ax
	#	    # find maximum peak
	#	    tunem[i,bnum]=tunes[np.argmax(amps)]
	#	    ampm[i,bnum]=np.max(amps)
	#	else:
	#	    tunem[i,bnum]=0.;
	#	    ampm[i,bnum]=0.;


	    if opt.MAKE:
		for iturn in range(nturns): plot_save_hist(np.arange(nbunch),pos[iturn,:],datetime.fromtimestamp(timestamp,tz=gmt).ctime()+', turn '+str(iturn),axhist,fighist,iturn+i*nturns,xlim=None,ylim=None,tit=dev);

	    # computes the average value vs turns for all bunches
	    posaveturns=np.average(np.abs(pos),0);
	    # Tried max instead of average -> worst
	    # Tried to take log -> much less clear
	    
	    # for the 2D plot
	    #time2D[i]=np.average((timestamp+toff+turns*Trev)/86400.);
	    time2D[i]=np.average((timestamp+toff+turns*Trev)/86400.);
	    data2D[i,:]=posaveturns;
	    
	    # extract those above threshold
	    if (opt.THRES!=None):
	    	ind=np.where(posaveturns>=opt.THRES);
		if (len(ind)>0):
		    for bun in ind[0]:
		    	tim=datetime.fromtimestamp(time2D[i]*86400-toff,tz=gmt)
		    	print >> fileunst, bun, "\t", tim.strftime("%H%M%S"),"\t", posaveturns[bun];
	    
    
	if (opt.BNUM!=None):
	    for bnum in opt.BNUM:
		figraw,axraw=init_figure();
		plot(time2D,data2D[:,bnum],'Slot '+str(bnum),'-','ADT pickups data (a.u.)',axraw,0,
	    	    xlab="Local time on "+dat.strftime("%Y-%m-%d"));
    		timeinterval=np.ceil((time2D[-1]-time2D[0])*86400/10.);
		axraw.set_title(dev);
    		axraw.xaxis_date(tz=None);
    		axraw.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
    		axraw.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
		if (opt.SAVE): end_figure(figraw,axraw,save=listname[0]+"_raw_"+dev+opt.OUT+'slot'+str(bnum));
		else: end_figure(figraw,axraw);
	
	    
	if opt.MAKE:
	    # make movie
	    make_movie(listname[0]+"_"+dev+opt.OUT+'.gif','_tmp',flagrm=True);	
	
    	
	if (len(listname)>1):
	    # 2D plot of average amplitudes, slot-by-slot
	    if not(opt.SAVE): plot2Dvsbunchestime(time2D,data2D,dat,optsave=None,tit='Average raw amplitude (abs), '+dev);
	    else: plot2Dvsbunchestime(time2D,data2D,dat,optsave=listname[0]+"_raw_"+dev+opt.OUT,tit='Average raw amplitude (abs), '+dev);

	    # 2D plot of maximum SUSSIX amplitudes, slot-by-slot
	    #if not(opt.SAVE): plot2Dvsbunchestime(time2D,ampm,dat,optsave=None,tit='Maximum Sussix line amplitude, '+dev);
	    #else: plot2Dvsbunchestime(time2D,ampm,dat,optsave=listname[0]+"_sussix_"+dev+opt.OUT,tit='Maximum Sussix line amplitude, '+dev);
	
	

    if not(opt.SAVE): pylab.show();

    sys.exit()

