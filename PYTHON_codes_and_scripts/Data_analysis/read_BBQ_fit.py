#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
from collimator_settings import concatenate_data_Timber
import pylab,re,dateutil,random,pytz
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from parser import *
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
from Timber import parseout
import math
import matplotlib
import matplotlib.dates
from plot_lib import init_figure,end_figure,cmap,plot
from io_lib import list_files
from read_Headtail_prt_fit import fit
from collimator_settings import concatenate_data_Timber
from datetime_lib import tt
#import subprocess
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--analysis",action="store_true",
                      help="Specify if a data analysis has to be performed",
                      metavar="ANA", default=False,dest="ANA")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify interval (on the measured signal) on which we do the fit (default= 6e7 & 1.5e9)",
                      metavar="BOUND", default=[6.e7,1.5e9],dest="BOUND")		      
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of seconds from the beginning after which we want to begin the fit",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-l", "--plot",type=int,
                      help="Specify 1 to plot also maxima, 2 to plot maxima and sliding average, 3 to plot sliding average (it also always plot raw data and fit)",
                      metavar="PLOT", default=0,dest="PLOT")
    parser.add_option("-m", "--multi",action="store_true",
                      help="Specify if there are several files saved with the multifile option in Timber",
                      metavar="MULT", default=False,dest="MULT")
    parser.add_option("-n", "--nbfile",type=int,
                      help="Specify maximum number of BBQ data files to take",
                      metavar="NFILE", default=10000,dest="NFILE")
    parser.add_option("-o", "--output",help="Specify some suffix for the output filenames (empty string by default)",
                      default="",dest="OUT")
    parser.add_option("-p", "--plane",
                      help="Specify plane: h or v",
                      metavar="PLANE", default=None,dest="PLANE")
    parser.add_option("-r", "--period",type=int,nargs=2,
                      help="Specify periods used for 1) extracting the envelop and 2) the sliding average (default=60 & 100)",
                      metavar="PERIOD", default=[60,100],dest="PERIOD")		      
    parser.add_option("-w", "--plotraw",type=int,
                      help="Specify some option for the plot of the raw data: 0 (default)= plot all turns vs. time, 1=plot all turns vs. number of turns, n=plot every n turns vs. time",
                      metavar="PLOTRAW", default=0,dest="PLOTRAW")
    parser.add_option("-x", "--extradata",action="callback",callback=multistring_parse,
                      help="Specify another TIMBER datafile and variable name (actually, only part of the full name, enough to identify it) to plot on the tune graph. File name can be a regular expression between quotes "" (for multifiles). Optionally, two additional arguments may be given (in this order): legend to put on the graph, and label for the upper x-axis.",
                      metavar="EXTRA", default=None,dest="EXTRA")
    parser.add_option( "--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    # Note: option -l, -d, -r and -g used only when -a option activated	      
    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args


#def plot(t,data,leg,dat,col,ylab,ax,lw=1.):
#    gmt=pytz.timezone('Europe/Amsterdam');
#    ax.plot_date(t,data,col+'-',label=leg,tz=gmt,lw=lw)
#    ax.set_xlabel("Time on "+dat.strftime("%Y-%m-%d"))
#    ax.set_ylabel(ylab);

def plotturns(BBQ,leg,dattim,col,ylab,ax):
    ax.plot(range(0,len(BBQ)),BBQ,col+'-',label=leg)
    ax.set_xlabel("Number of turns from "+dattim.strftime("%Y-%m-%d %H:%M:%S"))
    ax.set_ylabel(ylab);


def read_BBQ_multi(listname,energy,beam,plane,timeinterval=0.08,nturn=1024):

    # read BBQ data in several files, and put all the data together in single arrays
    # NOTE: - timeinterval is the time interval between lines in the data;
    # it has changed from 0.4s to 0.08s between 2011 and 2012,
    #	    - nturn is the number of turns at each line; it has changed from 
    # 8192 to 1024 between 2011 and 2012.
    
    # LHC revolution frequency and period
    print 'Energy: ',energy,' GeV';
    gamma=energy/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    LHCfrev=beta*299792458./26658.883; # the later is the LHC circumference
    LHCTrev=1./LHCfrev;
    # number of turns assumed to be present at each line
    #nturn=8192;
    # parameter for precise concatenation
    param=2;
    
    # concatenate datafiles in one single data set
    dataBBQ=[];timeBBQ=[];
    for j,name in enumerate(listname):
        sys.stdout.write('%s: '% name);
	if name.endswith('tsv'): sepa='\t';
	else: sepa=',';
        data=parseout(name,sep=sepa);
	for v in data['datavars']:
	    if (v.startswith('LHC.BQBBQ.')) and (v.endswith(beam+':ACQ_DATA_'+plane)): var=v;
	# note: in timber data column 0 is always time and column 1 is always the list 
	# of all the other columns
	#print len(data[var][1][0]) # here 1 is column number, 0 is line number
	
	ind1=-1;
	if (j>0):
	    cnt=data[var][1][0].count(dataBBQ[-1]);
	    if (cnt>0):
		ind1=data[var][1][0].index(dataBBQ[-1]);
		flag=(data[var][1][0][ind1-1]==dataBBQ[-2]);
		# ind1 should be around nturn-(timeinterval*LHCfrev)
		while (not flag)and(abs(ind1-nturn+(timeinterval*LHCfrev))>param):
	    	    ind1=data[var][1][0][ind1+1:].index(dataBBQ[-1])+ind1+1;
		    flag=(data[var][1][0][ind1-1]==dataBBQ[-2]);
		    print ind1,flag
		
		#print ind1,dataBBQ[-1],data[var][1][0][ind1],dataBBQ[-2],data[var][1][0][ind1-1];
		#sys.stdout.write('%lf\n' %timeBBQ[-1]);
		
	    else: print "Missing line in ",name," around ",data[var][0][0]
	
	dataBBQ.extend(data[var][1][0][ind1+1:]);
	timebeg=pylab.date2num(tt(data[var][0][0][:-4]))*86400+float(data[var][0][0][-4:]);
	timeBBQ.extend([timebeg+k*LHCTrev for k in range(ind1+1,nturn)]);
        #sys.stdout.write('%lf %lf %d\n'% (timebeg+(ind1+1)*LHCTrev,timeBBQ[-1],len(timeBBQ)));

	
	# connect lines between them
	for k,tps in enumerate(data[var][0][1:]):
	    #if (k<=10):print k,dataBBQ[-1]
	    #print k,dataBBQ[-1]
	    cnt=data[var][1][k+1].count(dataBBQ[-1]);
	    if (cnt>0):
		ind1=data[var][1][k+1].index(dataBBQ[-1]);
		flag=(data[var][1][k+1][ind1-1]==dataBBQ[-2]);
		#print ind1, flag
		# ind1 should be around nturn-(timeinterval*LHCfrev)
		kk=1;
		while ((not flag)and(abs(ind1-nturn+(timeinterval*LHCfrev))>param))and(kk<=cnt-1):
	    	    ind1=data[var][1][k+1][ind1+1:].index(dataBBQ[-1])+ind1+1;
		    flag=(data[var][1][k+1][ind1-1]==dataBBQ[-2]);
		    kk=kk+1;
		    #print ind1,flag
		if (not flag)or(abs(ind1-nturn+(timeinterval*LHCfrev))>param): ind1=-1; print "Missing line in ",name," around ",data[var][0][k+1],", ind1=",ind1,", data at ind1=",data[var][1][k+1][ind1];
		    
	    else: ind1=-1;print "Missing line in ",name," around ",data[var][0][k+1]

	    #print ind1,dataBBQ[-1],data[var][1][k+1][ind1],dataBBQ[-2],data[var][1][k+1][ind1-1];
	    dataBBQ.extend(data[var][1][k+1][ind1+1:]);
	    timebeg=pylab.date2num(tt(data[var][0][k+1][:-4]))*86400+float(data[var][0][k+1][-4:]);
	    timeBBQ.extend([timebeg+kk*LHCTrev for kk in range(ind1+1,nturn)]);
            #sys.stdout.write('%lf %lf %d\n'% (timebeg+(ind1+1)*LHCTrev,timeBBQ[-1],len(range(ind1+1,nturn))));


    # convert data to arrays
    tp=np.array(timeBBQ);
    BBQ=np.array([float(k) for k in dataBBQ])
    
    return tp,BBQ,LHCfrev


def plotextradata(optextra,ax,t0,tend,flagtime=0,patcol='-r',datet=None):

    # plot extra data set on axis ax
    # filename where to find data in optextra[0], part of variable name in optextra[1]
    # and (if present) legend in optextra[2] and x-axis label in optextra[3]
    # initial time in t0, final time in tend
    # flagtime is 0 if one wants to plot vs time, frev if it is vs turn numbers
    # patcol is the pattern+color for the plot
    
    # read file
    # create list of filenames
    listnameEXT=list_files([optextra[0]]);

    # legend
    if (len(optextra)>=3): leg=optextra[2];
    else: leg=optextra[1];
    # second x or y label
    if (len(optextra)>=4): label=optextra[3];
    else: label=optextra[1];

    # concatenate datafiles in one single data set
    data=concatenate_data_Timber(listnameEXT);
    for v in data['datavars']:
    	if (v.find(optextra[1])!=-1):
	    var=v;
    	    tpsother=np.array([pylab.date2num(tt(j[:-4])) for j in data[v][0]])*86400.+np.array([float(j[-4:]) for j in data[v][0]]);
	    dataother=np.array([float(j[0]) for j in data[v][1]]);
    # find times corresponding to the BBQ data
    #print var;
    #beg=np.where(tpsother>=tps[0]+opt.BEG/frev);beg=beg[0];
    beg=np.where(tpsother>=t0);beg=beg[0];
    if (len(beg)>0): beg=beg[0]-(beg[0]!=0); # beg[0]-1 if beg[0] is not 0, otherwise 0
    else: beg=0;
    end=np.where(tpsother>=tend);end=end[0];
    if (len(end)>0)and(end[0]<len(tpsother)-1): end=end[0]+1;
    else: end=-1;
    tpsnew=tpsother[beg:end];datanew=dataother[beg:end];
    print dataother,beg,end
    # plot
    ylim=ax.get_ylim();
    ax2=pylab.twinx(ax=ax); # second y axis
    if (flagtime==0):
	if (datet==None): labx='Time'
	else: labx="Time on "+datet.strftime("%Y-%m-%d")
    	plot(tpsnew/86400.,datanew,leg,patcol,label,ax2,0,xlab=labx);
    else:
    	frev=flagtime;
	plot((tpsnew-t0)*frev,datanew,leg,patcol,label,ax2,0,xlab="Number of turns");
    ax2.legend(loc=2);
    ax.set_ylim(ylim);
    
    return tpsnew,datanew,var,leg,label,ax2;
    


if __name__ == "__main__":
    opt,args=parsse(); 

    if opt.BEAM=="1":
    	beam="B1";
    elif opt.BEAM=="2":
    	beam="B2";
    else: print "specify beam 1 or 2"; sys.exit()

    if (opt.PLANE=="h") or (opt.PLANE=="H") or (opt.PLANE=="x"):
    	plane="H";
    elif (opt.PLANE=="v") or (opt.PLANE=="V") or (opt.PLANE=="y"):
    	plane="V";
    else: print "specify plane h or v"; sys.exit()
    

    # construct list of filenames
    root=opt.FILE[:-21]+'*.csv'
    #print root;
    if (opt.MULT):
    	#subprocess.Popen(['ls',root]); # doesn't work
	listname=glob.glob(root);
    else:
	listname=[];
	listname.append(opt.FILE);
	
    listname.sort();
    ind=listname.index(opt.FILE);
    listname=listname[ind:min(opt.NFILE+ind,len(listname))]
    for j in listname: print j;
	
    # read all the BBQ data in listname
    tp,BBQ,frev=read_BBQ_multi(listname,opt.EN,beam,plane);
    
    # initialize plot
    fig,ax=init_figure(); 

    # convert time to sth compatible with pylab.plot_date -> USELESS
    #tplot=np.array([pylab.num2date(j/86400.) for j in tp]);
    #tplot=pylab.num2date(tp/86400.);
    # find date
    dat=pylab.num2date(tp[0]/86400.).date();
    dattim=pylab.num2date(tp[0]/86400.);
    
    # plot raw data
    #plot(tplot,BBQ,'Measurement '+beam+' '+plane,dat,'b',"BBQ (arbitrary unit)",ax);
    #if (opt.PLOTRAW==0): plot(tp/86400.,BBQ,'BBQ '+beam+' '+plane,dat,'b',"BBQ (arbitrary unit)",ax);
    #elif (opt.PLOTRAW==1): plotturns(BBQ,'BBQ '+beam+' '+plane,dattim,'b',"BBQ (arbitrary unit)",ax);
    #else: plot(tp[::opt.PLOTRAW]/86400.,BBQ[::opt.PLOTRAW],'BBQ '+beam+' '+plane,dat,'b',"BBQ (arbitrary unit)",ax);
    if (opt.PLOTRAW==0): plot(tp/86400,BBQ,'BBQ '+beam+plane,'b',"BBQ (arbitrary unit)",ax,0,xlab="Time on "+dat.strftime("%Y-%m-%d"),plotevery=1);
    elif (opt.PLOTRAW==1): plotturns(BBQ,'BBQ '+beam+' '+plane,dattim,'b',"BBQ (arbitrary unit)",ax);
    else: plot(tp/86400,BBQ,'BBQ '+beam+plane,'b',"BBQ (arbitrary unit)",ax,0,xlab="Time on "+dat.strftime("%Y-%m-%d"),plotevery=opt.PLOTRAW);
    #file=open('out','w');
    #for i,j in enumerate(BBQ): print >> file, tplot[i],j;
    #file.close();
	    
    
    if (opt.ANA):
	# data analysis
	# take out averagek,j in enumerate(BBQ)
	aver=np.average(BBQ);
	print "BBQ average: ",aver;
	BBQ=BBQ-aver;
	
	if (opt.PLOTRAW!=1):
	    tau=fit(tp,BBQ,ax,True,"BBQ "+beam+plane,"BBQ (arbitrary unit)",1./frev,opt.BOUND,
		opt.PLOT,opt.PERIOD,opt.BEG,0,col='g',datet=dat);
	else:
	    tau=fit((tp-tp[0])*frev,BBQ,ax,True,"BBQ "+beam+plane,"BBQ (arbitrary unit)",1./frev,opt.BOUND,
		opt.PLOT,opt.PERIOD,opt.BEG,0,col='g');
	
	filetau=open(listname[0]+'_raw_tau_'+beam+plane+opt.OUT+'.txt','w');
    	print >> filetau, "\ttau[s]"
	print "Rise time in sec.: ", tau;
	print >> filetau, "\t", tau;
	filetau.close();

	
    if (opt.EXTRA!=None):
	# plot another data file
	tpsnew,datanew,var,leg,label,ax2=plotextradata(opt.EXTRA,ax,tp[0],tp[-1],flagtime=int((opt.PLOTRAW==1)*frev),datet=dat);

    if (opt.PLOTRAW!=1):
	ax.set_xlabel="Time on "+dat.strftime("%Y-%m-%d")
	timeinterval=np.ceil((tp[-1]-tp[0])/5.);
	print "Time interval for plot: ", timeinterval
	ax.xaxis_date(tz=None);
	if (timeinterval<60):
    	    ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
	elif (timeinterval<3600):
    	    ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=int(np.floor(timeinterval/60.))))
	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
	else:
    	    ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=int(np.floor(timeinterval/3600.))))
	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))


    ax.legend(loc=1);
    if opt.EXTRA!=None: end_figure(fig,ax2); 
    if opt.SAVE:
	end_figure(fig,ax,legpos=3*(opt.EXTRA!=None),save=listname[0]+'_raw_'+beam+plane+opt.OUT);
    else:
	end_figure(fig,ax,legpos=3*(opt.EXTRA!=None));
	pylab.show();

    sys.exit()

