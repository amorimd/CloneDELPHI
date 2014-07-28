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
from multiturndata_riccardodemaria_modifNico import *
from plot_lib import init_figure,cmap,plot,end_figure,build_colors,plot2D,set_fontsize
from io_lib import list_files
from read_Headtail_prt_fit import slideaver,envelop,fit
from datetime_lib import tt,tb,plotdate


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--analysis",action="store_true",
                      help="Specify if a data analysis has to be performed",
                      metavar="ANA", default=False,dest="ANA")
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 0) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the raw data plot (several -b options possible). Put a space between -b and list or number.",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify boundaries for the fit of individual bunches (default=[70,30000])",
                      metavar="BOUND", default=[70,30000],dest="BOUND")		      
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the ADT (.data) name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the fit",
                      metavar="BEG", default=2000,dest="BEG")
    parser.add_option("-i", "--ratio",type=float,
                      help="Specify ratio above which sudden single-point peaks should be suppressed (when using pre-treatment) (default=1.3)",
                      metavar="RATIO", default=1.3,dest="RATIO")
    parser.add_option("-l", "--plot",type=int,
                      help="Specify 1 to plot also maxima, 2 to plot maxima and sliding average, 3 to plot sliding average (it also always plots raw data and fit)",
                      metavar="PLOT", default=0,dest="PLOT")
    parser.add_option("-m", "--boundbeam",type=float,nargs=2,
                      help="Specify boundaries for the fit of the whole beam (default=[50,3000])",
                      metavar="BOUNDBM", default=[50,3000],dest="BOUNDBM")		      
    parser.add_option("-n", "--oneplot",action="store_true",
                      help="Specify if we plot all files in one single plot (vs. time). Then we also do a 2D colour plot with average absolute amplitude of each bunch vs time.",
                      metavar="ONEPLOT", default=False,dest="ONEPLOT")
    parser.add_option("-o", "--output",help="Specify the end of the output filename for rise times (default=tau.txt)",
                      default="tau.txt",dest="OUT")
    parser.add_option("-p", "--plotlog",type=int,
                      help="Specify if we want to plot in linear scale (0), in semilogx (1), in semilogy (2) or in loglog (3) (default=linear scale)",
                      metavar="LOG", default=0,dest="LOG")
    parser.add_option("-r", "--period",type=int,nargs=2,
                      help="Specify periods used for 1) extracting the envelop and 2) the sliding average (default=[60,10])",
                      metavar="PERIOD", default=[60,10],dest="PERIOD")		      
    parser.add_option("-s", "--slideaver",type=int,
                      help="Specify period used for the sliding average subtracted BEFORE analysis, when using pre-treatment (to get rid of constant offset + some low frequency oscillations). If 0, no sliding average done (default=20).",
                      metavar="SLIDE", default=20,dest="SLIDE")		      
    parser.add_option("-t", "--nopretreatment",action="store_false",
                      help="Specify if a pre-treatment of the raw data SHOULD NOT be performed",
                      metavar="PRE", default=True,dest="PRE")
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot (and fit with -a option) the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-z", "--rolloffset",type=int,
                      help="Specify the additional offset used in the data rolling (such that instability at the end), when using pre-treatment. -1 means no roll (default=10).",
                      metavar="ROLL", default=10,dest="ROLL")		      
    # Note: options -c, -e, -o, -l, -d, -r and -g used only when -a option activated
    # -i and -s options used only if -t option activated	      
    (opt, args) = parser.parse_args()
    #print "Selected files:", opt.FILE
    return opt, args


def takeout_sudden_peak(turns,pos,ratio,nbunch,deviceNames):

    # take out any sudden (single-point) peak of array 'pos'

    # first take out sliding average
    periodaver=20; # number of points on which we average
    posaver=[];
    pos2=[];
    for j,x in enumerate(pos):
	posaver.append(np.transpose(np.array([np.average(x[:,k:min(k+periodaver+1,len(x[0]))],1) for k,x0 in enumerate(x[0,:])])));
	pos2.append(pos[j]-posaver[j]);

    # extraction of maxima, then find points with big difference w.r.t previous and next ones
    for j,x in enumerate(pos2):
	takeout=[];takeoutcnt=[];
	for bnum in range(nbunch):
	    # extract envelop (maxima)
	    turnsenv,env=envelop(turns,np.abs(x[bnum]),10);
	    # find the sudden peaks
	    ind=np.where(np.abs(np.diff(env)/env[:-1])>ratio);
	    #ind2=np.where(np.abs(np.diff(env[1:])/env[:-2])>4.);
	    if (len(ind[0])>=1):
		for i,ind0 in enumerate(ind[0]):
		    if (ind0<len(env)-2):
			if (sign(env[ind0+2]-env[ind0+1])==-sign(env[ind0+1]-env[ind0]))and(abs((env[ind0+2]-env[ind0+1])/env[ind0+2])>ratio):
			    if (turnsenv[ind0+1] not in takeout): 
				takeout.append(turnsenv[ind0+1]);
				takeoutcnt.append(1);
			    else:
				k=takeout.index(turnsenv[ind0+1]);
				takeoutcnt[k]=takeoutcnt[k]+1;
	# take away only such turns that are identical in different bunches
	notakeout=[];
	for k,i in enumerate(takeout):
	    if (takeoutcnt[k]<=1): notakeout.append(i);
	for k in notakeout: takeout.remove(k);
	print deviceNames[j],', points taken out: ',takeout
	#if (len(takeout)>1): print "Warning: pre-treatment: ",len(takeout)," points taken out for ",deviceNames[j]
	takeout.sort();
	# take out those turn numbers and replace by average of two nearest values
	y=pos[j];
	for turn in takeout: 
	    if (turn>0)and(turn<len(turns)-1): y[:,turn]=(y[:,turn-1]+y[:,turn+1])/2.;
	    elif (turn==0): y[:,turn]=y[:,turn+1];
	    else: y[:,turn]=y[:,turn-1];
	pos[j]=y;
	
	
    return pos

                  
def substract_slideaverage(turns,pos,periodaver,nbunch):

    # sliding average, then substraction of the sliding average
    # periodaver is the number of points on which we average
    posaver=[];
    turnsaver=slideaver(turns,periodaver);
    for j,x in enumerate(pos):
	posaver.append(np.transpose(np.array([np.average(x[:,k:min(k+periodaver+1,len(x[0]))],1) for k,x0 in enumerate(x[0,:])])));
	y=[];
	for bnum in range(nbunch):
	    y.append(np.interp(turnsaver,turns,pos[j][bnum])-posaver[j][bnum]);
	pos[j]=np.array(y);
    turns=turnsaver;
    
    return turns,pos


def roll_data(turns,pos,nbunch,offset=10):

    # roll the ADT data such that the instability is always at the end
    for j,x in enumerate(pos):
	y=[];
	for bnum in range(nbunch):
	    maxi=int(turns[np.argmax(np.abs(x[bnum]))]);
	    y.append(np.roll(x[bnum,:],len(turns)-maxi-offset));
	pos[j]=np.array(y);
	
    return pos
    

def readADT(filename):
    	
    # extract various data from the ADT file
    
    # read file
    if filename.startswith('ADT'):
    	# Multiturn application (from Verena) kind of data
    	adt=MultiTurnData(filename);
	timestamp=adt.timestamp;
    else:
    	print 'Problem in the data filename: not recognized';
	sys.exit();
	
    # beam, number of bunches, of turns and bunch bucket numbers (25ns buckets) 
    # for the selected bunches
    #beam=adt.beam; # for Verena kind of data only; useless anyway
    nbunch=adt.bunches;
    nturns=adt.turns;

    # extract position data
    pos=[];
    for dev in adt.deviceNames: pos.append(getattr(getattr(adt,dev),'pos'));
    turns=np.arange(float(nturns));

    return turns,pos,nbunch,nturns,adt,timestamp
	

def convert_time(ti):

    # useless function
    
    # convert time from ADT file timestamp, to correct the year (add 1970=2012-43+1 years)
    #print date(1970,1,1);
    return ti+86400*pylab.date2num(date(1970,1,1));
	

if __name__ == "__main__":
    opt,args=parsse(); 

    gmt=pytz.timezone('Europe/Amsterdam');
    
    # number of seconds to convert to gmt time (for plot...)
    toff=7200;

    # revolution frequency and period (assumes protons)
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the circumference
    Trev=1./frev;
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    
    for i,filename in enumerate(listname):
    
    	# string for the title
        fil=filename.replace("_"," ").replace(".data","").replace("/"," ").replace("-"," ");
    	print "Selected file:", filename
	
    	# read file
	turns,pos,nbunch,nturns,adt,timestamp=readADT(filename);
	
        if (i==0):
	    # initial date and time
	    dat=datetime.fromtimestamp(timestamp,tz=gmt).date();
	    t0=timestamp;
	    print dat,datetime.fromtimestamp(t0,tz=gmt).ctime();
	    # color pattern
	    col=build_colors(nbunch);
	    if (opt.ONEPLOT):
    		# data for 2D plot: initialization
		data2D=[];
    		for dev in adt.deviceNames: data2D.append(np.zeros((len(listname),nbunch)));
		time2D=np.zeros(len(listname));
		xaveturns=np.zeros(nbunch);
		    

	# bunch bucket numbers (25ns buckets) for the selected bunches
	if (opt.BNUM==None):
    	    bunches=range(nbunch);
	else:
    	    bunches=np.array(opt.BNUM);
    	    bunches.sort();
	print "Selected bunches (begin at number 0):", bunches
	bunchesADT=[adt.bunchNumbers[j] for j in range(nbunch)];

	if (not(opt.ONEPLOT))or(i==0):
	    # initialize plots
	    ax=[];fig=[];axbeam=[];figbeam=[];
	    for dev in adt.deviceNames:
		fig1,ax1=init_figure()
		fig.append(fig1);ax.append(ax1);
		if (opt.AVER):
	            fig1,ax1=init_figure()
        	    figbeam.append(fig1);axbeam.append(ax1);

	# initialize tau, rise time files and plot
	if (opt.ANA):
	    filetau=[];filetaubeam=[];
	    for dev in adt.deviceNames:
		filetau.append(open(filename[:-5]+dev+opt.OUT,'w'));
		print >> filetau[-1], "Bunch\ttau[s]"
		if (opt.AVER):
	            filetaubeam.append(open(filename[:-5]+'_aver'+dev+opt.OUT,'w'));
		    print >> filetaubeam[-1], "tau[s]"
    	    tau=np.zeros(nbunch);
    	    figtau,axtau=init_figure();
	

	# pre-treatment:
	if (opt.PRE):
	    # take out any sudden (single-point) peak
	    pos=takeout_sudden_peak(turns,pos,opt.RATIO,nbunch,adt.deviceNames)		

	    # take out sliding average
            if (opt.SLIDE>0): turns,pos=substract_slideaverage(turns,pos,opt.SLIDE,nbunch)
		
	    # roll the data such that the instability is always at the end
	    if (opt.ROLL!=-1): pos=roll_data(turns,pos,nbunch,opt.ROLL)
		
	# end of data pre-treatment
	

        # main loop: analyse and plot for each device (4 in general)
	for j,x in enumerate(pos):

	    for bnum in range(nbunch-1,-1,-1):
				              
	    	x1=x[bnum];
	    
    	    	if (bnum in bunches):
		    # plot raw data for the bunches chosen
		    if (opt.ONEPLOT):
		        time=(timestamp+toff+turns*Trev)/86400.;
			if (i==0):
		            plotdate(time,x1,adt.deviceNames[j].replace("ADTBpos","")+', bunch '+str(bunchesADT[bnum]),dat,'-',"Average position of the bunch [a.u.]",ax[j],opt.LOG,plotevery=len(listname),colr=col[bnum]);
			else:
		            plotdate(time,x1,'',dat,'-',"Average position of the bunch [a.u.]",ax[j],opt.LOG,plotevery=len(listname),colr=col[bnum]);
			    
		    else:
    		    	plot(turns,x1,adt.deviceNames[j].replace("ADTBpos","")+', bunch '+str(bunchesADT[bnum]),'-',"Average position of the bunch [a.u.]",ax[j],opt.LOG);
			

		if (opt.AVER):
	    	    if (bnum==nbunch-1):
			xave=x1;
		    else:
			xave=x1+xave;	    

		# data analysis (fit)
    		if (opt.ANA):
		    flag=(bnum in bunches);
		    tau[bnum]=fit(turns,x1,ax[j],flag,adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]),"Average position of the bunch [a.u.]",Trev,opt.BOUND,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG);
		    print >> filetau[j], bnum, "\t", tau[bnum];
		    
		# for the 2D plot
		if (opt.ONEPLOT):
		    # computes the average absolute value vs turns for this bunch
		    xaveturns[bnum]=np.average(np.abs(x1));
		

	    # for the average over all bunches
	    if (opt.AVER):
		xave/=float(nbunch);	    
		# plot raw data (average)
		if (opt.ONEPLOT):
		    time=(timestamp+toff+turns*Trev)/86400.;
		    if (i==0):
		    	plotdate(time,xave,adt.deviceNames[j].replace("ADTBpos",""),dat,'.',"Average position of the "+str(nbunch)+"bunches [a.u.]",axbeam[j],opt.LOG,plotevery=len(listname),colr=col[bnum]);
		    else:
		    	plotdate(time,xave,'',dat,'.',"Average position of the "+str(nbunch)+"bunches [a.u.]",axbeam[j],opt.LOG,plotevery=len(listname),colr=col[bnum]);
		else:
    		    plot(turns,xave,adt.deviceNames[j].replace("ADTBpos",""),'.',"Average position of the "+str(nbunch)+"bunches [a.u.]",axbeam[j],opt.LOG);
		axbeam[j].set_title(fil);

		if (opt.ANA):
	    	    # fit
		    taubeam=fit(turns,xave,axbeam[j],True,adt.deviceNames[j],"Average position of the "+str(nbunch)+" bunches [a.u.]",Trev,opt.BOUNDBM,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG);
		    print adt.deviceNames[j]+": Average rise times of the "+str(nbunch)+" bunches in sec.: ", taubeam;
		    print >> filetaubeam[j], taubeam;

	    # for the 2D plot
	    if (opt.ONEPLOT):			
		time2D[i]=np.average((timestamp+toff+turns*Trev)/86400.);
		data2D[j][i,:]=xaveturns;
		


	    # finalize plots
	    if (opt.ANA):
    		filetau[j].close()
		axtau.plot(range(0,nbunch),tau,label=adt.deviceNames[j].replace("ADTBpos","")+' '+r'$\tau$');
		axtau.set_title(fil);
		axtau.set_xlabel("Bunch number");
		axtau.set_ylabel("Rise time [s]");
		axtau.legend(loc=0);
		set_fontsize(figtau,'large');

	    #ax[j].set_title(fil);
	    ax[j].legend(loc=0);
	    if (not(opt.ONEPLOT)): ax[j].set_title(fil);
	    #set_fontsize(fig[j],'large');
	    #set_fontsize(fig[j],'xx-large');
	    set_fontsize(fig[j],'medium');
	    if (opt.AVER):
		axbeam[j].legend(loc=0);
		if (not(opt.ONEPLOT)): axbeam[j].set_title(fil);
	        set_fontsize(figbeam[j],'large');

    
    if (opt.ONEPLOT):
	timeinterval=np.ceil((timestamp+nturns*Trev-t0)/8.);
	print "Time interval for plot: ", timeinterval
	for j,dev in enumerate(adt.deviceNames):
	    if (timeinterval<60):
    		ax[j].xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
		ax[j].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
	    elif (timeinterval<3600):
    		ax[j].xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=int(np.floor(timeinterval/60.))))
		ax[j].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
	    else:
    		ax[j].xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=int(np.floor(timeinterval/3600.))))
		ax[j].xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
		
	    # 2D plot of amplitudes, bunch-by-bunch
	    # re-interpolate data on a regular mesh (25ns slots for x axis)
	    timeintervalinterp=np.min(np.diff(time2D))/2;
	    newtime2D=np.arange(time2D[0],time2D[-1],timeintervalinterp);
	    newdata2D=np.zeros((len(newtime2D),3564));
	    for bnum,bunch in enumerate(bunchesADT):
	    	newdata2D[:,bunch]=np.interp(newtime2D,time2D,data2D[j][:,bnum]);
	    # plot
	    fig2D,ax2D=init_figure();
	    plot2D(newdata2D,0,3563,newtime2D[0],newtime2D[-1],'25ns slot number',
		'Local time on '+dat.strftime("%Y-%m-%d"),dev,ax2D);
	    ax2D.yaxis_date(tz=gmt);
    	    ax2D.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval))
	    ax2D.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
	    end_figure(fig2D,ax2D);
	
	

    pylab.show();

    sys.exit()

