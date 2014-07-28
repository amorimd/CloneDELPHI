#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);

import pylab,re,random,pytz
import numpy as np
from string import split, replace
from parser_lib import *
import math
import matplotlib
from multiturndata_riccardodemaria_modifNico import *
from SussixNM import *
from read_ADT_fit import roll_data,takeout_sudden_peak,substract_slideaverage,readADT
from plot_lib import init_figure,cmap,plot,end_figure,plot2D,make_movie
from io_lib import list_files
from read_Headtail_prt_fit import slideaver,envelop,fit
from read_Headtail_prt_sussix import extract,plot_spec,collect_and_fit,intercalate
from read_Headtail_prt_coupledbunch import gettune
from numpy import fft
from scipy import interpolate
from datetime import time,datetime,date


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 0) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the raw data plot (several -b options possible). Put a space between -b and list or number.",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify bounds for the vertical scale for the spectrum plot (default=[1.e-3,1.e3])",
                      metavar="BOUND", default=[1.e-3,1.e3],dest="BOUND")		      
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the ADT (.data) name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the analysis (default=2000)",
                      metavar="BEG", default=2000,dest="BEG")
    parser.add_option("-i", "--ratio",type=float,
                      help="Specify ratio above which sudden single-point peaks should be suppressed (when using pre-treatment) (default=1.3)",
                      metavar="RATIO", default=1.3,dest="RATIO")
    parser.add_option("-k", "--makemovie",action="store_true",
                      help="Specify if we make tune spectra plots and movies with them (now movies are deactivated)",
                      metavar="MAKE", default=False,dest="MAKE")
    parser.add_option("-l", "--clim",type=float,nargs=2,
                      help="Specify the limits of the color scale for the 2D plot (default: autoscale, which can be different for each plot)",
                      metavar="CLIM", default=None,dest="CLIM")
    parser.add_option("-n", "--end",type=int,
                      help="Specify the number of turns (from the beginning) below which we want to analyse (default=32768)",
                      metavar="END", default=32768,dest="END")
    parser.add_option("-o", "--output",help="Specify the end of the output filename for rise times (default=tausussix.txt)",
                      default="tausussix.txt",dest="OUT")
    parser.add_option("-p", "--nopretreatment",action="store_false",
                      help="Specify if a pre-treatment of the raw data SHOULD NOT be performed",
                      metavar="PRE", default=True,dest="PRE")
    parser.add_option("-q", "--tunes",type=float,nargs=2,
                      help="Specify the factional part of the tunes (x and y) (default=[0.28,0.31])",
                      metavar="TUNES", default=[0.28,0.31],dest="TUNES")
    parser.add_option("-r", "--fourier",action="store_true",
                      help="Specify if we use the FFT instead of Sussix",
                      metavar="FFT", default=False,dest="FFT")
    parser.add_option("-s", "--slideaver",type=int,
                      help="Specify period used for the sliding average subtracted BEFORE analysis (to get rid of constant offset + some low frequency oscillations). If 0, no sliding average done (default=20)..",
                      metavar="SLIDE", default=20,dest="SLIDE")	      
    parser.add_option("-t", "--turns",type=int,nargs=2,
                      help="Specify number of turns to calculate the spectra on, and number of turns between two successive windows (default=[8192,500]). NOTE: if first argument=number of turns recorded on each data set, then no tune nor rise time values are extracted (2D spectra, and 1D maximum amplitudes and tunes plotted along all files) (in that case 2nd number does not matter but should not be 0, and -g option should be 0)",
                      metavar="TURNS", default=[8192,500],dest="TURNS")
    parser.add_option("-u", "--boundmin",type=float,
                      help="Specify minimum amplitude of spectral lines to begin the fit (default=5)",
                      metavar="BMIN", default=5,dest="BMIN")
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot (and fit with -a option) the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--width",type=float,
                      help="Specify the half width of the window into which we should find the tune (default=0.05)",
                      metavar="WIDTH", default=0.05,dest="WIDTH")
    parser.add_option("-z", "--rolloffset",type=int,
                      help="Specify the additional offset used in the data rolling (such that instability at the end). -1 means no roll (default=10).",
                      metavar="ROLL", default=10,dest="ROLL")	      
    parser.add_option("--nodatetuneamp",action="store_true",
                      help="Specify if we do NOT use the dates for the x-axis for the tune and amplitude plots vs time (use turns instead). Useful when python gets stuck because of dates on the x-axis of the tune plot (unknown reason).",
                      metavar="NODATE", default=False,dest="NODATE")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    parser.add_option("--toff",type=int,
                      help="Specify time offset in seconds (default=7200 seconds), to get correct local time on the graphs.",
                      metavar="TOFF", default=7200,dest="TOFF")
    # Note: options -c, -e, -o, -l, -d, -r and -g used only when -a option activated
    # -i and -s options used only if -t option activated	      
    (opt, args) = parser.parse_args()
    #print "Selected files:", opt.FILE
    return opt, args



if __name__ == "__main__":
    opt,args=parsse(); 

    gmt=pytz.timezone('Europe/Amsterdam');
    
    # number of seconds to convert to gmt time (for plot...)
    toff=opt.TOFF;

    # revolution frequency and period (assumes protons)
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the circumference
    Trev=1./frev;
    
    Qx=opt.TUNES[0]
    Qy=opt.TUNES[1]
    Qx=Qx-floor(Qx)
    Qy=Qy-floor(Qy)
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if (opt.MAKE):
	figsp=[];axsp=[];
    	if (opt.AVER):
	    figbeamsp=[];axbeamsp=[];
    
    
    for ifile,filename in enumerate(listname):
    
    	# string for the title
        fil=filename.replace("_"," ").replace(".data","").replace("/"," ").replace("-"," ");
    	print "Selected file:", filename
	
    	# read file
	turns,pos,nbunch,nturns,adt,timestamp=readADT(filename);
	
	# bunch bucket numbers (25ns buckets) for the selected bunches
	if (opt.BNUM==None):
    	    bunches=range(nbunch);
	else:
    	    bunches=np.array(opt.BNUM);
    	    bunches.sort();bunches=list(bunches);

	print "Selected bunches (begin at number 0):", bunches
	bunchesADT=[adt.bunchNumbers[j] for j in range(nbunch)];


    	if ((ifile==0)and(len(listname)>0))and(opt.TURNS[0]==nturns):
	    # initialize 2D data
	    ntunes=ceil(2.*opt.WIDTH/4.e-6); # tune precision chosen for 2D plot=4e-6
	    data2D=np.zeros((len(pos),len(bunches),len(listname),ntunes),dtype=float);
	    t2D=np.zeros(len(listname),dtype=float);
	    tune1D=np.zeros((len(pos),len(bunches),len(listname)),dtype=float);
	    amp1D=np.zeros((len(pos),len(bunches),len(listname)),dtype=float);
	    fig2D=[];ax2D=[];

	if (opt.AVER):
	    # initialize plots for average (amplitude - tunes vs turns)
       	    if (opt.TURNS[0]!=nturns):
		figbeamamp,axbeamamp=init_figure()
            	figbeamtu,axbeamtu=init_figure()


       	if (opt.TURNS[0]!=nturns):
	    # initialize tau, rise time files and plot
	    filetau=[];filetaubeam=[];
	    for dev in adt.deviceNames:
		if (opt.AVER):
	            filetaubeam.append(open(filename[:-5]+'_aver'+dev+opt.OUT,'w'));
	            print >> filetaubeam[-1], "\ttau[s]\ttune"
		filetau.append(open(filename[:-5]+dev+opt.OUT,'w'));
		print >> filetau[-1], "Bunch\ttau[s]\ttune\tsigma"
    	
	    tau=np.zeros(nbunch);
	    tune=np.zeros(nbunch);
	    sigma=np.zeros(nbunch);
    	    figtau,axtau=init_figure();
    	    figtune,axtune=init_figure();
	

	# pre-treatment:
	if (opt.PRE):
	    # take out any sudden (single-point) peak
	    pos=takeout_sudden_peak(turns,pos,opt.RATIO,nbunch,adt.deviceNames)		

	    # take out sliding average
            if (opt.SLIDE>0): turns,pos=substract_slideaverage(turns,pos,opt.SLIDE,nbunch)
		
	    # roll the data such that the instability is always at the end
	    if (opt.ROLL!=-1): pos=roll_data(turns,pos,nbunch,opt.ROLL)
		
	

        # main loop: analyse (sussix or fft) and plot for each device (4 in general)
	for j,x in enumerate(pos):
	
	    print 'Device name: ',adt.deviceNames[j]

	    # choose appropriate plane (horizontal or vertical) and tune
	    if (adt.deviceNames[j].find('Hor')!=-1):
	    	Q=Qx;plane='x';
	    else: 
	    	Q=Qy;plane='y';
	    
	    if (opt.TURNS[0]!=nturns):
	    	#initialize tune shift and amplitude vs turns plots
    	    	figamp,axamp=init_figure();
    	    	figtu,axtu=init_figure();
	    
	    if (opt.TURNS[0]<nturns): turnmax=min(nturns-opt.TURNS[0],opt.END);
	    else: turnmax=1;

	    # loop on bunches
	    for bnum in range(0,nbunch):
				              
	    	x1=x[bnum];
		
		print 'bunch no ',bnum
		
		if (opt.AVER):
	    	    if (bnum==0):
			xave=x1;
		    else:
			xave=x1+xave;	    

		if (bnum in bunches):
		
		    bn=bunches.index(bnum);
		    
		    if (opt.MAKE)and(ifile==0):
			# initialize spectra plots
			figsptmp,axsptmp=init_figure();
			figsp.append(figsptmp);axsp.append(axsptmp);

    		    if ((ifile==0)and(len(listname)>0))and(opt.TURNS[0]==nturns):
			# initialize 2D plots (amps vs tunes and turns) ("hump buster" like)
			fig2Dtmp,ax2Dtmp=init_figure();
			fig2D.append(fig2Dtmp);ax2D.append(ax2Dtmp);

		    # initialize tune and amplitude
		    tunem=[];ampm=[];turnm=[];


		    # loop on turns, to make a sliding window (analysed by sussix or FFT)
		    for turn in range(opt.BEG,turnmax,opt.TURNS[1]):

			x2=x1[turn:turn+opt.TURNS[0]]

			if (opt.FFT):
			    y2=fft.fft(x2);
			    tunes,amps=extract(np.arange(len(x2))/float(len(x2)),np.abs(y2),Q,opt.WIDTH)
			    tunes2D=tunes;amps2D=amps
			
			else:
			    suss=gettune(x2,np.zeros(len(x2)),Q,0.,0.07)
			    # sort and extract only tunes at 0.05 from nominal tune Q
			    tunes,amps=extract(suss.ox,suss.ax,Q,opt.WIDTH)
			    # intercalate some zero values spaced by opt.WIDTH/2. if tunes & amps are too "empty"
			    # (too avoid strange colored bands in 2D plots)
			    tunes2D,amps2D=intercalate(tunes,amps,np.min(suss.ax),opt.WIDTH/2.,xmin=Q-opt.WIDTH,xmax=Q+opt.WIDTH);

			if (opt.MAKE):
			    # plot spectrum for the bunches chosen
			    #plot_spec(tunes,amps,fil+' '+adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]),'-',filename+'_'+adt.deviceNames[j]+'_b'+str(bnum),axsp[j*len(bunches)+bn],figsp[j*len(bunches)+bn],Q,turn,opt.BOUND,flagleg=False,flaglog=True);
			    plot_spec(tunes,amps,adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]),'-','',axsp[j*len(bunches)+bn],figsp[j*len(bunches)+bn],Q,turn,opt.BOUND,flagleg=True,flaglog=True,flagclear=False,leg=datetime.fromtimestamp(timestamp+toff,tz=gmt).ctime());

			if (len(listname)>0)and(opt.TURNS[0]==nturns):
	    		    # construct data for 2D plot
			    t2D[ifile]=timestamp;
			    if (len(amps)>0):
	    			data2D[j,bn,ifile,:]=np.interp(Q-opt.WIDTH+2.*opt.WIDTH*np.arange(ntunes)/float(ntunes),tunes2D,amps2D);
				tune1D[j,bn,ifile]=tunes[np.argmax(amps)];
				amp1D[j,bn,ifile]=np.max(amps);
			    else:
			    	tune1D[j,bn,ifile]=nan;amp1D[j,bn,ifile]=nan;

			if (opt.TURNS[0]!=nturns)and(len(amps)>0):
			    # find maximum peak
			    tunem.append(tunes[np.argmax(amps)])
			    ampm.append(np.max(amps))
			    turnm.append(turn)


		    #if (opt.MAKE):
			#make_movie(filename[:-5]+'_'+adt.deviceNames[j]+'_b'+str(bunchesADT[bnum])+'.gif','_tmp');
			#pylab.close(figsp)


		    # collect tunes and amplitudes
		    if (opt.TURNS[0]!=nturns)and(len(tunem)>0):
			tune[bnum],tau[bnum],sigma[bnum]=collect_and_fit(tunem,ampm,turnm,'bunch '+str(bunchesADT[bnum]),
				plane,opt.BMIN,(bnum in bunches),axtu,axamp,Trev,'Tune (fractional part)',
				'Amplitude of the tune line',firstturn=0,lastturn=None,flagaver=True)
			print >> filetau[j], bunchesADT[bnum], "\t", tau[bnum]*Trev, "\t",tune[bnum], "\t", sigma[bnum];

		

	    # for the average over all bunches
	    if (opt.AVER):
	    
	        print 'Beam average'
		xave/=float(nbunch);	    
		
                if (ifile==0)and(opt.MAKE):
		    figbeamsptmp,axbeamsptmp=init_figure();
		    figbeamsp.append(figbeamsptmp);axbeamspappend(axbeamsptmp);
		
		# initialize tune and amplitude
		tunem=[];ampm=[];turnm=[];

		# loop on turns, to make a sliding window (analysed by sussix or FFT)
		for turn in range(opt.BEG,turnmax,opt.TURNS[1]):

		    xave2=xave[turn:turn+opt.TURNS[0]]

		    if (opt.FFT):
			yave2=fft.fft(xave2);
			tunes,amps=extract(np.arange(len(xave2))/float(len(xave2)),np.abs(yave2),Q,opt.WIDTH)

		    else:
			suss=gettune(xave2,np.zeros(len(xave2)),Q,0.,0.07)
			# sort and extract only tunes at 0.05 from nominal tune Q
			tunes,amps=extract(suss.ox,suss.ax,Q,opt.WIDTH)

		    if (opt.MAKE):
			# plot spectrum for the whole beam
			#plot_spec(tunes,amps,fil+' '+adt.deviceNames[j]+', average over all bunches','-',filename+'_'+adt.deviceNames[j]+'_aver',axbeamsp[j],figbeamsp[j],Q,turn,opt.BOUND,flagleg=False,flaglog=True);
			plot_spec(tunes,amps,adt.deviceNames[j]+', average over all bunches','-','',axbeamsp[j],figbeamsp[j],Q,turn,opt.BOUND,flagleg=True,flaglog=True,flagclear=False,leg=datetime.fromtimestamp(timestamp+toff,tz=gmt).ctime());

		    if (opt.TURNS[0]!=nturns)and(len(amps)>0):
			# find maximum peak
			tunem.append(tunes[np.argmax(amps)])
			ampm.append(np.max(amps))
			turnm.append(turn)
			

		# collect tunes and amplitudes
		if (opt.TURNS[0]!=nturns)and(len(tunem)>0):
		    tunebeam,taubeam,sigmabeam=collect_and_fit(tunem,ampm,turnm,adt.deviceNames[j],
			    plane,opt.BMIN,True,axbeamtu,axbeamamp,Trev,'Tune (fractional part)',
			    'Amplitude of the tune line',firstturn=0,lastturn=None,flagaver=True)
		    print adt.deviceNames[j]+": Average rise time of the "+str(nbunch)+" bunches in sec.: ", taubeam*Trev;
		    print adt.deviceNames[j]+": Average tune of the "+str(nbunch)+" bunches: ", tunebeam;
		    print >> filetaubeam[j], "\t", taubeam*Trev, "\t", tunebeam, "\t", sigmabeam;



	    if (opt.TURNS[0]!=nturns):
		# finalize plots
    		filetau[j].close()

		# plot tau vs bunches
		axtau.plot(bunchesADT,tau*Trev,label=adt.deviceNames[j]);
		axtau.set_title(fil);
		axtau.set_xlabel("Bunch number");
		axtau.set_ylabel("Rise time [s]");
		axtau.legend(loc=0);
		end_figure(figtau,axtau,save=opt.SAVE*(filename.replace(".data","")+adt.deviceNames[j]+"_tau_vs_bunches"));

		# plot tune vs bunches
		axtune.plot(bunchesADT,tune,label=adt.deviceNames[j]);
		axtune.set_title(fil);
		axtune.set_xlabel("Bunch number");
		axtune.set_ylabel("Tune");
		axtune.legend(loc=0);
		end_figure(figtune,axtune,save=opt.SAVE*(filename.replace(".data","")+adt.deviceNames[j]+"_tune_vs_bunches"));

		axtu.set_title(fil+' '+adt.deviceNames[j]+', '+plane);
		axamp.set_title(fil+' '+adt.deviceNames[j]+', '+plane);
		axtu.legend(loc=0);axamp.legend(loc=0);
		end_figure(figamp,axamp,save=opt.SAVE*(filename.replace(".data","")+adt.deviceNames[j]+"_amp_bunch"+str(bunchesADT[bnum])));
		end_figure(figtu,axtu,save=opt.SAVE*(filename.replace(".data","")+adt.deviceNames[j]+"_tune_bunch"+str(bunchesADT[bnum])));


	if (opt.AVER)and(opt.TURNS[0]!=nturns):
	    # finalize plots for average
	    axbeamtu.set_title(fil+', '+plane+', average over all bunches');
	    axbeamamp.set_title(fil+', '+plane+', average over all bunches');
	    axbeamamp.legend(loc=0);axbeamtu.legend(loc=0);
	    end_figure(figbeamamp,axbeamamp,save=opt.SAVE*(filename.replace(".data","")+"_amp_aver"));
	    end_figure(figbeamtu,axbeamtu,save=opt.SAVE*(filename.replace(".data","")+"_tune_aver"));

    
    if (len(listname)>0)and(opt.TURNS[0]==nturns):
    
	# make the 2D plots ("hump buster" like: turns vs tunes with color=amplitudes)
	# for each device and bunch

	dat=datetime.fromtimestamp(t2D[0],tz=gmt).date();
	t1=pylab.date2num(datetime.fromtimestamp(t2D[0],tz=gmt));
	t2=pylab.date2num(datetime.fromtimestamp(t2D[-1],tz=gmt));
	tim1D=[pylab.date2num(datetime.fromtimestamp(t0,tz=gmt))+toff/86400. for t0 in t2D];
	
	for j,x in enumerate(pos):
	
	    # choose appropriate plane (horizontal or vertical) and tune
	    if (adt.deviceNames[j].find('Hor')!=-1): Q=Qx;plane='x';
	    else: Q=Qy;plane='y';

	    for bn,bnum in enumerate(bunches):
	    
		# 2D "water-fall" plots
		f=interpolate.interp1d(t2D,data2D[j,bn],axis=0);
		data2Dnew=f(np.arange(t2D[0],t2D[-1],1.));
		plot2D(data2Dnew,Q-opt.WIDTH,Q+opt.WIDTH-2.*opt.WIDTH/float(ntunes),
    			t1+toff/86400.,t2+toff/86400.,'Tune','Local time on '+dat.strftime("%Y-%m-%d"),
			adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]),ax2D[j*len(bunches)+bn],
			colorlabel='Spectrum amplitude',colorlim=opt.CLIM,fig=fig2D[j*len(bunches)+bn]);
		ax2D[j*len(bunches)+bn].yaxis_date(tz=gmt);
		timeinterval=int(np.ceil((t2D[-1]-t2D[0])/8.));
		ax2D[j*len(bunches)+bn].yaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval))
		ax2D[j*len(bunches)+bn].yaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'));
    	
    		end_figure(fig2D[j*len(bunches)+bn],ax2D[j*len(bunches)+bn],save=opt.SAVE*(listname[0].replace(".data","")+adt.deviceNames[j]+"_waterfall_bunch"+str(bunchesADT[bnum])));
		
		# 1D tune and max amplitude plot
		if (len(~np.isnan(tune1D[j,bn,:]))>0):
		    timeinterval=int(np.ceil((t2D[-1]-t2D[0])/5.));
		    fig1,ax1=init_figure();fig2,ax2=init_figure();
		    if not(opt.NODATE):
			plot(tim1D,tune1D[j,bn,:],'','-xb','Tune',ax1,0,lw=4.,xlab='Local time on '+dat.strftime("%Y-%m-%d"));
			plot(tim1D,amp1D[j,bn,:],'','-xb','Amplitude of main tune line',ax2,0,lw=4.,xlab='Local time on '+dat.strftime("%Y-%m-%d"));
			ax1.set_title(adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]))
			ax2.set_title(adt.deviceNames[j]+', bunch '+str(bunchesADT[bnum]))
			ax1.xaxis_date(tz=gmt);
			ax1.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval));
			ax1.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'));
			ax2.xaxis_date(tz=gmt);
			ax2.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval));
			ax2.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'));
		    else:
			plot((np.array(tim1D)-tim1D[0])*86400*frev,tune1D[j,bn,:],'','-xb','Tune',ax1,0,lw=4.,xlab='Number of turns');
			plot((np.array(tim1D)-tim1D[0])*86400*frev,amp1D[j,bn,:],'','-xb','Amplitude of main tune line',ax2,0,lw=4.,xlab='Number of turns');
		    end_figure(fig1,ax1,save=opt.SAVE*(listname[0].replace(".data","")+adt.deviceNames[j]+"_tune_bunch"+str(bunchesADT[bnum])));
		    end_figure(fig2,ax2,save=opt.SAVE*(listname[0].replace(".data","")+adt.deviceNames[j]+"_amp_bunch"+str(bunchesADT[bnum])));
		    	
    
#    if (opt.MAKE):
#	# make movies of the spectra, for each bunch and each device
#	for j,x in enumerate(pos):
#	    
#	    for bnum in range(0,nbunch):
#		make_movie(adt.deviceNames[j]+'_b'+str(bunchesADT[bnum])+'.gif','_tmp'+filename+'_'+adt.deviceNames[j]+'_b'+str(bnum));
#		pylab.close(figsp)
#
#	    if (opt.AVER):
#		# make movies of the spectra, for averages
#		make_movie(adt.deviceNames[j]+'_aver.gif','_tmp'+filename+'_'+adt.deviceNames[j]+'_aver');
#		pylab.close(figbeamsp)

    
    
    if (opt.MAKE):
	for j,x in enumerate(pos):
	    for bn,bnum in enumerate(bunches):
    		end_figure(figsp[j*len(bunches)+bn],axsp[j*len(bunches)+bn],save=opt.SAVE*(listname[0].replace(".data","")+adt.deviceNames[j]+"_spectrum_bunch"+str(bunchesADT[bnum])));
	    
	    if (opt.AVER): end_figure(figbeamsp[j],axbeamsp[j],save=opt.SAVE*(listname[0].replace(".data","")+adt.deviceNames[j]+"_spectrum_aver"));
	    

    if not(opt.SAVE): pylab.show();

    sys.exit()

