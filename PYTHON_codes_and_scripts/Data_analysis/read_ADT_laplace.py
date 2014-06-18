#!/usr/bin/python2.6

import sys
#sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
import pylab,re,random
import numpy as np
from string import split, replace
from parser import *
import math
import matplotlib
from multiturndata_riccardodemaria_modifNico import *
from SussixNM import *
from plot_lib import init_figure,cmap,plot,set_fontsize
from read_ADT_fit import roll_data,takeout_sudden_peak,substract_slideaverage,readADT
from read_Headtail_prt_fit import slideaver,envelop,fit
from read_Headtail_prt_laplace import Laplace,find_peak_lapl,plot_spec_comp
from read_Headtail_prt_sussix import plot_spec

def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 0) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the plot and analysis (several -b options possible). Put a space between -b and list or number.",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--end",type=int,
                      help="Specify the number of turns from the beginning below which we want to do the analysis (default=32768)",
                      metavar="END", default=32768,dest="END")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the ADT .data name of the file to analyse (several files possible)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the analysis (default=2000)",
                      metavar="BEG", default=2000,dest="BEG")
    parser.add_option("-i", "--ratio",type=float,
                      help="Specify ratio above which sudden single-point peaks should be suppressed (when using pre-treatment) (default=1.3)",
                      metavar="RATIO", default=1.3,dest="RATIO")
    parser.add_option("-n", "--nmodes",type=int,
                      help="Specify number of modes to extract (default=5)",
                      metavar="NMOD", default=5,dest="NMOD")
    parser.add_option("-o", "--output",help="Specify the end of the output filename for rise times (default=tauLaplace.txt)",
                      default="tauLaplace.txt",dest="OUT")
    parser.add_option("-p", "--nopretreatment",action="store_false",
                      help="Specify if a pre-treatment of the raw data SHOULD NOT be performed",
                      metavar="PRE", default=True,dest="PRE")
    parser.add_option("-q", "--tunes",type=float,nargs=2,
                      help="Specify the factional part of the tunes (x and y) (default=[0.28,0.31])",
                      metavar="TUNES", default=[0.28,0.31],dest="TUNES")
    parser.add_option("-s", "--slideaver",type=int,
                      help="Specify period used for the sliding average subtracted BEFORE analysis, when using pre-treatment (to get rid of constant offset + some low frequency oscillations). If 0, no sliding average done (default=20).",
                      metavar="SLIDE", default=20,dest="SLIDE")		      
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also analyse the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--width",type=float,
                      help="Specify the half width of the window into which we should find the tune (default=0.05)",
                      metavar="WIDTH", default=0.05,dest="WIDTH")
    parser.add_option("-z", "--rolloffset",type=int,
                      help="Specify the additional offset used in the data rolling (such that instability at the end), when using pre-treatment. -1 means no roll (default=10).",
                      metavar="ROLL", default=10,dest="ROLL")		      
    # Note: options -c, -e, -o, -l, -d, -r and -g used only when -a option activated
    # -i and -s options used only if -t option activated	      
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    return opt, args



if __name__ == "__main__":
    opt,args=parsse(); 

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
    
    for filename in opt.FILE:
    
    	# string for the title
        fil=filename.replace("_"," ").replace(".data","").replace("/"," ").replace("-"," ");
	
    	# read file
	turns,pos,nbunch,nturns,adt,timestamp=readADT(filename);
	
	# bunch bucket numbers (25ns buckets) for the selected bunches
	if (opt.BNUM==None):
    	    bunches=range(nbunch);
	else:
    	    bunches=np.array(opt.BNUM);
    	    bunches.sort();

	print "Selected bunches (begin at number 0):", bunches
	bunchesADT=[adt.bunchNumbers[j] for j in range(nbunch)];


	if (opt.AVER):
	    # initialize plots for average(spectrum)
            figbeamsp,axbeamsp=init_figure()


	# initialize tau, rise time files and plots
	filetaubeam=[];filetau=[];
	for dev in adt.deviceNames:
	    if (opt.AVER):
		filetaubeam.append(open(filename[:-5]+'_aver'+dev+opt.OUT,'w'));
		print >> filetaubeam[-1], "\ttau[s]\ttune\tAmp"
	    
	    filetau.append(open(filename[:-5]+dev+opt.OUT,'w'));
	    print >> filetau[-1], "Bunch\ttau[s]\ttune\tAmp"

    	
	tau=np.zeros(nbunch);
	tune=np.zeros(nbunch);
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
		
	

        # main loop: analyse and plot for each device (4 in general)
	for j,x in enumerate(pos):
	
	    print 'Device name: ',adt.deviceNames[j]

	    # choose appropriate plane (horizontal or vertical) and tune
	    if (adt.deviceNames[j].find('Hor')!=-1):
	    	Q=Qx;plane='x';
	    else: 
	    	Q=Qy;plane='y';
	    
	    # loop on bunches
	    for bnum in range(0,nbunch):
				              
	    	x1=x[bnum];
		
		print 'bunch no ',bnum
		
		begt=np.where(turns>=opt.BEG);begt=begt[0]; # for some reason begt itself is a single-element array of array
		endt=np.where(turns>=opt.END);
		if (len(endt[0])>0): endt=endt[0]; # for some reason endt itself is a single-element array of array
		else: endt=[len(x1)];
		x2=x1[begt[0]:endt[0]];
		#print begt[0],len(turns),len(x1),len(x2)
		
		if (opt.AVER):
	    	    if (bnum==0):
			xave=x2;
		    else:
			xave=x2+xave;	    

		
		if (bnum in bunches):
		
		    # complex frequency spectral analysis
		    tunes,rates,amps=Laplace(x2,opt.NMOD);

		    # find main spectral line
		    tune[bnum],ra,amp=find_peak_lapl(tunes,rates,amps,Q,opt.WIDTH*2.,0)
		    tune[bnum]=tune[bnum]+Q;
		    if (ra!=0): tau[bnum]=1./ra;
		    print >> filetau[j], bunchesADT[bnum], "\t", tau[bnum]*Trev, "\t", tune[bnum], "\t", abs(amp);


    	    	#if (bnum in bunches):
		    # plot spectrum for the bunches chosen
		    #figsp,axsp=init_figure();
	            #axsp.plot(np.sort(tunes),np.abs(amps[np.argsort(tunes)]),'-',lw=2.5);
			

		

	    # for the average over all bunches
	    if (opt.AVER):
	    
	        print 'Beam average'
		xave/=float(nbunch);	    
		
		# complex frequency spectral analysis
		tunes,rates,amps=Laplace(xave,opt.NMOD);

		# find main spectral line
		tunebeam,ra,ampbeam=find_peak_lapl(tunes,rates,amps,Q,opt.WIDTH*2.,0)
		if (ra!=0): taubeam=1./ra;
		else: taubeam=0.;
		print >> filetaubeam[j], "\t", taubeam*Trev, "\t", tunebeam+Q, "\t", abs(ampbeam);

		# plot spectrum for the whole beam
        	figbeamsp,axbeamsp=init_figure()
	        plot_spec_comp(tunes,rates,amps,axbeamsp,"Average, "+adt.deviceNames[j]);

		print adt.deviceNames[j]+": Average rise time of the "+str(nbunch)+" bunches in sec.: ", taubeam*Trev;
		print adt.deviceNames[j]+": Average tune of the "+str(nbunch)+" bunches: ", tunebeam;
		
		filetaubeam[j].close()


	    # finalize
    	    filetau[j].close()
	    
	    # plot tau vs bunches
	    axtau.plot(bunchesADT,tau*Trev,label=adt.deviceNames[j]);
	    axtau.set_title(fil);
	    axtau.set_xlabel("Bunch number");
	    axtau.set_ylabel("Rise time [s]");
	    axtau.legend(loc=0);
	    set_fontsize(figtau,'xx-large');

	    # plot tune vs bunches
	    axtune.plot(bunchesADT,tune,label=adt.deviceNames[j]);
	    axtune.set_title(fil);
	    axtune.set_xlabel("Bunch number");
	    axtune.set_ylabel("Tune");
	    axtune.legend(loc=0);
	    set_fontsize(figtune,'xx-large');



    pylab.show();

    sys.exit()

