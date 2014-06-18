#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
import pylab,re,dateutil,random,pytz
import numpy as np
from string import split, replace
from parser import *
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
from Timber import parseout
import math
import matplotlib
import matplotlib.dates
#import subprocess
import glob
from SussixNM import *
from plot_lib import init_figure,end_figure,cmap,plot,plot2D,make_movie,set_fontsize,build_colors
from io_lib import list_files,write_Timber
from datetime_lib import set_axisdate
from read_Headtail_prt_coupledbunch import gettune
from read_Headtail_prt_fit import slideaver,envelop,fit
from read_Headtail_prt_sussix import extract,plot_spec,collect_and_fit,intercalate
from read_BBQ_fit import read_BBQ_multi,plotextradata
from numpy import fft
from datetime import time,datetime,date
from collimator_settings import concatenate_data_Timber
from scipy import optimize
#from time import time


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--averaging",type=int,
                      help="Specify on how many successive points we average the tunes obtained (default=1, i.e. no averaging). Note: choosing a value different from 1 deactivates the 2D spectrum plot.",
                      metavar="AVER", default=1,dest="AVER")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify bounds for the vertical scale for the spectrum plot (default=[1.e4,1.e9])",
                      metavar="BOUND", default=[1.e4,1.e9],dest="BOUND")		      
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the fit (default=0)",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-i", "--interval",type=float,nargs=2,
                      help="Specify the maximum dispersion of the tune in an interval, and the minimum number of turns of such an interval to allow the average tune computation (this is e.g. to evaluate the tune shift related to certain events such as collimator movements). Default = 5e-5 and 300000.",
                      metavar="INTER", default=[5.e-5,300000],dest="INTER")
    parser.add_option("-j", "--legend",action="store_true",
                      help="Specify if add the legend (turn number) in spectrum plot",
                      metavar="LEG", default=False,dest="LEG")
    parser.add_option("-k", "--makemovie",action="store_true",
                      help="Specify if we make a movie with the sliding spectra (now not anymore movie but several spectra on same plot)",
                      metavar="MAKE", default=False,dest="MAKE")
    parser.add_option("-l", "--plot",type=int,
                      help="Specify 1 to plot also maxima, 2 to plot maxima and sliding average, 3 to plot sliding average (it also always plots raw data and fit)",
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
    parser.add_option("-q", "--tune",type=float,
                      help="Specify the fractional part of the tune (default=0.31).",
                      metavar="TUNE", default=0.31,dest="TUNE")
    parser.add_option("-r", "--fourier",action="store_true",
                      help="Specify if we use the FFT instead of Sussix",
                      metavar="FFT", default=False,dest="FFT")
    parser.add_option("-s", "--stopturn",type=int,
                      help="Specify the turn at which we stop the fit (default=stop at the end)",
                      metavar="STOP", default=None,dest="STOP")
    parser.add_option("-t", "--turns",type=int,nargs=2,
                      help="Specify number of turns to calculate the spectra on, and number of turns between two successive windows (default=8192 and 500)",
                      metavar="TURNS", default=[8192,500],dest="TURNS")
    parser.add_option("-u", "--boundmin",type=float,
                      help="Specify minimum amplitude of spectral lines to begin the fit (default=1e8)",
                      metavar="BMIN", default=1e8,dest="BMIN")
    parser.add_option("-v", "--subtractaverage",action="store_true",
                      help="Specify if we subtract the average value before computing the spectrum",
                      metavar="SUBAVER", default=False,dest="SUBAVER")
    parser.add_option("-w", "--width",type=float,
                      help="Specify the half width of the window into which we should find the tune (default=0.05)",
                      metavar="WIDTH", default=0.05,dest="WIDTH")
    parser.add_option("-x", "--extradata",action="callback",callback=multistring_parse,
                      help="Specify another TIMBER datafile and variable name (actually, only part of the full name, enough to identify it) to plot on the tune graph. File name can be a regular expression between quotes "" (for multifiles). Optionally, three additional arguments may be given (in this order): legend to put on the graph, label for the upper x-axis, and flag to say if we use this data to get the boundaries where to compute tune averages (with -i option) (if no 5th option => default=False, otherwise put 1=True).",
                      metavar="EXTRA", default=None,dest="EXTRA")
    parser.add_option("-y", "--precision",type=float,
                      help="Specify precision of the spectral lines interpolation for 2D plot (default=1e-5)",
                      metavar="PRECISION", default=1e-5,dest="PRECISION")
    parser.add_option("-z", "--erase",action="store_true",
                      help="Specify if we erase all _tmp*.png files after making the movie",
                      metavar="RM", default=False,dest="RM")
    parser.add_option("--chroma",type=float,nargs=2,
                      help="With -x option only, giving radial trims: specify if we try to recalculate chromaticity using the radial trims given in -x option. 2 parameters should be provided: frequency of the trim in Hz (usually 2.5 Hz) and amplitude of them in delta_p/p units (usually 2e-4).",
                      metavar="CHROMA", default=None,dest="CHROMA")
    parser.add_option("--clim",type=float,nargs=2,
                      help="Specify the limits of the color scale for the 2D plot (default: autoscale, which can be different for each plot)",
                      metavar="CLIM", default=None,dest="CLIM")
    parser.add_option("--cut",action="store_true",
                      help="Specify if we cut the FFT or Sussix analysis between the BEG and STOP turns (resp. -g and -s options).",
                      metavar="CUT", default=False,dest="CUT")
    parser.add_option("--nodatetuneamp",action="store_true",
                      help="Specify if we do NOT use the dates for the x-axis for the tune and amplitude plots (use turns instead). Useful when python gets stuck because of dates on the x-axis of the tune plot (unknown reason).",
                      metavar="NODATE", default=False,dest="NODATE")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    parser.add_option("--split",action="callback",callback=multifloatint_parse,
                      help="Specify if we split the data in several parts (vs time) which all have different factional tunes. It then cshould be in the form ""--split 200000 0.312 500000 0.311"" meaning that from the beginning till 200000 turns we use Q=opt.TUNE (see above), then from 200000 turns till 500000 turns Q=0.312, and finally Q=0.311 for the rest. PLEASE NOTE: USE THIS DIRTY WAY OF FILTERING ONLY IF ABSOLUTELY NECESSARY!",
                      metavar="SPLIT", default=None,dest="SPLIT")

    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args


def filteredfft(x,low,up,exparray=None):

    # computes the fft only from the indices low to up
    # NOTE: this is VERY SLOW (even when exparray is given) unless 'up' is VERY CLOSE to 'low'.
    y=np.zeros(up-low+1,dtype=complex);
    N=float(len(x));
    if (exparray==None):
    	exparray=np.zeros((up-low+1,len(x)),dtype=complex);
	j=np.arange(N);
        for ik,k in enumerate(np.arange(low,up+1,dtype=float)):
	    exparray[ik,:]=np.exp(-2.*pi*1j*j*k/N);
	
    for ik,k in enumerate(range(low,up+1)):
        y[ik] = np.sum(exparray[ik,:]*x);
    
    return y;
    
    
def scan_tunes(tunes,turns,maxdif,nturns,axtu,filename,addstring='',databreaks=None,turnsbreaks=None,xscaleplot=1.,xoffplot=0.):

    # scan all the tunes to find some interval where an average value can be computed,
    # then plot those average values
    # xscaleplot and xoffplot are resp. a scaling factor and an offset for the xaxis, in case of plot 
    # -> we multiply turns by xscaleplot, then add xoffplot (for the plot only)
    
    if (databreaks==None):
    	databreaks=tunes;
	turnsbreaks=turns;
    # find break points such that in between those points the dispersion of the data
    # in 'databreaks' is not larger than maxdif
    breaks=find_breaks(databreaks,turnsbreaks,turns,maxdif);
    
    # find which are the intervals spanning more than nturns turns, compute
    # the mean tune and standard deviation on these, write them in a file, and plot
    filetune=open(filename+'_tunes'+addstring+'.txt','w');
    print >> filetune, "\ttune\tsigma\tfirst_turn\tlast_turn"
    k=0;lab='Average tune value';tuneshifts=[];vartuneshifts=[];
    for ibr,br in enumerate(breaks):
    	if (ibr>0)and(turns[br]-turns[breaks[ibr-1]]>=nturns):
	    tuneaver=np.average(tunes[breaks[ibr-1]:br]);
	    sig=np.sqrt(np.average(tunes[breaks[ibr-1]:br]*tunes[breaks[ibr-1]:br])-tuneaver**2);
	    print >> filetune, "\t", tuneaver, "\t", sig, "\t", turns[breaks[ibr-1]], "\t", turns[br];
	    if (k>0):
	    	lab='';
		tuneshifts.append(abs(tuneaver-oldtuneaver)); # tuneshift in absolute value
		vartuneshifts.append((sig**2+oldsig**2)/2.); # variance
		
	    plot(turns[breaks[ibr-1]:br]*xscaleplot+xoffplot,tuneaver*np.ones(len(turns[breaks[ibr-1]:br])),lab,'-g','Tune (fractional part)',axtu,0);
	    k+=1;oldtuneaver=tuneaver;oldsig=sig;
	    
    
    filetune.close();
    
    # write final average tune shift
    filetuneshift=open(filename+'_tuneshift'+addstring+'.txt','w');
    tuneshifts=np.array(tuneshifts);vartuneshifts=np.array(vartuneshifts);
    tuneshift=np.average(tuneshifts);sigtuneshift=np.sqrt(np.average(vartuneshifts));
    print >> filetuneshift, "\tabs(tuneshift)\tsigma"
    print >> filetuneshift, "\t", tuneshift, "\t", sigtuneshift;
    filetuneshift.close();
	    
    return;

    
def find_breaks(databreaks,turnsbreaks,turns,maxdif):

    # find break points such that in between those points the dispersion of 
    # 'databreaks' is not larger than maxdif
    maxi=-1.e50;mini=1.e50;breaks=[0];
    for idata,data in enumerate(databreaks):
	if (data>maxi): maxi=data;
	if (data<mini): mini=data;
	if ((maxi-mini)>maxdif):
	    maxi=-1.e50;mini=1.e50;
	    breaks.append(idata);
    
    breaks.append(len(databreaks)-1);
    #print breaks
    
    # convert those breaks such that turnsbreaks[breaks]=turns[newbreaks]
    newbreaks=[];
    for br in breaks:
    	l=locate(turns,turnsbreaks[br]);
	if (turns[l-1]==turnsbreaks[br]): l -= 1;
        newbreaks.append(l);
    #print newbreaks
    
    return newbreaks;
  

def locate(data,point):

    # search a table ordered in ascending order 'data' (inspired from Numerical Recipes) 		    ***
    # Gives the position lprov (integer) in data, such that
    # data[lprov-1]<=point<data[lprov]

    il=0;iu=len(data);
    while ((iu-il)>1):
    	im=(iu+il)/2; # midpoint
   	if (point>=data[im]): il=im;
   	else: iu=im;

    if (point<=data[0]): lprov=1;
    elif (point>=data[-1]): lprov=len(data)-1; 
    else: lprov=iu;

    return lprov;


def diff_with_sin(param,x,y,amp):

    # computes difference between signal and a sine wave
    # absissae in x, signal in y, frequency (Hz) in freq, amplitude (y units) in amp.
    # Parameters (used to fit) are
    # 	param[0]: phase (in x units)
    # 	param[1]: offset (in y units)
    
    freq-param[2];
    eps=sum(1e8*(y-param[1]-amp*np.sin(2.*np.pi*freq*(x+param[0])))**2);
    return eps;
    
    
def chroma(t_tunes,tunes,t_trims,trims,t0,frev,freq,amp,ax,legtrim,labeltrim,beam,plane,output,nevery=10000):

    # plot tunes vs. radial trims and computes Q' and Q''
    # uses:
    # 	- (t_tunes, tunes): tunes vs time [s],
    # 	- (t_trims, trims): radial trims (deltap_p/p) vs time [s],
    #	- t0 & frev: initial time and revolution frequency
    #	- freq & amp: frequency and amplitude of the trims (we fit them with a sine),
    #	- ax, leg & label: axes on which to plot the fitted trims, legend and ylabel,
    #	- beam & plane: beam & plane studied,
    # 	- output: output file name for the chromaticity and Q'' (put data in a Timber-like way),
    #	- nevery: fit every this number of points of tunes
    
    # fit the radial trims with a sine wave
    offsetini=np.average(trims[:2*(len(trims)/2)]);# initial value for offset
    xini=[0.,offsetini,freq]
    #bnds=((0, 2./freq), (None, None), (freq*0.9, freq*1.1)) # bounds for minimization
    #xsol,fl,dic=optimize.fmin_l_bfgs_b(diff_with_sin,xini,approx_grad=True,args=(t_trims,trims,amp),
    #	bounds=bnds,maxfun=20000,factr=100.,m=100,iprint=0,pgtol=1.e-15,epsilon=1.e-12)
    #xsol=optimize.fmin(diff_with_sin,xini,args=(t_trims,trims,amp),maxfun=20000)
    xsol=optimize.fmin_powell(diff_with_sin,xini,args=(t_trims,trims,amp),maxfun=20000,xtol=1.e-10,ftol=1.e-10)
    phase=xsol[0];offset=xsol[1];freq=xsol[2];print phase,offset,offsetini,freq
    # recompute trims with the sine fit and t_tunes
    trims_fit=offset+amp*np.sin(2.*np.pi*freq*(t_tunes+phase));
    
    # plot the fitted trims vs time (same figure as tunes)
    plot((t_tunes-t0)*frev,trims_fit,legtrim+' (sine fit)','-g',labeltrim,ax,0);
    
    # plot (on separate figure) tune vs trims
    figvs,axvs=init_figure();
    #plot(trims_fit,tunes,'','.b','Tune (fractional part), '+beam+plane,axvs,0,xlab=legtrim);
    
    # 2nd order polynomial fit of tunes vs trims every "nevery" turns
    tQ=np.zeros(len(t_tunes[::nevery]));
    Q=np.zeros(len(t_tunes[::nevery]));
    Qprime=np.zeros(len(t_tunes[::nevery]));
    Qsec=np.zeros(len(t_tunes[::nevery]));
    col=build_colors(len(range(0,len(t_tunes),nevery)))
    for it,t in enumerate(range(0,len(t_tunes),nevery)):
    	#print len(trims_fit[it:it+nevery]),nevery;
	p=np.polyfit(trims_fit[t:t+nevery],tunes[t:t+nevery],2)
	tQ[it]=t_tunes[t];Q[it]=p[2];
	Qprime[it]=-p[1]; # WARNING : I don't know why the minus sign !
	Qsec[it]=2.*p[0];
    	plot(trims_fit[t:t+nevery],tunes[t:t+nevery],'','x','Tune (fractional part), '+beam+plane,axvs,0,xlab=legtrim,colr=col[it]);
    	plot(trims_fit[t:t+nevery],p[0]*trims_fit[t:t+nevery]**2+p[1]*trims_fit[t:t+nevery]+p[2],'','-','Tune (fractional part), '+beam+plane,axvs,0,xlab=legtrim,colr=col[it]);
	
    # write in a file Q, Q' and Q'' (Timber-like way)
    fileQ=open(output,'w');
    write_Timber(tQ,Q,output,"LHC.CALC.Q_"+beam+"_"+plane)
    write_Timber(tQ,Qprime,output,"LHC.CALC.Qprime_"+beam+"_"+plane)
    write_Timber(tQ,Qsec,output,"LHC.CALC.Qsec_"+beam+"_"+plane)

    end_figure(figvs,axvs);

    return tQ,Q,Qprime,Qsec;
    

if __name__ == "__main__":
    opt,args=parsse();
    
    if (opt.EXTRA!=None)and(len(opt.EXTRA)<2): 
    	print "too few arguments for -x option";
	sys.exit();

    # revolution frequency and period (assumes protons)
    #print 'Energy: ',opt.EN,' GeV';
    #gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    #beta=np.sqrt(1.-1./(gamma*gamma));
    #frev=beta*299792458./opt.CIRC; # the later is the circumference
    #Trev=1./frev;

    # define maximum number of files put in a single array
    #if (opt.AVER==1):
    datasplit=10;
    #else: datasplit=30; # choose higher number for tune shift computation

    Q=opt.TUNE
    Q=Q-floor(Q)
    #gmt=pytz.timezone('Europe/Amsterdam');
    
    if (opt.SPLIT!=None):
	Qsplit=[Q];
    	nsplit=[0];
    	for iarg,arg in enumerate(opt.SPLIT):
	    if (mod(iarg,2)==0): nsplit.append(arg);
	    else: Qsplit.append(arg);
	nsplit=np.array(nsplit);


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
    #for j in listname: print j;
    
    ndata=int(ceil(len(listname)/float(datasplit)));
    print ndata;
	
    filetaubeam=open(listname[0]+'_tau_'+beam+plane+opt.OUT+'.txt','w');
    print >> filetaubeam, "\ttau[s]\ttune\tsigma(tune)"
	
	
    # read all the BBQ data in listname
    tp=[];BBQ=[];lenBBQ=0;turnoffset=[];
    for idata in range(0,ndata):

    	if (idata<(ndata-1)):
	    #take one more datafile than datasplit (datasplit+1 files), to have some margin at the end
	    tpstmp,BBQtmp,frev=read_BBQ_multi(listname[idata*datasplit:(idata+1)*datasplit+1],opt.EN,beam,plane);
	
	else: tpstmp,BBQtmp,frev=read_BBQ_multi(listname[idata*datasplit:],opt.EN,beam,plane);

	if (idata==0): t0=tpstmp[0];
	# tps is the time in seconds
	tptmp=(tpstmp-t0)*frev; # now it's in turn numbers
	tp.append(tptmp);
	BBQ.append(BBQtmp);
	# first turn of each set
	turnoffset.append(int(tp[idata][0]));
	print idata,tp[idata][0],turnoffset[idata];
   
    Trev=1./frev;   	
    # total length
    lenBBQ=turnoffset[ndata-1]+len(BBQ[ndata-1]);
    turnoffset=np.array(turnoffset);
    print lenBBQ, turnoffset 
    
    
    if (opt.STOP==None): stopturn=lenBBQ;
    else: stopturn=opt.STOP;
    

    #initialize tune shift and amplitude vs turns plots
    figamp,axamp=init_figure();
    figtu,axtu=init_figure();

    # initialize tune and amplitude
    tunem=[];ampm=[];

    if (opt.MAKE):
	# initialize spectra plots
	figsp,axsp=init_figure();

    if (opt.AVER==1):
	# initialize 2D plots (amps vs tunes and turns) ("hump buster" like)
	# and initialize 2D data
	fig2D,ax2D=init_figure(axes=[0.16,0.1,0.8,0.8]);
	ntunes=ceil(2.*opt.WIDTH/opt.PRECISION);
	if (opt.CUT): data2D=np.zeros((len(range(opt.BEG,stopturn-opt.TURNS[0]+1,opt.TURNS[1])),ntunes),dtype=float);
        else: data2D=np.zeros((len(range(0,lenBBQ-opt.TURNS[0]+1,opt.TURNS[1])),ntunes),dtype=float);
 
    

#    if (opt.FFT):
#	# parameters for filtered FFT
#	N=float(opt.TURNS[0]);
#	low=floor((Q-opt.WIDTH)*N);
#	up=ceil((Q+opt.WIDTH)*N);
#    	exparr=np.zeros((up-low+1,N),dtype=complex);
#	j=np.arange(N);
#	for ik,k in enumerate(np.arange(low,up+1,dtype=float)):
#	    exparr[ik,:]=np.exp(-2.*pi*1j*j*k/N);

    if (opt.FFT):
	# parameters for fast sliding FFT
	N=float(opt.TURNS[0]);
	low=floor((Q-opt.WIDTH)*N);
	up=ceil((Q+opt.WIDTH)*N);
	k=np.arange(low,up+1,dtype=float);
	p=opt.TURNS[1];
	exparray=np.exp(2.*pi*1j*k*float(p)/N);
	j=np.arange(opt.TURNS[1],dtype=float);
    	exparray2=np.zeros((up-low+1,p),dtype=complex);
        for ik,k0 in enumerate(k):
	    exparray2[ik,:]=np.exp(-2.*pi*1j*j*k0/N);
    

    if (opt.CUT): turns=np.arange(opt.BEG,stopturn-opt.TURNS[0]+1,opt.TURNS[1]);
    else: turns=np.arange(0,lenBBQ-opt.TURNS[0]+1,opt.TURNS[1]);

    if (opt.FFT): print "   Sliding FFT ...";
    else: print "   Sliding Sussix ...";
    # loop on turns, to make a sliding window (analysed by sussix or sliding FFT)
    #for iturn,turn in enumerate(range(opt.BEG,stopturn-opt.TURNS[0],opt.TURNS[1])):
    for iturn,turn in enumerate(turns):

	if (not(opt.FFT))and(mod(turn-turns[0],10000)==0): print "Turn",turn
	
	if (opt.SPLIT==None): Qextr=Q;
	else:
	    # find around which tune we look for this turn (dirty filtering option)
    	    isplit=np.searchsorted(nsplit,turn,side='right');
	    if (isplit>0): isplit-=1;
	    Qextr=Qsplit[isplit];

	# find which sets this turn belongs to
	idata=np.where(turnoffset<=turn);idata=idata[0][-1];
	#print turn,idata,iturn;

	x2=BBQ[idata][turn-turnoffset[idata]:turn-turnoffset[idata]+opt.TURNS[0]];
	
	if (opt.SUBAVER): x2=x2-mean(x2)

	
	if (opt.FFT):

	    #if (mod(turn-opt.BEG,1000*opt.TURNS[1])==0)or(opt.TURNS[1]>opt.TURNS[0]/10.):
	    if (mod(turn,1000*opt.TURNS[1])==0)or(opt.TURNS[1]>opt.TURNS[0]/10.):
	        # do the normal fft from time to time to avoid error accumulation
		#st=time();
        	y2=fft.fft(x2);
		#print "turn %d, do a normal fft" % turn
        	#print "t1: ",time()-st
        	# sort and extract only tunes at opt.WIDTH from Qextr
        	tunes,amps=extract(np.arange(len(x2))/float(len(x2)),np.abs(y2),Qextr,opt.WIDTH)
		y2old=y2[low:up+1];
		x2old=x2;

	    else:
	        #st=time();
        	#y2=filteredfft(x2,low,up,exparray=exparr); # this is VERY SLOW
	        # fast sliding FFT - filtered over interesting tunes
		y2=exparray*(y2old+np.dot(exparray2,x2[-p:]-x2old[:p]));
		#print "t2: ",time()-st
        	#print max(np.abs(y3-y2[low:up+1])/np.abs(y2[low:up+1]))
        	tunes,amps=extract(np.arange(low,up+1)/float(len(x2)),np.abs(y2),Qextr,opt.WIDTH)
        	#print max(np.abs(amps2-amps)/np.abs(amps))
		y2old=y2;
		x2old=x2;
		
	    tunes2D=tunes;amps2D=amps;
	    
	else :
	    if (opt.AVER==1):
	        # 2D plots -> we need many lines in the spectrum
	    	suss=gettune(x2,np.zeros(len(x2)),Qextr,0.,0.07)
	    	# sort and extract only tunes at opt.WIDTH from Qextr
	    	tunes,amps=extract(suss.ox,suss.ax,Qextr,opt.WIDTH)
		# intercalate some zero values spaced by opt.WIDTH/2. if tunes & amps are too "empty"
		# (too avoid strange colored bands in 2D plots)
	    	tunes2D,amps2D=intercalate(tunes,amps,np.min(suss.ax),opt.WIDTH/2.,xmin=Qextr-opt.WIDTH,xmax=Qextr+opt.WIDTH);

	    else:
	    	# only tune (max. spectral line) is of interest -> take much less harmonics in spectrum
	    	nharm=0;flag=True;amps=0.;tunes=[-1.];
		while flag:
		    tuneold=tunes[np.argmax(amps)]
		    nharm+=5;#print "nharm=",nharm
		    suss=gettune(x2,np.zeros(len(x2)),Qextr,0.,0.07,nharm=nharm)
	    	    # sort and extract only tunes at opt.WIDTH from Qextr
	    	    tunes,amps=extract(suss.ox,suss.ax,Qextr,opt.WIDTH)
		    flag=False;
		    if (len(amps)==0):
		    	flag=True;amps=0.;tunes=[-1.];
		    elif (tunes[np.argmax(amps)]!=tuneold): flag=True;
		    
		if (nharm>20): print "nharm=",nharm;

 	if (opt.MAKE):
	    # plot spectrum
	    tim=pylab.num2date((t0+turn/frev)/86400.);
	    #plot_spec(tunes,amps,'Measurement '+beam+' '+plane+', '+tim.strftime("%Y-%m-%d %H:%M:%S"),'-',beam+plane+tim.strftime("%Y-%m-%d_%H:%M:%S"),axsp,figsp,Q,turn,opt.BOUND,width=opt.WIDTH,flagleg=opt.LEG);
	    plot_spec(tunes,amps,'Measurement '+beam+' '+plane+', '+tim.strftime("%Y-%m-%d %H:%M:%S"),'-',beam+plane+tim.strftime("%Y-%m-%d_%H:%M:%S"),axsp,figsp,Q,turn,opt.BOUND,lw=2.5,width=opt.WIDTH,flagleg=opt.LEG,flaglog=True,flagclear=False);

	if (opt.AVER==1):
	    # construct data for 2D plot
	    data2D[iturn,:]=np.interp(Q-opt.WIDTH+2.*opt.WIDTH*np.arange(ntunes)/float(ntunes),tunes2D,log(amps2D));

	# find maximum peak, and average with previous tune values obtained
	tunem.append(tunes[np.argmax(amps)])
	ampm.append(np.max(amps))
	    
    print "   Done";
			

    dat=pylab.num2date(t0/86400.).date();
    if (opt.AVER==1):
	# make the 2D plot ("water fall" or "hump buster" like: turns vs tunes with color=amplitudes)
	plot2D(data2D,Q-opt.WIDTH,Q+opt.WIDTH-2.*opt.WIDTH/float(ntunes),
    	    (t0+turns[0]/frev)/86400.,(t0+turn/frev)/86400.,'Tune',
	    'Local time on '+dat.strftime("%Y-%m-%d"),'',ax2D,
	    colorlabel="Log(spectrum amplitude)",colorlim=opt.CLIM,fig=fig2D);
    	timeinterval=ceil(((t0+turn/frev)-(t0+turns[0]/frev))/10.);
	set_axisdate(ax2D,'y',timeinterval);
    
    #if (opt.MAKE):
    	# make a movie with the spectra
	#make_movie(listname[0]+'.gif','_tmp',flagrm=opt.RM);
	#pylab.close(figsp)

    timeinterval=ceil(((t0+turn/frev)-(t0+turns[0]/frev))/5.);
    print "   Sliding averages ...";
    # do a sliding average
    #turns=slideaver(np.arange(opt.BEG,stopturn-opt.TURNS[0],opt.TURNS[1]),opt.AVER);turns=turns[:-opt.AVER];
    #turns=slideaver(np.arange(0,len(BBQ)-opt.TURNS[0],opt.TURNS[1]),opt.AVER);turns=turns[:-opt.AVER];
    turns=slideaver(turns+opt.TURNS[0]/2.,opt.AVER);turns=turns[:-opt.AVER];
    tunemaver=slideaver(tunem,opt.AVER);tunemaver=tunemaver[:-opt.AVER];
    ampmaver=slideaver(ampm,opt.AVER);ampmaver=ampmaver[:-opt.AVER];
    print "   Done";
    # collect tunes and amplitudes
    print "   Collecting tunes and tau ...";
    if not(opt.NODATE):
	tunebeam,taubeam,sigmabeam=collect_and_fit(tunemaver,ampmaver,turns,'BBQ '+beam+plane,
	    plane,opt.BMIN,True,axtu,axamp,Trev,'Tune (fractional part)',
	    'Amplitude of the tune line',firstturn=opt.BEG,lastturn=stopturn,
	    flagaver=True,xlab='Local time on '+dat.strftime("%Y-%m-%d"),xscaleplot=Trev/86400.,xoffplot=t0/86400.)
	set_axisdate(axtu,'x',timeinterval);
	set_axisdate(axamp,'x',timeinterval);
    else:
	tunebeam,taubeam,sigmabeam=collect_and_fit(tunemaver,ampmaver,turns,'BBQ '+beam+plane,
	    plane,opt.BMIN,True,axtu,axamp,Trev,'Tune (fractional part)',
	    'Amplitude of the tune line',firstturn=opt.BEG,lastturn=stopturn,
	    flagaver=True)
    print "Rise time in sec.: ", taubeam*Trev;
    print "Tune: ", tunebeam;
    print >> filetaubeam, "\t", taubeam*Trev, "\t", tunebeam, "\t", sigmabeam;
    filetaubeam.close();
    print "   Done";
    

    if (opt.EXTRA!=None):
	# plot another data file
	if (opt.CHROMA==None): patcol='-r';
	else: patcol='.r';
	
	if not(opt.NODATE):
	    tpsnew,datanew,var,leg,label,axtu2=plotextradata(opt.EXTRA,axtu,t0+turns[0]/frev,t0+turns[-1]/frev,flagtime=0,patcol=patcol);
	    set_axisdate(axtu2,'x',timeinterval);
	else:
	    tpsnew,datanew,var,leg,label,axtu2=plotextradata(opt.EXTRA,axtu,t0+turns[0]/frev,t0+turns[-1]/frev,flagtime=frev,patcol=patcol);
	
	if (opt.CHROMA!=None):
	    freq=opt.CHROMA[0];amp=opt.CHROMA[1];
	    tQ,Q,Qprime,Qsec=chroma(t0+turns/frev,tunemaver,tpsnew,datanew,t0+turns[0]/frev,frev,freq,amp,axtu2,leg,label,beam,plane,listname[0]+'_chroma_'+beam+plane+opt.OUT+'.txt',nevery=1000);
	    
	if (opt.AVER==1):
	    # also add extra plot on 2D tune plot
	    xlim=ax2D.get_xlim();
	    ylim=ax2D.get_ylim();
	    ax2D2=pylab.twiny(ax=ax2D); # second x axis
	    if (var.find('BCTDC')!=-1):
	    	if var.find('.B1'): leg='B1 intensity';
	    	elif var.find('.B2'): leg='B2 intensity';
	    	plot(datanew/1e13,tpsnew/86400.,leg,'-k','',ax2D2,0,lw=4.,xlab=r"Beam intensity $/10^{13}$");
	    #elif (var.find('BETASTAR')!=-1):
	    #	plot(datanew,tpsnew/86400.,opt.EXTRA[1],'-k','',ax2D2,0,lw=4.,xlab=r" $\beta^*$ [cm]");
	    else:
	    	plot(datanew,tpsnew/86400.,leg,'-k','',ax2D2,0,lw=4.,xlab=label);
	    
	    ax2D2.legend(loc=3);
	    set_axisdate(ax2D2,'y',timeinterval);
    	    ax2D.set_xlim(xlim);
    	    ax2D.set_ylim(ylim);
    	    ax2D2.set_ylim(ylim);
	
    	
    if (opt.EXTRA!=None)and((len(opt.EXTRA)>=5)and(opt.EXTRA[4]=="1")):
    	# extract tunes in intervals delimited by large variations of the additional variable
	# (typically collimator half-gap, but it can be something else)
    	print "   Scan ",opt.EXTRA[1]," to get interval without large variations, and get tune averages over these intervals...";
	scan_tunes(tunemaver,turns,opt.INTER[0],opt.INTER[1],axtu,listname[0],addstring='_'+beam+plane+opt.OUT,databreaks=datanew,turnsbreaks=(tpsnew-t0)*frev);
    else:  
        # scan the tune to get averages over intervals where it does not move too much
    	print "   Scan the tunes to get interval averages ...";
        scan_tunes(tunemaver,turns,opt.INTER[0],opt.INTER[1],axtu,listname[0],addstring='_'+beam+plane+opt.OUT)
    print "   Done";
    
    
    if opt.SAVE:
    
	if (opt.AVER==1): end_figure(fig2D,ax2D,save=listname[0]+'_waterfall_'+beam+plane+opt.OUT);
	end_figure(figamp,axamp,save=listname[0]+'_amp_'+beam+plane+opt.OUT);
	end_figure(figtu,axtu,legpos=1,save=listname[0]+'_tune_'+beam+plane+opt.OUT);
	if (opt.MAKE): end_figure(figsp,axsp,save=listname[0]+'_spec_'+beam+plane+opt.OUT);
    else:
        
	if (opt.AVER==1): end_figure(fig2D,ax2D);
	end_figure(figamp,axamp);
	end_figure(figtu,axtu,legpos=1);
	if (opt.MAKE): end_figure(figsp,axsp);

	pylab.show();

    sys.exit()

