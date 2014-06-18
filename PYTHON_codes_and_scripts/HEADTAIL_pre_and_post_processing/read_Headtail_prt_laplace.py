#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
from parser import *
from string import *
import numpy as np
from numpy import fft
import pylab,os,re
from scipy import optimize
from read_cfg import read_cfg
from plot_lib import plot,init_figure,end_figure,plot2D
from io_lib import list_files
from read_Headtail_prt import extract_Headtail_param
from string_lib import takeout_common,takeout_spaces
from read_Headtail_prt_fit import read_prt_file,slideaver


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the spectrum plot (several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-e", "--every",type=int,
                      help="Specify the number of bunches to skip in the analysis- i.e. we analyse only one every 'this number' bunch (default=1)",
                      metavar="EVERY", default=1,dest="EVERY")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the fit (default=0)",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-i", "--bunchmin",type=int,
                      help="Specify the minimum bunch number to analyse (we analyse the bunches from nbunch - last bunch - to this bunch number) (default=1)",
                      metavar="NBMIN", default=1,dest="NBMIN")
    parser.add_option("-l", "--legend",action="append",
                      help="Specify the legend for the plot, for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-m", "--mminmmax",type=int,nargs=2,
                      help="Specify minimum and maximum headtail mode to find (default=[0,1])",
                      metavar="MMAX", default=[0,1],dest="MMAX")
    parser.add_option("-n", "--nmodes",type=int,nargs=2,
                      help="Specify number of modes to extract (x and y) (default=[10,10])",
                      metavar="NMOD", default=(10,10),dest="NMOD")
    parser.add_option("-o", "--output",help="Specify output suffix for rise times and tune shifts (default=tau.txt)",
                      default="tau.txt",dest="OUT")
    parser.add_option("-p", "--plot",action="store_true",
                      help="Specify if we plot the successive signals after each subtraction of harmonic, and the final spectra",
                      metavar="PLOT", default=False,dest="PLOT")
    parser.add_option("-s", "--stopturn",type=int,
                      help="Specify the turn at which we stop analysing (default=10000000 - if it's higher than the total number of turns, then we stop analysis at the maximum of the raw data)",
                      metavar="STOP", default=10000000,dest="STOP")
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot (and analyse with -a option) the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--width",type=float,
                      help="Specify the spectrum width around the unperturbed tune (default=0.02)",
                      metavar="WIDTH", default=0.02,dest="WIDTH")
    (opt, args) = parser.parse_args()
    return opt, args


#def plot2D(data,nxmin,nxmax,nymin,nymax,xlab,ylab,tit,ax):
#
#    # Creating image
#    # NOTE: imshow plots the matrix "as it looks like when you read it" ->
#    # this means that the lines with higher indices are at the BOTTOM of the graph
#    # (contrary to what would be more logical, i.e. higher indices at the top). Columns
#    # follow a more logical behaviour, i.e. higher indices at the right.
#    # Therefore, one needs to flip "upside-down" the date in the vertical directions,
#    # to reverse the orders of the lines -> flipud command
#
#   # NOTE 2: data is of the format [iy,ix] (i.e. y axis in first dimension)
#    
#    ax.imshow(np.flipud(data), aspect='auto', extent=[nxmin,nxmax,nymin,nymax])
#    # Plotting contour lines
#    #ax.contour(data, aspect='auto', extent=[nxmin,nxmax,nymin,nymax])
#    ax.set_xlabel(xlab)
#    ax.set_ylabel(ylab)
#    ax.set_title(tit)



def fftshifted(alpha,y):

    # compute the fft of y*exp(2*pi*1j*alpha*[0:length(y)-1])
    y2=y*np.exp(2*np.pi*1j*alpha*np.arange(len(y)));
    return fft.fft(y2);


def maxfftshifted(alpha,y,imin,imax):

    # find the max in absolute value of the spectrum computed by fftshifted
    # between indices imin and imax, and returns its opposite
    x=np.abs(fftshifted(alpha,y));
    return -np.max(x[imin:imax]);
    

def weightedlaplacespectrum(sigma,x,y):

    # see also "Laplace spectrum exponential decomposition and pole-zero estimation", by M.J. Corinthios,
    # IEEE Proc.-Vis. Image Signal Process., Vol. 148, No. 5, October 2001
    y2=np.exp(-sigma*x)*y;
    y3=fft.fft(y2);
    eps=sum(np.abs(y2)**2)/float(len(y2));
    return np.abs(y3)**2/eps;


def weightedlaplacespectrum_max(sigma,x,y):

    return -np.max(weightedlaplacespectrum(sigma,x,y));


def weightedlaplacespectrum2D(param,x,y):

    # see also "Laplace spectrum exponential decomposition and pole-zero estimation", by M.J. Corinthios,
    # IEEE Proc.-Vis. Image Signal Process., Vol. 148, No. 5, October 2001
    sigma=param[0];alpha=param[1];
    y2=np.exp(-sigma*x)*y;
    y3=fftshifted(alpha,y2);
    eps=sum(np.abs(y2)**2)/float(len(y2));
    return np.abs(y3)**2/eps;


def weightedlaplacespectrum2D_max(param,x,y):

    return -np.max(weightedlaplacespectrum2D(param,x,y));


def plot_weightedlaplacespectrum2D(x,y,sigmin=0.,sigmax=1e-4,sigmaint=1e-7):

    # 2D plot of a weighted spectrum vs sigma and alpha
    sigma=np.arange(sigmin,sigmax,sigmaint);
    alpmin=0.;alpmax=1.;
    #alpha=np.arange(alpmin,alpmax,1e-7);
    w=np.zeros((len(sigma),len(x[::50])));
    for isig,sig in enumerate(sigma):
    	#for ialp,alp in enumerate(alpha):
	    #w[isig,ialp]=weightedlaplacespectrum2D_max((sig,alp),x,y);
	wtmp=weightedlaplacespectrum(sig,x,y);
	w[isig,:]=wtmp[::50];

    fig,ax=init_figure();
    plot2D(w,alpmin,alpmax,sigmin,sigmax,"Tune"," $ \sigma $ ","2D weighted spectrum",ax);
    end_figure(fig,ax);


def plot_spectrum2D(x,y,sigmin=0.,sigmax=1e-4,sigmaint=1e-6,tunemin=0.305,tunemax=0.325,tuneint=1e-4):

    # 2D plot of a spectrum vs sigma and tune
    sigma=np.arange(sigmin,sigmax,sigmaint);
    tunes=np.arange(tunemin,tunemax,tuneint);
    w=np.zeros((len(sigma),len(tunes)));
    for isig,sig in enumerate(sigma):
    	for itune,tune in enumerate(tunes):
	    w[isig,itune]=abs(amp_complex_freq((sig,tune),x,y));

    fig,ax=init_figure();
    plot2D(w,tunemin,tunemax,sigmin,sigmax,"Tune"," $ \sigma $ ","2D spectrum",ax);
    end_figure(fig,ax);


def absfft0_minusexp(param,x,y):

    y2=y-(param[0]+1j*param[1])*np.exp(x*(param[2]+2.*np.pi*1j*param[3]));
    y3=fftshifted(-param[3],y2);
    return np.abs(y3[0])
   
 
def amp_complex_freq(param,x,y):

    # compute complex amplitude of an exponential term with a complex frequency
    rate=param[0]; # also called sigma in a few places
    tune=param[1];

    y2=y*np.exp((-rate-2.*np.pi*1j*tune)*x);
    amp=sum(y2)/float(len(x)); # Note: sum(y2)=fft.fft(y2)[0]
    
    return amp;


def quick_fit_exp(y):

    # fit by an exponential a signal. we take only points above the noise level.

    N=float(len(y));
    noiselevel=2.*noise(y);
    stdy=slide_std_deviation(y,2000);
    beg=np.where(stdy>noiselevel);beg=beg[0];

    if (len(beg)>0):
	p=np.polyfit(np.arange(beg[0],N),np.log(np.abs(y[beg[0]:])),1);
	growthrate=p[0];amp=np.exp(p[1]);
	#fig,ax=init_figure();ax.plot(np.abs(y),'.b');
	#ax.plot(stdy,'g',lw=3.);
	#ax.plot(np.arange(beg[0],N),amp*np.exp(growthrate*np.arange(beg[0],N)),'r',lw=3.);
    else:
    	growthrate=0.;amp=0.;
	
    return growthrate,amp;


def find_Laplace(y):

    # complex frequency spectral analysis (still experimental...)
    # see "Laplace spectrum exponential decomposition and pole-zero estimation", by M.J. Corinthios,
    # IEEE Proc.-Vis. Image Signal Process., Vol. 148, No. 5, October 2001
    
    # find the largest growing exponential term, and subtract it
    N=float(len(y));
    
    # first guess for sigma (highest growth rate)
    #sigma=optimize.fminbound(weightedlaplacespectrum_max,0,1.,args=(np.arange(N),y),xtol=1e-12,maxfun=2000);
    rate,amp1=quick_fit_exp(y)
        
    # first guess for shift to be applied before doing fft to get accurate frequency of highest term
    #alpha=optimize.fminbound(maxfftshifted,0,1./N,args=(y*np.exp(-sigma*np.arange(N)),0,N),xtol=1e-12);
    alpha=optimize.fminbound(maxfftshifted,0,1./N,args=(y*np.exp(-rate*np.arange(N)),0,N),xtol=1e-12);
    
    # 2D minimization to get more accurate sigma and alpha
    #x0=optimize.fmin(weightedlaplacespectrum2D_max,[sigma[0],alpha[0]],args=(np.arange(N),y),maxfun=20000,maxiter=200)
    #print "sigma[0]=",sigma[0],", alpha[0]=",alpha[0],", x0=",x0,", rate fitted=",rate;
    x0=(rate,alpha[0]);
    
    # frequency of highest term
    tune=float(np.argmax(np.abs(fftshifted(x0[1],y*np.exp(-x0[0]*np.arange(N))))))/N-x0[1];
    
    # complex amplitude of the exponential term found
    amp=amp_complex_freq((x0[0],tune),np.arange(N),y);
    
    # 4D optimization (WORSE WHEN USING IT)
    #x1=optimize.fmin(absfft0_minusexp,[amp.real,amp.imag,x0[0],tune],args=(np.arange(len(y)),y),maxfun=200000,maxiter=10000,xtol=1e-12,ftol=1e-12);
    
    # subtract
    y3=y-amp*np.exp(np.arange(N)*(x0[0]+2.*np.pi*1j*tune));
    #y3=y-(x1[0]+1j*x1[1])*np.exp(np.arange(N)*(x1[2]+2.*np.pi*1j*x1[3]));
    
    # return rise time (in turns), tune, array after subtraction, and amplitude
    return x0[0],tune,y3,amp
    

def Laplace(y,nharm,noiselevel=None,flagplot=False):

    # find harmonics with complex frequencies
    # we subtract the harmonics obtained at each step
    
    if (noiselevel==None):

	# find nharm harmonics
	y2=y;rates=[];tunes=[];amps=[];
	for i in range(nharm):
	    rate,tune,y3,amp=find_Laplace(y2);
	    y2=y3;
	    rates.append(rate);
	    tunes.append(tune);
	    amps.append(amp);

	    # to plot some weighted spectrum
	    #spec=[];nmin=-1000.;nmax=1000;
	    #for sigma in np.arange(nmin,nmax)/1000000.:
	    #    z2=np.exp(-sigma*np.arange(len(y3)))*y3;
	    #    z3=fft.fft(z2);
	    #    eps=sum(np.abs(z2)**2)/float(len(z2));
	    #    spec.append(np.abs(z3)**2/eps);
	    #pylab.figure();pylab.plot(np.arange(nmin,nmax)/1000000.,spec);pylab.show();
	    
    else:
    
	# find harmonics until we reach the noise level
	y2=y;rates=[];tunes=[];amps=[];
	k=0;col=['b','r','g','m','k','c','y'];
	if (flagplot): fig,ax=init_figure();

	# test plot of 2D spectrum
	#plot_weightedlaplacespectrum2D(np.arange(len(y)),y);
	#plot_spectrum2D(np.arange(len(y)),y);pylab.show();
	
	while (k<10)and((len(pylab.mlab.find(np.isnan(y2)))==0)and(np.max(np.abs(y2))>noiselevel)):
	    rate,tune,y3,amp=find_Laplace(y2);
	    if (flagplot): ax.plot(np.arange(len(y)),np.real(y2),'.'+col[k%7]);
	    k=k+1;y2=y3;
	    rates.append(rate);
	    tunes.append(tune);
	    amps.append(amp);
	
	    # test plot of 2D spectrum
	    #plot_weightedlaplacespectrum2D(np.arange(len(y)),y2,sigmin=0.,sigmax=0.02,sigmaint=1e-5);
    
	if (flagplot):
	    ax.plot(np.arange(len(y)),np.real(y2),'.'+col[k%7]);
	    ax.set_ylim((-np.max(np.abs(y)),np.max(np.abs(y))));
	    
	#pylab.show();
    
    return tunes,rates,amps;
    

def noise(y,npts=2000):

    # compute some noise level of a signal
    # here it is defined as the minimum along the signal of a "sliding standard deviation"
    # of the signal y, over npts
    
    std=slide_std_deviation(y,npts);
    # test plot
    #fig,ax=init_figure();ax.plot(np.arange(len(y)),np.real(y),'.b');
    #ax.plot(np.arange(len(y)),std,'r',lw=3.);pylab.show();
    
    return np.min(std);


def slide_std_deviation(data,period):
    # do a sliding standard deviation (over nb points=period)
    
    # might be very slow with large period 
    return np.array([np.sqrt(np.var(data[k:min(k+period,len(data))])) for k,j in enumerate(data[:-period])]);


def extract_lapl(tunes,rates,amps,Q,width):

    # sort and extract only tunes at width from the tune Q
    dic1=dict(np.transpose([tunes,amps]))
    dic2=dict(np.transpose([tunes,rates]))
    tunessort=np.sort(tunes)
    m=np.searchsorted(tunessort,Q-width)
    p=np.searchsorted(tunessort,Q+width,side='right')
    newtunes=tunessort[m:p]
    newamps=np.array([dic1[i] for i in newtunes]);
    newrates=np.array([dic2[i] for i in newtunes]);
   
    return newtunes,newrates,newamps;
    

def find_peak_lapl(tunes,rates,amps,Q,Qs,m):

    # find the headtail mode m from the spectrum of amplitude 'amp' vs. the tunes 'tunes'
    tu,ra,am=extract_lapl(tunes,rates,amps,Q+m*Qs,Qs/2.)
    if len(tu)==0:
        tuneshiftm=0;ratem=0;ampm=0;
    else:
	# tune shift of the peak from the zero intensity value Q+mQs
	tuneshiftm=tu[np.argmax(np.abs(am))]-(Q+m*Qs)
	# peak amplitude
	ampm=np.max(np.abs(am))
	# peak growth rate
	ratem=ra[np.argmax(np.abs(am))]

    return tuneshiftm,ratem,ampm
	

def plot_spec_comp(tunes,rates,amps,ax,tit=""):

    # plot complex frequency spectrum (marker is larger with higher growth rate)
    norm=np.max(rates);
    for i,t in enumerate(tunes):
	ax.plot(tunes[i],np.abs(amps[i]),'.',lw=2.5,ms=rates[i]*30./norm);
	
    ax.set_xlabel('Tune');
    ax.set_ylabel('Amplitude');
    ax.set_title(tit+'size of the spot proportional to growth rate');


if __name__ == "__main__":


    opt,args=parsse();
    
    #coefficent to apply to the noise level (we are converged when we are below this * noiselevel)
    coefnoise=2.;
    
    mmax=opt.MMAX[1]
    mmin=opt.MMAX[0]

    bunches=np.array(opt.BNUM);
    bunches.sort();
    print "Selected bunches for the spectrum plot:", bunches

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    # create list of associated legends (either with opt.LEG or 
    # the list of file names taking out all the common parameters
    # in the names)
    if (opt.LEG!=None):
        listleg=opt.LEG;
    else:
        listleg=takeout_common(listname);
	listleg=takeout_spaces(listleg);
    
    # find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
    gamma,circ,Trev,nbunch=extract_Headtail_param(listname[0][:-8])
    

    # initialize plots and variables for tune shifts and rise times
    #figtau,axtau=init_figure();
    #figtuneshift,axtuneshift=init_figure();

    taux=np.zeros((nbunch,-mmin+mmax+1));
    tauy=np.zeros((nbunch,-mmin+mmax+1));
    tuneshiftx=np.zeros((nbunch,-mmin+mmax+1));
    tuneshifty=np.zeros((nbunch,-mmin+mmax+1));
    ampx=np.zeros((nbunch,-mmin+mmax+1));
    ampy=np.zeros((nbunch,-mmin+mmax+1));
    if (opt.AVER):
    	taubeamx=np.zeros(-mmin+mmax+1);
	taubeamy=np.zeros(-mmin+mmax+1);
    	tuneshiftbeamx=np.zeros(-mmin+mmax+1);
	tuneshiftbeamy=np.zeros(-mmin+mmax+1);
    	ampbeamx=np.zeros(-mmin+mmax+1);
	ampbeamy=np.zeros(-mmin+mmax+1);
	    
    width=opt.WIDTH; # width chosen for the spectrum 


    # loop on filenames
    for ifile,filename in enumerate(listname):
    
    	print filename;
    
	# read from cfg file gamma, circumference, beta functions and tunes
	
	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    
	cfgfile=filename[:-8]+'.cfg';
	betax=float(read_cfg(cfgfile,"Horizontal_beta_function"));
	betay=float(read_cfg(cfgfile,"Vertical_beta_function"));  
	Qx=float(read_cfg(cfgfile,"Horizontal_tune"));
	Qy=float(read_cfg(cfgfile,"Vertical_tune"));  
	Qs=float(read_cfg(cfgfile,"Synchrotron_tune"));
	Qx=Qx-np.floor(Qx);Qy=Qy-np.floor(Qy); 
	print 'Betax=',betax,', Betay=',betay;
	print 'Qx=',Qx,', Qy=',Qy,', Qs=',Qs;

	#filex=open(filename[:-4]+'_Laplace_x'+opt.OUT,'w');
	#filey=open(filename[:-4]+'_Laplace_y'+opt.OUT,'w');
	#print >> filex, "Bunch\tHeadtail_mode_x\ttaux[s]\tTune_shiftx"
	#print >> filey, "Bunch\tHeadtail_mode_y\ttauy[s]\tTune_shifty"
	
	if (opt.AVER)and(opt.PLOT):
            # initialize spectrum plots
	    figspecx,axspecx=init_figure();
	    figspecy,axspecy=init_figure();
	    #figratex,axratex=init_figure();
	    #figratey,axratey=init_figure();

	file=[];
	for j,m in enumerate(range(mmin,mmax+1)):
	    file.append(open(filename[:-4]+'_Laplace_m'+str(m)+opt.OUT,'w'));
	    print >> file[j], "Bunch\ttaux[s]\ttauy[s]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	    
	filebunchmost=open(filename[:-4]+'_Laplace_most_'+opt.OUT,'w')
	filebunchmost1=open(filename[:-4]+'_Laplace_most_amp_'+opt.OUT,'w')
	filebunchmost2=open(filename[:-4]+'_Laplace_most_rate_'+opt.OUT,'w')
	print >> filebunchmost, "Bunch\tGrowthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	print >> filebunchmost1, "Bunch\tGrowthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	print >> filebunchmost2, "Bunch\tGrowthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"

	if (opt.AVER):
	    filebeam=open(filename[:-4]+'_aver_Laplace_'+opt.OUT,'w')
	    print >> filebeam, "Headtail_mode\ttaux[s]\ttauy[s]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	    filebeammost=open(filename[:-4]+'_aver_Laplace_most_'+opt.OUT,'w')
	    filebeammost1=open(filename[:-4]+'_aver_Laplace_most_amp_'+opt.OUT,'w')
	    filebeammost2=open(filename[:-4]+'_aver_Laplace_most_rate_'+opt.OUT,'w')
	    print >> filebeammost, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	    print >> filebeammost1, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
	    print >> filebeammost2, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"

    	# for the legend
        #fil=filename.replace("_"," ").replace("_prt.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read prt file
	x,xp,y,yp=read_prt_file(filename,nbunch);
	
	for bnum in range(nbunch,opt.NBMIN-1,-opt.EVERY):
	
	    print "bunch no ",bnum
	    
	    x1=x[bnum-1::nbunch];
	    xp1=xp[bnum-1::nbunch];
    	    y1=y[bnum-1::nbunch];
    	    yp1=yp[bnum-1::nbunch];

	    
	    if (bnum in bunches)and(opt.PLOT):
		# initialize spectrum plots
		figx,axx=init_figure();
		figy,axy=init_figure();

	    
	    if (opt.AVER):
	    	if (bnum==nbunch):
		    xave=x1;
		    yave=y1;
		    xpave=xp1;
		    ypave=yp1;
		else:
		    xave=x1+xave;
		    yave=y1+yave;	    
		    xpave=xp1+xpave;	    
		    ypave=yp1+ypave;	    


	    if (opt.STOP>=len(x1)):
	    	if (np.argmax(x1)>opt.BEG+100):
	    	    x2=x1[opt.BEG:np.argmax(x1)];
		    xp2=xp1[opt.BEG:np.argmax(x1)];
		else: 
	    	    x2=x1[opt.BEG:];
		    xp2=xp1[opt.BEG:];
		    
	    	if (np.argmax(y1)>opt.BEG+100):
		    y2=y1[opt.BEG:np.argmax(y1)];
	    	    yp2=yp1[opt.BEG:np.argmax(y1)];
		else:
	    	    y2=y1[opt.BEG:];
		    yp2=yp1[opt.BEG:];
		    
	    else:
		x2=x1[opt.BEG:opt.STOP];y2=y1[opt.BEG:opt.STOP]
		xp2=xp1[opt.BEG:opt.STOP];yp2=yp1[opt.BEG:opt.STOP]
	    
	    x3=x2-1j*betax*xp2;y3=y2-1j*betay*yp2; # valid only because betaprime=0 in Headtail
	    
	    # noise estimate of raw signals
	    noisex=coefnoise*noise(x3);
	    noisey=coefnoise*noise(y3);
	    
	    # complex frequency spectral analysis
	    tunesx,ratesx,ampsx=Laplace(x3,opt.NMOD[0],noiselevel=noisex);
	    tunesy,ratesy,ampsy=Laplace(y3,opt.NMOD[1],noiselevel=noisey);

	    # simply print the first tuneshift and rates extracted
	    if (len(ratesx)>0)and(len(ratesy)>0):
	    	print >> filebunchmost, bnum, "\t", ratesx[0]/Trev, "\t", ratesy[0]/Trev, "\t", tunesx[0]-Qx, "\t", tunesy[0]-Qy, "\t", abs(ampsx[0]), "\t", abs(ampsy[0]);

	    # find spectral lines associated with headtail modes
	    for j,m in enumerate(range(mmin,mmax+1)):
		tuneshiftx[bnum-1,j],ra,ampx[bnum-1,j]=find_peak_lapl(tunesx,ratesx,ampsx,Qx,Qs,m)
		if (ra!=0): taux[bnum-1,j]=1./ra;
		tuneshifty[bnum-1,j],ra,ampy[bnum-1,j]=find_peak_lapl(tunesy,ratesy,ampsy,Qy,Qs,m)
		if (ra!=0): tauy[bnum-1,j]=1./ra;
		print >> file[j], bnum, "\t", taux[bnum-1,j]*Trev, "\t", tauy[bnum-1,j]*Trev, "\t", tuneshiftx[bnum-1,j]-Qx, "\t", tuneshifty[bnum-1,j]-Qy, "\t", abs(ampx[bnum-1,j]), "\t", abs(ampy[bnum-1,j]);
		
	    if (bnum in bunches)and(opt.PLOT):
		# plot spectra
		plot_spec_comp(tunesx,ratesx,ampsx,axx,tit='Horizontal, bunch '+str(bnum));
		plot_spec_comp(tunesy,ratesy,ampsy,axy,tit='Vertical, bunch '+str(bnum));
	    
	    # find spectral line of highest amplitude (within width/2 from the tune)
	    tumostx,ramostx,ampmostx=find_peak_lapl(tunesx,ratesx,ampsx,Qx,width,0)
	    tumosty,ramosty,ampmosty=find_peak_lapl(tunesy,ratesy,ampsy,Qy,width,0)
	    print >> filebunchmost1, bnum, "\t", ramostx/Trev, "\t", ramosty/Trev, "\t", tumostx-Qx, "\t", tumosty-Qy, "\t", abs(ampmostx), "\t", abs(ampmosty);
	    # find spectral line of highest growth rate (within width/2 from the tune)
	    # (to do this we just need to invert 'amps' and 'rates' in find_peak_lapl)
	    tumostx,ampmostx,ramostx=find_peak_lapl(tunesx,ampsx,ratesx,Qx,width,0)
	    tumosty,ampmosty,ramosty=find_peak_lapl(tunesy,ampsy,ratesy,Qy,width,0)
	    print >> filebunchmost2, bnum, "\t", ramostx/Trev, "\t", ramosty/Trev, "\t", tumostx-Qx, "\t", tumosty-Qy, "\t", abs(ampmostx), "\t", abs(ampmosty);
	
    	for j,m in enumerate(range(mmin,mmax+1)): file[j].close();
	#filex.close();filey.close();
	
	filebunchmost1.close();
	filebunchmost2.close();
			    
	pat=['-','--','.','x','o','+','d','s','<','>','t']
	#for j,m in enumerate(range(mmin,mmax+1)):
	    #axtau.plot(range(nbunch,0,-opt.EVERY),taux[::-opt.EVERY,j]*Trev,pat[j%8],label='mode m='+str(m)+', '+fil+' '+r'$\tau_x$',lw=2.5);
	    #axtau.plot(range(nbunch,0,-opt.EVERY),tauy[::-opt.EVERY,j]*Trev,pat[j%8],label='mode m='+str(m)+', '+fil+' '+r'$\tau_y$',lw=2.5);
	    #axtuneshift.plot(range(nbunch,0,-opt.EVERY),tuneshiftx[::-opt.EVERY,j],pat[j%8],label='mode m='+str(m)+', '+fil+' '+', x',lw=2.5);
	    #axtuneshift.plot(range(nbunch,0,-opt.EVERY),tuneshifty[::-opt.EVERY,j],pat[j%8],label='mode m='+str(m)+', '+fil+' '+', y',lw=2.5);


	if (opt.AVER):
	    
	    print "Beam average"
	    
	    # compute average first
	    for bnum in range(nbunch,0,-1):
	
		x1=x[bnum-1::nbunch];
		xp1=xp[bnum-1::nbunch];
    		y1=y[bnum-1::nbunch];
    		yp1=yp[bnum-1::nbunch];

	    	if (bnum==nbunch):
		    xave=x1;
		    yave=y1;
		    xpave=xp1;
		    ypave=yp1;
		else:
		    xave=x1+xave;
		    yave=y1+yave;	    
		    xpave=xp1+xpave;	    
		    ypave=yp1+ypave;	    

	    xave/=float(nbunch);
	    yave/=float(nbunch);
	    xpave/=float(nbunch);
	    ypave/=float(nbunch);
	    

	    if (opt.STOP>=len(xave)):
	    	xave2=xave[opt.BEG::];yave2=yave[opt.BEG::];
	    	xpave2=xpave[opt.BEG::];ypave2=ypave[opt.BEG::];
	    else:
		xave2=xave[opt.BEG:opt.STOP];yave2=yave[opt.BEG:opt.STOP]
		xpave2=xpave[opt.BEG:opt.STOP];ypave2=ypave[opt.BEG:opt.STOP]
	    
	    xave3=xave2-1j*betax*xpave2;
	    yave3=yave2-1j*betay*ypave2;
	    
	    # noise estimate of average raw signals
	    noisex=coefnoise*noise(xave3);
	    noisey=coefnoise*noise(yave3);
	    
	    # complex frequency spectral analysis
	    tunesx,ratesx,ampsx=Laplace(xave3,opt.NMOD[0],noiselevel=noisex,flagplot=opt.PLOT);
	    tunesy,ratesy,ampsy=Laplace(yave3,opt.NMOD[1],noiselevel=noisey,flagplot=opt.PLOT);

	    # simply print tuneshifts and rates in their natural order
	    for i in range(min(len(tunesx),len(tunesy))):
		print >> filebeammost, ratesx[i]/Trev, "\t", ratesy[i]/Trev, "\t", tunesx[i]-Qx, "\t", tunesy[i]-Qy, "\t", abs(ampsx[i]), "\t", abs(ampsy[i]);

	    # find spectral lines associated with headtail modes
	    for j,m in enumerate(range(mmin,mmax+1)):
		tuneshiftbeamx[j],ra,ampbeamx[j]=find_peak_lapl(tunesx,ratesx,ampsx,Qx,Qs,m)
		if (ra!=0): taubeamx[j]=1./ra;
		tuneshiftbeamy[j],ra,ampbeamy[j]=find_peak_lapl(tunesy,ratesy,ampsy,Qy,Qs,m)
		if (ra!=0): taubeamy[j]=1./ra;
		
		#print fil, ", mode m=",m," - Average rise times of the beam (x & y) in sec.: ", taubeamx[j]*Trev, taubeamy[j]*Trev;
		#print fil, ", mode m=",m," - Average tune shifts of the beam (x & y): ", tuneshiftbeamx[j], tuneshiftbeamy[j];
		print >> filebeam, m, "\t", taubeamx[j]*Trev, "\t", taubeamy[j]*Trev,"\t", tuneshiftbeamx[j]-Qx, "\t", tuneshiftbeamy[j]-Qy, "\t", abs(ampbeamx[j]), "\t", abs(ampbeamy[j]);
		
	    # find spectral lines of highest amplitudes (within width/2 from the tune)
	    indx=np.flipud(np.argsort(np.abs(ampsx)));indy=np.flipud(np.argsort(np.abs(ampsy)));
	    #print tunesx,ratesx,ampsx,indx;
	    #print tunesy,ratesy,ampsy,indy;
	    for i in range(min(len(indx),len(indy))):
		print >> filebeammost1, ratesx[indx[i]]/Trev, "\t", ratesy[indy[i]]/Trev, "\t", tunesx[indx[i]]-Qx, "\t", tunesy[indy[i]]-Qy, "\t", abs(ampsx[indx[i]]), "\t", abs(ampsy[indy[i]]);
	    
	    # find spectral line of highest growth rate (within width/2 from the tune)
	    indx=np.argsort(ratesx);indy=np.argsort(ratesy);
	    for i in range(min(len(indx),len(indy))):
		print >> filebeammost2, ratesx[indx[i]]/Trev, "\t", ratesy[indy[i]]/Trev, "\t", tunesx[indx[i]]-Qx, "\t", tunesy[indy[i]]-Qy, "\t", abs(ampsx[indx[i]]), "\t", abs(ampsy[indy[i]]);

	#    indx=np.where(taubeamx>0);
	#    if (len(indx[0])>0):
	#	rax=1./np.min(taubeamx[indx[0]]);ix=indx[0][np.argmin(taubeamx[indx[0]])];
	#    else:
	#	rax=0.;ix=0;
	#    indy=np.where(taubeamy>0);
	#    if (len(indy[0])>0):
	#	ray=1./np.min(taubeamy[indy[0]]);iy=indy[0][np.argmin(taubeamy[indy[0]])];
	#    else:
	#	ray=0.;iy=0;
	#    print >> filebeammost, rax, "\t", ray,"\t", tuneshiftbeamx[ix], "\t", tuneshiftbeamy[iy], "\t", abs(ampbeamx[ix]), "\t", abs(ampbeamy[iy]);
	    
	    if (opt.PLOT):
		# plot spectra
		plot_spec_comp(tunesx,ratesx,ampsx,axspecx,tit='Horizontal, average');
		plot_spec_comp(tunesy,ratesy,ampsy,axspecy,tit='Vertical, average');
	    #pat=['.','x','o','+','d','s','<','>','^','v','p','*','h']
	    #tunesx=np.array(tunesx);ampsx=np.array(ampsx);
	    #tunesy=np.array(tunesy);ampsy=np.array(ampsy);
	    #axspecx.plot(np.sort(tunesx),np.abs(ampsx[np.argsort(tunesx)]),'-',lw=2.5);
	    #axspecy.plot(np.sort(tunesy),np.abs(ampsy[np.argsort(tunesy)]),'-',lw=2.5);
	    #axratex.plot(tunex[bnum-1,:],1./taux[bnum-1,:],pat[bnum/opt.EVERY%len(pat)],label='bunch '+str(bnum),lw=2.5);
	    #axratey.plot(tuney[bnum-1,:],1./tauy[bnum-1,:],pat[bnum/opt.EVERY%len(pat)],label='bunch '+str(bnum),lw=2.5);

   	    filebeam.close();
	    filebeammost.close();
	    filebeammost1.close();
	    filebeammost2.close();
			


    #axtau.set_xlabel("Bunch number");
    #axtau.set_ylabel("Rise time [s]");
    #axtau.legend(loc=0);
    #axtuneshift.set_xlabel("Bunch number");
    #axtuneshift.set_ylabel("Tune shift");
    #axtuneshift.legend(loc=0);
	

    if (opt.PLOT): pylab.show();
    sys.exit()
