#!/usr/bin/python2.6

import sys,os
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
from SussixNM import *
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,re,dateutil,random,pytz
import numpy as np
import scipy.io as sio
from string import split, replace
from optparse import OptionParser
from Timber import parseout
import math
import matplotlib
from plot_lib import plot,set_fontsize,init_figure,end_figure
from io_lib import list_files,test_and_create_dir
from tables_lib import select_and_average
from read_Headtail_prt_sussix import gettune,collect_and_fit,extract,find_peak
from read_Headtail_prt_fit import fit


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the .mat name of the file. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      #help="Specify the TIMBER .csv name of the file (several -f options possible -> several files); it can also be regular expressions encompassing several files at a time (with e.g. '*', etc.)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-d", "--damper",action="store_false",
                      help="Specify if we do not distinguish openloop damper to closedloop damper. If we distinguish, separate plots and analysis are done for openloop/closedloop cases. Default=distinguish the two.",
                      metavar="DAMPER", default=True,dest="DAMPER")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify the kinetic energy (GeV). Default=1.4 GeV",
                      metavar="ENERGY", default=1.4,dest="ENERGY")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of ms (from the injection) on which we compute the intensity average and after which we fit the amplitude of the tune line, for -q option (default=0).",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-i", "--intensity",action="store_true",
                      help="Specify if we plot and save the intensity data. Default=no intensity plot.",
                      metavar="INT", default=False,dest="INT")
    parser.add_option("-o", "--output",help="Specify output filenames additional suffix. Default=no suffix.",
                      default='',dest="OUT")
    parser.add_option("-r", "--raw",action="store_true",
                      help="Specify if we plot and save raw Qmeter data (from the BQSB class), and do a tune and rate analysis. Default=no raw data plot.",
                      metavar="RAW", default=False,dest="RAW")
    parser.add_option("-s", "--scope",action="store_true",
                      help="Specify if we plot and save scope traces (from the SCOPE31 class). Default=no scope data plot.",
                      metavar="SCOPE", default=False,dest="SCOPE")
    parser.add_option("-u", "--boundmin",type=float,
                      help="Specify minimum amplitude of spectral lines to begin the fit, for -q option (default=1e-7)",
                      metavar="BMIN", default=1.e7,dest="BMIN")
    parser.add_option("--scopepar",type=int,nargs=2,
                      help="Specify scope parameter: number of ns per division (10 divisions on the whole). Default=10.",
                      metavar="SCOPEPAR", default=10,dest="SCOPEPAR")
    parser.add_option("--tunes",type=float,nargs=3,
                      help="Specify the unperturbed tunes x & y (fractional part) and the width of the spectrum around them (for analysis with -q option). Default=(0.1,0.45,0.05).",
                      metavar="TUNES", default=(0.1,0.45,0.02),dest="TUNES")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


if __name__ == "__main__":


    opt,args=parsse();
    
    # revolution frequency and period (assumes protons)
    print 'Kinetic energy: ',opt.ENERGY,' GeV';
    E0=0.938272;
    energy_total=opt.ENERGY+E0
    circumference=200*np.pi; # PS circumference in m
    gamma=energy_total/E0 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./circumference;
    Trev=1./frev; # revolution period in s
    print 'beta=',beta,', gamma=',gamma,', Trev=',Trev*1e6,'us';

    Qx=opt.TUNES[0]; # fractional part of horizontal tune
    Qy=opt.TUNES[1]; # fractional part of vertical tune
    width=opt.TUNES[2]; # width chosen for the spectrum
    injection=170; # injection time in ms (for intensity)

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if (opt.INT):
        # initialize list with all intensity data (over all files) along the cycle
	intensity=[];

    if (opt.SCOPE):
        # initialize list with all scope data (over all files)
	scope=[[],[],[]];

    if (opt.RAW):

        # initialize lists with all raw Qmeter data (over all files) along the cycle
	raw=[[],[]];

    	# initialize arrays (tunes, their std deviation, growth rates and intensities for each file)
	tunesx=np.zeros(len(listname));
	sigtunesx=np.zeros(len(listname));
	ratesx=np.zeros(len(listname));
	ratesfitx=np.zeros(len(listname));
	tunesy=np.zeros(len(listname));
	sigtunesy=np.zeros(len(listname));
	ratesy=np.zeros(len(listname));
	ratesfity=np.zeros(len(listname));
	inten=np.zeros(len(listname));
	
	# initialize damper openloop array
	openloop=np.zeros(len(listname));
	# initialize lists with tunes (over all files) along the cycle
	tunebx=[];tuneby=[];
    
    for ifile,filename in enumerate(listname):
    
	dirname=filename[:-4];
	# create (if not already there) directory to put plots
	test_and_create_dir(dirname);
	print "Analysing ",filename;
    
	# read MATLAB kind of file and extract data structure
	mat=sio.loadmat(filename,struct_as_record=False,squeeze_me=True);
	data=mat['data'];
	
	loopstr='';
	
	if (opt.DAMPER):
	    # openloop = 0 if damper is in closed loop (on), 1 if it is in open loop (off)
	    loop=['_closedloop','_openloop'];
	    
	    openloop[ifile]=data.PA.TFB_DSPU_H.SwitchesAcq.openLoop;
	    if (openloop[ifile]!=0)and(openloop[ifile]!=1): print "    Pb with openloop: ",openloop[ifile];openloop[ifile]=nan;
	    else:
		print "    Openloop=",openloop[ifile];
		loopstr=loop[int(openloop[ifile])];


	if (opt.INT):
	
	    if (not(isscalar(data.PR.STRTR.Samples.samples)))and(len(data.PR.STRTR.Samples.samples)>0):
		# extract and plot intensity along the cycle
		intensity.append(data.PR.STRTR.Samples.samples);
		timesI=data.PR.STRTR.Samples.firstSampleTime+np.arange(len(intensity[-1])); # times in ms
		if ('oldtimesI' in locals()):
		    if (timesI[0]!=oldtimesI[0])or(len(timesI)!=len(oldtimesI)): print "Pb: not the same intensity times in all files !", timesI,oldTimesI;sys.exit();
		oldtimesI=timesI;

		fig,ax=init_figure();
		plot(timesI,intensity[-1],'','-b','Intensity [ $ 10^{10} $ p+/b]',ax,0,lw=3.,xlab='Time along the cycle [ms]');
		end_figure(fig,ax,save=dirname+'/'+dirname+'_intensity'+loopstr+opt.OUT);

	    else:
	    	print "    No intensity data";
		intensity.append(nan);
	    
	
	if (opt.SCOPE):
	
	    if ((not(isscalar(data.PR.SCOPE31.CH01.Acquisition.value)))and(not(isscalar(data.PR.SCOPE31.CH02.Acquisition.value))))and(not(isscalar(data.PR.SCOPE31.CH03.Acquisition.value))):

		# extract and plot traces from the scope
		name=['sum','x','y'];
		# number of points and number of traces
		N=data.PX.N_MEAS.Delay.delay;
		if ('oldN' in locals()):
		    if (N==-1): N=oldN;
		    if(N!=oldN): print "Pb: not the same nb of points in scope data in all files !",N,oldN;sys.exit();
		if (N==-1): N=200; # value by default
		oldN=N;
		ntraces=len(data.PR.SCOPE31.CH01.Acquisition.value)/N;

		# position in the cycle in ms at which traces are taken
		delay=data.PX.C_STRT_OAS.Delay.delay;
		# traces taken every 'inter' turns
		inter=data.PX.TREV_INTER.Delay.delay
		# sum signal
		scope[0].append(reshape(data.PR.SCOPE31.CH01.Acquisition.value,(ntraces,N)).transpose());
		# x signal
		scope[1].append(reshape(data.PR.SCOPE31.CH02.Acquisition.value,(ntraces,N)).transpose());
		# y signal
		scope[2].append(reshape(data.PR.SCOPE31.CH03.Acquisition.value,(ntraces,N)).transpose());

		for i in range(3):
	    	    fig,ax=init_figure();
		    # always take out the third trace (synchronization pb)
		    ind=concatenate((np.arange(1),np.arange(2,ntraces)));
		    x=np.arange(float(N))*10*opt.SCOPEPAR/float(N);
		    # always assuming 10 divisions on the whole
		    plot(x,scope[i][-1][:,ind],'','-b','Scope '+name[i]+' signal [a.u.]',ax,0,lw=2,xlab='Longitudinal position in the bunch [ns]');
		    # plot in red one single trace
		    plot(x,scope[i][-1][:,ntraces/2],'One single trace','-r','Scope '+name[i]+' signal [a.u.]',ax,0,lw=3,xlab='Longitudinal position in the bunch [ns]');
		    # plot in light grey a sine wave of frequency 20 MHz=0.02 GHz
		    ampsine=ceil(np.max(scope[i][-1])/5000.)*5000;
		    plot(x,ampsine*(1+0.5*np.sin(2.*np.pi*0.02*x)),'Sine of freq. 20 MHz','-','Scope '+name[i]+' signal [a.u.]',ax,0,lw=3,xlab='Longitudinal position in the bunch [ns]',colr=[0.8,0.8,0.8]);
		    ax.set_title(str(ntraces)+' traces, every '+str(inter)+' turn(s), starting at '+str(delay)+' ms in the cycle');
		    end_figure(fig,ax,save=dirname+'/'+dirname+'_scope_'+name[i]+'_'+str(delay)+'ms'+loopstr+opt.OUT);

	    else:
	        print "    No scope data";
		for i in range(3): scope[i].append(nan);


	if (opt.RAW):
	
	    # intensity (average over opt.BEG first ms)
	    intensitytmp=data.PR.STRTR.Samples.samples;
	    if (not(isscalar(intensitytmp)))and(len(intensitytmp)>0):
		timesItmp=data.PR.STRTR.Samples.firstSampleTime+np.arange(len(intensitytmp)); # times in ms
		ind=pylab.mlab.find((timesItmp>injection)*(timesItmp<injection+opt.BEG));
		inten[ifile]=np.average(intensitytmp[ind]);
	    else: inten[ifile]=nan;print "  warning: no intensity data.";
		
	    # extract raw Qmeter data
	    name=['Hor','Ver'];
	    # horizontal
	    raw[0].append(data.PR.BQSB.Acquisition.rawDataH);
	    # vertical
	    raw[1].append(data.PR.BQSB.Acquisition.rawDataV);
	    
	    if ((not(isscalar(raw[0][-1])))and(data.PR.BQSB.SamplerAcquisition.nbOfMeas!=-1))and(not(isscalar(data.PR.BQSB.SamplerAcquisition.measStamp))):
	    
		# parameters
		nacq=data.PR.BQSB.SamplerAcquisition.nbOfMeas; # nb of acquisitions
		nturns=len(raw[0][-1])/nacq; # nb of turns in each acquisition
		# times (ms) at the end of each acquisition
		timesAcq=np.array(data.PR.BQSB.SamplerAcquisition.measStamp[:nacq],dtype=float);
		if ('oldnturns' in locals()):
		    if (nturns!=oldnturns)or(len(timesAcq)!=len(oldtimesAcq)): print "Pb: not the same nb of points and/or times in raw data in all files !",nturns,oldnturns,timesAcq,oldtimesAcq;sys.exit();
		    elif (np.abs(timesAcq-oldtimesAcq)>1).any(): print "Pb: not the same times in raw data in all files !",np.abs(timesAcq-oldtimesAcq);sys.exit();
		oldnturns=nturns;oldtimesAcq=timesAcq;
		
		# reconstruct times for the complete raw data (in ms)
		timesraw=[];
		for i in range(nacq):
	            timesraw.extend([timesAcq[i]-Trev*j*1000 for j in range(nturns-1,-1,-1)])
		timesraw=np.array(timesraw);

		# injection time (estimate)
		injtime=timesAcq[0]-nturns*Trev*1000;

	        # plot raw Qmeter data
		figraw=[];axraw=[];tau=[];
		for i in range(2):
	    	    fig,ax=init_figure();figraw.append(fig);axraw.append(ax);
		    plot(timesraw,raw[i][-1],'Q meter','.b','Qmeter raw '+name[i]+'. data [a.u.]',axraw[i],0,lw=2.5,xlab='Time along the cycle [ms]',plotevery=20);
		    
		    # do an exponential fit of the envelop
		    tau=fit(timesraw,raw[i][-1],axraw[0],True,'Q meter','Qmeter raw '+name[i]+'. data [a.u.]',0.001,[1e8,1.e9],0,(60,10),opt.BEG+injtime,0,col='r');
		    axraw[i].set_xlabel('Time along the cycle [ms]');
		    end_figure(figraw[i],axraw[i],save=dirname+'/'+dirname+'_Qmeter_raw_'+name[i]+loopstr+opt.OUT);
		    
		    if (i==0): ratesfitx[ifile]=1./tau;
		    else: ratesfity[ifile]=1./tau;

	    	# tune and growth rate analysis

		#initialize figures for tune shift and amplitude vs turns plots
    		figampx,axampx=init_figure();
    		figampy,axampy=init_figure();
    		figtux,axtux=init_figure();
    		figtuy,axtuy=init_figure();

		tunebx.append([]);ampbx=[];tuneby.append([]);ampby=[];

		for it,time in enumerate(timesAcq):

		    suss=gettune(raw[0][-1][it*nturns:(it+1)*nturns],np.zeros(nturns),raw[1][-1][it*nturns:(it+1)*nturns],np.zeros(nturns),1.,1.,Qx,Qy,0.01,narm=300);

		    # sort and extract only tunes at "width" from nominal tunes
		    tunetmpx,amptmpx=extract(suss.ox,suss.ax,Qx,width)
		    tunetmpy,amptmpy=extract(suss.oy,suss.ay,Qy,width)

		    # find spectral line of highest amplitude
		    tu,am=find_peak(tunetmpx,amptmpx,Qx,width,0)
		    tunebx[-1].append(tu+Qx);ampbx.append(am);
		    tu,am=find_peak(tunetmpy,amptmpy,Qy,width,0)
		    tuneby[-1].append(tu+Qy);ampby.append(am);

		# average acquisition times (over nturns)
		newtimesAcq=timesAcq-nturns*Trev*1000/2.;

		# find first time for fit
		first=opt.BEG+injtime;
		ifirst=np.where(newtimesAcq>first);ifirst=ifirst[0];
		if (len(ifirst)>0): ifirst=ifirst[0];
		else: ifirst=0;

		# last turns for fit (at the maximum of the amplitude)
		lastx=newtimesAcq[np.argmax(ampbx[ifirst:])+ifirst];
		lasty=newtimesAcq[np.argmax(ampby[ifirst:])+ifirst];
		#print ifirst,first,lastx,lasty;

		# collect tunes and amplitudes of most unstable mode, and fit
		tunesx[ifile],taux,sigtunesx[ifile]=collect_and_fit(tunebx[-1],ampbx,
			newtimesAcq,'SUSSIX, horizontal','x',opt.BMIN,True,axtux,axampx,
			0.001,'Tune','Amplitude of the tune line',firstturn=first,
			lastturn=lastx,flagaver=True,xlab="Time along the cycle [ms]")
		tunesy[ifile],tauy,sigtunesy[ifile]=collect_and_fit(tuneby[-1],ampby,
			newtimesAcq,'SUSSIX, vertical','y',opt.BMIN,True,axtuy,axampy,
			0.001,'Tune','Amplitude of the tune line',firstturn=first,
			lastturn=lasty,flagaver=True,xlab="Time along the cycle [ms]")
		ratesx[ifile]=1./(taux*0.001);
		ratesy[ifile]=1./(tauy*0.001);		

		# finalization
	        end_figure(figtux,axtux,save=dirname+'/'+dirname+'_tunesx'+loopstr+opt.OUT);
	        end_figure(figtuy,axtuy,save=dirname+'/'+dirname+'_tunesy'+loopstr+opt.OUT);
	        end_figure(figampx,axampx,save=dirname+'/'+dirname+'_ampsx'+loopstr+opt.OUT);
	        end_figure(figampy,axampy,save=dirname+'/'+dirname+'_ampsy'+loopstr+opt.OUT);


	    else:
	        print "    No raw Qmeter data";
		raw[0][-1]=nan;raw[1][-1]=nan;
		tunesx[ifile]=nan;tunesy[ifile]=nan;
		ratesx[ifile]=nan;ratesy[ifile]=nan;
		ratesfitx[ifile]=nan;ratesfity[ifile]=nan;
		sigtunesx[ifile]=nan;sigtunesy[ifile]=nan;
		tunebx.append(nan);tuneby.append(nan);
		
	    # output
	    filetunes=open(dirname+'/'+dirname+'_Qmeter_tunes'+loopstr+opt.OUT+'.txt','w')
	    print >> filetunes, "Intensity[10^10]\tGrowthratex[s-1]\tGrowthratey[s-1]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]\tTune_x\tTune_y\tSigma_tune_x\tSigma_tune_y"
	    print >> filetunes, inten[ifile], "\t", ratesx[ifile], "\t", ratesy[ifile],"\t", ratesfitx[ifile], "\t", ratesfity[ifile],"\t", tunesx[ifile], "\t", tunesy[ifile],"\t", sigtunesx[ifile], "\t", sigtunesy[ifile];
	    filetunes.close();
		

    # for final averages
    if (opt.DAMPER): strdamper=['_closedloop','_openloop'];legdamper=['Damper in closed loop','Damper in open loop'];
    else: strdamper=[''];legdamper=[''];

    # define colors for average plots (closed loop = first color, openloop = second color);
    col=['-r','-b'];
    
    if (opt.INT):
	# compute and plot the average intensity along the cycle
	
	fig,ax=init_figure();
	
	for iloop in range(1,-1,-1):

	    # openloop is an array with 0 if damper is in closed loop for the file, 1 if it is in open loop so
	    # - when iloop=0: take only files where loop is closed, i.e. openloop=0
	    # - when iloop=1: take only files where loop is open, i.e. openloop=1

	    intensityav,a,b=select_and_average(intensity,openloop,iloop);
	    
	    if not(isnan(intensityav).all()): plot(timesI,intensityav,legdamper[iloop],col[iloop],'Average intensity [ $ 10^{10} $ p+/b]',ax,0,lw=3.,xlab='Time along the cycle [ms]');
	
	end_figure(fig,ax,save='average_intensity'+opt.OUT);

    
    if (opt.SCOPE):
	# compute and plot the average scope signal
	name=['sum','x','y'];
	
	for i in range(3):

	    fig,ax=init_figure();
	    
	    for iloop in range(1,-1,-1):

		scopeav,a,b=select_and_average(scope[i],openloop,iloop);
		if not(isnan(scopeav).all()):
		    N=len(scopeav[:,0]);
		    ind=concatenate((np.arange(1),np.arange(2,ntraces)));
		    plot(np.arange(float(N))*10*opt.SCOPEPAR/float(N),scopeav[:,ind[:-1]],'',col[iloop],'Average scope '+name[i]+' signal [a.u.]',ax,0,lw=2,xlab='Longitudinal position in the bunch [ns]');
		    # last trace with legend
		    plot(np.arange(float(N))*10*opt.SCOPEPAR/float(N),scopeav[:,ind[-1]],legdamper[iloop],col[iloop],'Average scope '+name[i]+' signal [a.u.]',ax,0,lw=2,xlab='Longitudinal position in the bunch [ns]');
		
	    ax.set_title(str(ntraces)+' traces, every '+str(inter)+' turn(s), starting at '+str(delay)+' ms in the cycle');
	    end_figure(fig,ax,save='average_scope_'+name[i]+'_'+str(delay)+'ms'+opt.OUT);

    
    if (opt.RAW):
	# compute and plot the average raw data signal
    	
	name=['Hor','Ver'];
	
	for i in range(2):
	    
	    figav,axav=init_figure();
	    figmax,axmax=init_figure();
	
	    for iloop in range(1,-1,-1):

		rawav,rawmax,b=select_and_average(raw[i],openloop,iloop,flagmax=True);
		if not(isnan(rawav).all()):
		    plot(timesraw,rawav,legdamper[iloop],col[iloop],'Average Qmeter raw '+name[i]+'. data [a.u.]',axav,0,lw=2.5,xlab='Time along the cycle [ms]',plotevery=20);
		    plot(timesraw,rawmax,legdamper[iloop],col[iloop],'Max. Qmeter raw '+name[i]+'. data [a.u.]',axmax,0,lw=2.5,xlab='Time along the cycle [ms]',plotevery=20);

	    end_figure(figav,axav,save='average_Qmeter_raw_'+name[i]+opt.OUT);
	    end_figure(figmax,axmax,save='max_Qmeter_raw_'+name[i]+opt.OUT);

   
 	# compute the average and sigmas over all files of the intensity, tunes and growth rates

	figx,axx=init_figure();
	figy,axy=init_figure();

	for iloop in range(1,-1,-1):

	    ind=pylab.mlab.find(openloop==iloop);
	    inten1=inten[ind];
	    ratesx1=ratesx[ind];
	    ratesy1=ratesy[ind];
	    ratesfitx1=ratesfitx[ind];
	    ratesfity1=ratesfity[ind];
	    tunesx1=tunesx[ind];
	    tunesy1=tunesy[ind];
	    sigtunesx1=sigtunesx[ind];
	    sigtunesy1=sigtunesy[ind];

	    # take out nan and infinities
	    inten1=inten1[pylab.mlab.find(isfinite(inten1))];
	    ratesx1=ratesx1[pylab.mlab.find(isfinite(ratesx1))];
	    ratesy1=ratesy1[pylab.mlab.find(isfinite(ratesy1))];
	    ratesfitx1=ratesfitx1[pylab.mlab.find(isfinite(ratesfitx1))];
	    ratesfity1=ratesfity1[pylab.mlab.find(isfinite(ratesfity1))];
	    sigtunesx1=sigtunesx1[pylab.mlab.find(isfinite(tunesx1))];
	    sigtunesy1=sigtunesy1[pylab.mlab.find(isfinite(tunesy1))];
	    tunesx1=tunesx1[pylab.mlab.find(isfinite(tunesx1))];
	    tunesy1=tunesy1[pylab.mlab.find(isfinite(tunesy1))];

	    # compute averages and standard deviation
	    inten_av=np.average(inten1);inten_sig=np.sqrt(np.var(inten1));
	    ratesx_av=np.average(ratesx1);ratesx_sig=np.sqrt(np.var(ratesx1));
	    ratesy_av=np.average(ratesy1);ratesy_sig=np.sqrt(np.var(ratesy1));
	    ratesfitx_av=np.average(ratesfitx1);ratesfitx_sig=np.sqrt(np.var(ratesfitx1));
	    ratesfity_av=np.average(ratesfity1);ratesfity_sig=np.sqrt(np.var(ratesfity1));
	    tunesx_av=np.average(tunesx1);
	    tunesy_av=np.average(tunesy1);
	    tunesx_sig=np.sqrt((np.sum(sigtunesx1**2)+np.sum(tunesx1**2))/float(len(tunesx1))-tunesx_av**2);
	    tunesy_sig=np.sqrt((np.sum(sigtunesy1**2)+np.sum(tunesy1**2))/float(len(tunesy1))-tunesy_av**2);

	    # output
	    filefinal=open('final_tunes_and_growthrates'+strdamper[iloop]+opt.OUT+'.txt','w')
	    print >> filefinal, "Intensity[10^10]\tSigma_intensity[10^10]\tGrowthratex[s-1]\tGrowthratey[s-1]\tSigma_rate_x[s-1]\tSigma_rate_y[s-1]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]\tSigma_rate_fitx[s-1]\tSigma_rate_fity[s-1]\tTune_x\tTune_y\tSigma_tune_x\tSigma_tune_y"
	    print >> filefinal, inten_av, "\t", inten_sig, "\t", ratesx_av, "\t", ratesy_av,"\t", ratesx_sig, "\t", ratesy_sig,"\t", ratesfitx_av, "\t", ratesfity_av,"\t", ratesfitx_sig, "\t", ratesfity_sig,"\t", tunesx_av, "\t", tunesy_av,"\t", tunesx_sig, "\t", tunesy_sig;
	    filefinal.close();

	    # now compute and plot average tunes along the cycle
	    tunexav,a,tunexvar=select_and_average(tunebx,openloop,iloop,flagvar=True);
	    tuneyav,a,tuneyvar=select_and_average(tuneby,openloop,iloop,flagvar=True);

	    if not(isnan(tunexav).all()): axx.errorbar(newtimesAcq,tunexav,yerr=tunexvar,fmt=col[iloop],label=legdamper[iloop],lw=3);
	    if not(isnan(tuneyav).all()): axy.errorbar(newtimesAcq,tuneyav,yerr=tuneyvar,fmt=col[iloop],label=legdamper[iloop],lw=3);


	axx.set_xlabel('Time along the cycle [ms]');
	axx.set_ylabel('Average hor. tune');
	axy.set_xlabel('Time along the cycle [ms]');
	axy.set_ylabel('Average ver. tune');
	end_figure(figx,axx,save='average_tunesx'+opt.OUT);
	end_figure(figy,axy,save='average_tunesy'+opt.OUT);
	    
