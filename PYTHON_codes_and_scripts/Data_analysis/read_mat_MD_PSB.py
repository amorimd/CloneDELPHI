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
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of ms (from the injection) on which we compute the intensity average and after which we fit the amplitude of the tune line, for -q option (default=0).",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-i", "--intensitygammafrev",action="store_true",
                      help="Specify if we plot and save the intensity, gamma and frev data. Default=no such plots.",
                      metavar="INT", default=False,dest="INT")
    parser.add_option("-n", "--nacq",type=int,
                      help="Specify the number of turns in each raw data acquisitions. Default is given by data.BR.BQSB.SamplerAcquisition.nbOfMeas, or, if not present, the highest power of 2 dividing the length of the raw data.",
                      metavar="NACQ", default=None,dest="NACQ")
    parser.add_option("-o", "--output",help="Specify output filenames additional suffix. Default=no suffix.",
                      default='',dest="OUT")
    parser.add_option("-r", "--raw",action="store_true",
                      help="Specify if we plot and save raw Qmeter data (from the BQSB class), and do a tune and rate analysis. Default=no raw data plot.",
                      metavar="RAW", default=False,dest="RAW")
    parser.add_option("-u", "--boundmin",type=float,
                      help="Specify minimum amplitude of spectral lines to begin the fit, for -q option (default=1e-7)",
                      metavar="BMIN", default=1.e7,dest="BMIN")
    parser.add_option("--frev",type=float,
                      help="Specify the revolution frequency in Hz (in case not found from files). Default = 1 MHz",
                      metavar="FREV", default=1e6,dest="FREV")
    parser.add_option("--bounds",type=float,nargs=2,
                      help="Specify the raw data fit boundaries (we fit when Qmeter between those two values). Default = [1e8,1e9].",
                      metavar="BOUNDS", default=[1e8,1e9],dest="BOUNDS")
    parser.add_option("--period",type=float,nargs=2,
                      help="Specify the number of points on which we compute the envelop (sliding maxima) and the number of points on which we do the sliding average, for the raw data fit. Default = (60,10).",
                      metavar="PERIOD", default=(60,10),dest="PERIOD")
    parser.add_option("--tunes",type=float,nargs=3,
                      help="Specify the unperturbed tunes x & y (fractional part) and the width of the spectrum around them (for analysis with -q option). Default=(0.1,0.45,0.05).",
                      metavar="TUNES", default=(0.2,0.35,0.02),dest="TUNES")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def find_power_2(n):

    # find the maximum power of 2 that is a divider of an integer n
    
    ndiv=n;
    while mod(ndiv,2)==0: ndiv /= 2;
	
    return n/ndiv;


if __name__ == "__main__":


    opt,args=parsse();
    
    # revolution frequency and period (assumes protons)
    #print 'Kinetic energy: ',opt.ENERGY,' GeV';
    E0=0.938272;
    #energy_total=opt.ENERGY+E0
    #circumference=200*np.pi; # PS circumference in m
    #gamma=energy_total/E0 # the later is the rest mass of a proton
    #beta=np.sqrt(1.-1./(gamma*gamma));
    #frev=beta*299792458./circumference;
    #Trev=1./frev; # revolution period in s
    #print 'beta=',beta,', gamma=',gamma,', Trev=',Trev*1e6,'us';

    Qx=opt.TUNES[0]; # fractional part of horizontal tune
    Qy=opt.TUNES[1]; # fractional part of vertical tune
    width=opt.TUNES[2]; # width chosen for the spectrum
    injection=275; # injection time in ms (for intensity)

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if (opt.INT):
        # initialize list with all intensity, gamma and frev data (over all files) along the cycle
	intensity=[];gamma=[];frev=[];


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
	

	if (opt.INT):
	    
	    if (not(isscalar(data.BR3.STRTR.Samples.samples+0)))and(len(data.BR3.STRTR.Samples.samples)>0):
		# extract and plot intensity along the cycle
		intensity.append(data.BR3.STRTR.Samples.samples);
		timesI=data.BR3.STRTR.Samples.firstSampleTime+np.arange(len(intensity[-1])); # times in ms
		if ('oldtimesI' in locals()):
		    if (timesI[0]!=oldtimesI[0])or(len(timesI)!=len(oldtimesI)): print "Pb: not the same intensity times in all files !", timesI,oldTimesI;sys.exit();
		oldtimesI=timesI;

		fig,ax=init_figure();
		plot(timesI,intensity[-1],'','-b','Intensity [ $ 10^{10} $ p+/b]',ax,0,lw=3.,xlab='Time along the cycle [ms]');
		end_figure(fig,ax,save=dirname+'/'+dirname+'_intensity'+loopstr+opt.OUT);

	    else:
	    	print "    No intensity data";
		intensity.append(nan);
	    

	    if (not(isscalar(data.BR.STGAMMA.Samples.samples+0)))and(len(data.BR.STGAMMA.Samples.samples)>0):
		# extract and plot intensity along the cycle
		gamma.append(data.BR.STGAMMA.Samples.samples);
		timesG=data.BR.STGAMMA.Samples.firstSampleTime+np.arange(len(gamma[-1])); # times in ms
		if ('oldtimesG' in locals()):
		    if (timesG[0]!=oldtimesG[0])or(len(timesG)!=len(oldtimesG)): print "Pb: not the same gamma times in all files !", timesG,oldTimesG;sys.exit();
		oldtimesG=timesG;

		fig,ax=init_figure();
		plot(timesG,gamma[-1],'','-b'," $ \gamma $ ",ax,0,lw=3.,xlab='Time along the cycle [ms]');
		end_figure(fig,ax,save=dirname+'/'+dirname+'_gamma'+opt.OUT);

	    else:
	    	print "    No gamma data";
		gamma.append(nan);
	    
	    
	    if (not(isscalar(data.BR.STFREVCALC.Samples.samples+0)))and(len(data.BR.STFREVCALC.Samples.samples)>0):
		# extract and plot intensity along the cycle
		frev.append(data.BR.STFREVCALC.Samples.samples);
		timesF=data.BR.STFREVCALC.Samples.firstSampleTime+np.arange(len(frev[-1])); # times in ms
		if ('oldtimesF' in locals()):
		    if (timesF[0]!=oldtimesF[0])or(len(timesF)!=len(oldtimesF)): print "Pb: not the same frevcalc times in all files !", timesF,oldTimesF;sys.exit();
		oldtimesF=timesF;

		fig,ax=init_figure();
		plot(timesF,frev[-1],'','-b',"frev",ax,0,lw=3.,xlab='Time along the cycle [ms]');
		end_figure(fig,ax,save=dirname+'/'+dirname+'_frev'+opt.OUT);

	    else:
	    	print "    No frevcalc data";
		frev.append(nan);
	    
	
	if (opt.RAW):
	
	    # extract raw Qmeter data
	    name=['Hor','Ver'];
	    # horizontal
	    raw[0].append(data.BR.BQSB.Acquisition.rawDataH);
	    # vertical
	    raw[1].append(data.BR.BQSB.Acquisition.rawDataV);
	    
	    if not(isscalar(raw[0][-1]+0)):
	    
		# nb of acquisitions
		if opt.NACQ==None:
		    try:
			nacq=data.BR.BQSB.SamplerAcquisition.nbOfMeas;
    		    except AttributeError:
		    	nacq=find_power_2(len(raw[0][-1]));
		else: nacq=opt.NACQ;
		
		nturns=len(raw[0][-1])/nacq; # nb of turns in each acquisition
		
		# times (ms) at the end of each acquisition
		try :
		    timesAcq=np.array(data.BR.BQSB.SamplerAcquisition.measStamp[:nacq],dtype=float);timelab='Time along the cycle [ms]'
		except AttributeError:
		    timesAcq=np.arange(nacq)*1000;timelab='Time along the cycle [a.u.]'
		
		# intensity (average before time of raw data)
		intensitytmp=data.BR3.STRTR.Samples.samples;
		if (not(isscalar(intensitytmp+0)))and(len(intensitytmp)>0):
		    timesItmp=data.BR3.STRTR.Samples.firstSampleTime+np.arange(len(intensitytmp)); # times in ms
		    ind=pylab.mlab.find((timesItmp>injection)*(timesItmp<timesAcq[0]+opt.BEG));
		    inten[ifile]=np.average(intensitytmp[ind]);
		else: inten[ifile]=nan;print "  warning: no intensity data.";

		# frev (at time of raw data)
		frevtmp=data.BR.STFREVCALC.Samples.samples;
		if (not(isscalar(frevtmp+0)))and(len(frevtmp)>0):
		    timesFtmp=data.BR.STFREVCALC.Samples.firstSampleTime+np.arange(len(frevtmp)); # times in ms
		    ind=pylab.mlab.find((timesFtmp>timesAcq[0])*(timesFtmp<timesAcq[0]+opt.BEG));
		    fr=np.average(frevtmp[ind])*1e3;
		else: fr=opt.FREV;print "  warning: no frevcalc data.";
		Trev=1/fr;
		print "fr=",fr;


		#print timesAcq
		if ('oldnturns' in locals()):
		    if (nturns!=oldnturns)or(len(timesAcq)!=len(oldtimesAcq)): print "Pb: not the same nb of points and/or times in raw data in all files !",nturns,oldnturns,timesAcq,oldtimesAcq;sys.exit();
		    elif (np.abs(timesAcq-oldtimesAcq)>1).any(): print "Pb: not the same times in raw data in all files !",np.abs(timesAcq-oldtimesAcq),timesAcq,oldtimesAcq;sys.exit();
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
		    plot(timesraw,raw[i][-1],'Q meter','.b','Qmeter raw '+name[i]+'. data [a.u.]',axraw[i],0,lw=2.5,xlab=timelab,plotevery=1);
		    
		    # do an exponential fit of the envelop
		    #print timesraw,injtime,opt.BEG
		    tau=fit(timesraw,raw[i][-1],axraw[i],True,'Q meter','Qmeter raw '+name[i]+'. data [a.u.]',0.001,opt.BOUNDS,0,opt.PERIOD,opt.BEG+injtime,0,col='r');
		    axraw[i].set_xlabel(timelab);
		    end_figure(figraw[i],axraw[i],save=dirname+'/'+dirname+'_Qmeter_raw_'+name[i]+loopstr+opt.OUT);
		    
		    if (i==0): ratesfitx[ifile]=1./tau;
		    else: ratesfity[ifile]=1./tau;

	    	if False:
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
			    lastturn=lastx,flagaver=True,xlab=timelab)
		    tunesy[ifile],tauy,sigtunesy[ifile]=collect_and_fit(tuneby[-1],ampby,
			    newtimesAcq,'SUSSIX, vertical','y',opt.BMIN,True,axtuy,axampy,
			    0.001,'Tune','Amplitude of the tune line',firstturn=first,
			    lastturn=lasty,flagaver=True,xlab=timelab)
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
	    #print >> filetunes, "Intensity[10^10]\tGrowthratex[s-1]\tGrowthratey[s-1]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]\tTune_x\tTune_y\tSigma_tune_x\tSigma_tune_y"
	    #print >> filetunes, inten[ifile], "\t", ratesx[ifile], "\t", ratesy[ifile],"\t", ratesfitx[ifile], "\t", ratesfity[ifile],"\t", tunesx[ifile], "\t", tunesy[ifile],"\t", sigtunesx[ifile], "\t", sigtunesy[ifile];
	    print >> filetunes, "Intensity[10^10]]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]";
	    print >> filetunes, inten[ifile], "\t", ratesfitx[ifile], "\t", ratesfity[ifile],;
	    filetunes.close();
		

    # for final averages
    strdamper=[''];legdamper=[''];

    # define colors for average plots (closed loop = first color, openloop = second color);
    col=['-r','-b'];
    
    if (opt.INT):
	# compute and plot the average intensity (of all files), along the cycle
	
	fig,ax=init_figure();
	
	data1=np.array([el for el in intensity if isfinite(el).any()])
	intensityav=np.average(data1,axis=0);
	    
	if not(isnan(intensityav).all()): plot(timesI,intensityav,'','b','Average intensity [ $ 10^{10} $ p+/b]',ax,0,lw=3.,xlab='Time along the cycle [ms]');
	
	end_figure(fig,ax,save='average_intensity'+opt.OUT);

    
    if (opt.RAW):
	# compute and plot the average raw data signal
    	
	name=['Hor','Ver'];
	
	for i in range(2):
	    
	    if len(raw[i])>0:
	    
		figav,axav=init_figure();
		figmax,axmax=init_figure();
	    
		data1=np.array([el for el in raw[i] if isfinite(el).any()])
		rawav=np.average(data1,axis=0);
		rawmax=np.max(np.abs(data1),axis=0);

		if not(isnan(rawav).all()):
		    plot(timesraw,rawav,'','b','Average Qmeter raw '+name[i]+'. data [a.u.]',axav,0,lw=2.5,xlab=timelab,plotevery=1);
		    plot(timesraw,rawmax,'','b','Max. Qmeter raw '+name[i]+'. data [a.u.]',axmax,0,lw=2.5,xlab=timelab,plotevery=1);

		end_figure(figav,axav,save='average_Qmeter_raw_'+name[i]+opt.OUT);
		end_figure(figmax,axmax,save='max_Qmeter_raw_'+name[i]+opt.OUT);

   
 	# compute the average and sigmas over all files of the intensity, tunes and growth rates

	figx,axx=init_figure();
	figy,axy=init_figure();

	# take out nan and infinities
	inten1=inten[pylab.mlab.find(isfinite(inten))];
	ratesx1=ratesx[pylab.mlab.find(isfinite(ratesx))];
	ratesy1=ratesy[pylab.mlab.find(isfinite(ratesy))];
	ratesfitx1=ratesfitx[pylab.mlab.find(isfinite(ratesfitx))];
	ratesfity1=ratesfity[pylab.mlab.find(isfinite(ratesfity))];
	#sigtunesx1=sigtunesx[pylab.mlab.find(isfinite(tunesx))];
	#sigtunesy1=sigtunesy[pylab.mlab.find(isfinite(tunesy))];
	#tunesx1=tunesx[pylab.mlab.find(isfinite(tunesx))];
	#tunesy1=tunesy[pylab.mlab.find(isfinite(tunesy))];

	# compute averages and standard deviation
	inten_av=np.average(inten1);inten_sig=np.sqrt(np.var(inten1));
	#ratesx_av=np.average(ratesx1);ratesx_sig=np.sqrt(np.var(ratesx1));
	#ratesy_av=np.average(ratesy1);ratesy_sig=np.sqrt(np.var(ratesy1));
	ratesfitx_av=np.average(ratesfitx1);ratesfitx_sig=np.sqrt(np.var(ratesfitx1));
	ratesfity_av=np.average(ratesfity1);ratesfity_sig=np.sqrt(np.var(ratesfity1));
	#tunesx_av=np.average(tunesx1);
	#tunesy_av=np.average(tunesy1);
	#tunesx_sig=np.sqrt((np.sum(sigtunesx1**2)+np.sum(tunesx1**2))/float(len(tunesx1))-tunesx_av**2);
	#tunesy_sig=np.sqrt((np.sum(sigtunesy1**2)+np.sum(tunesy1**2))/float(len(tunesy1))-tunesy_av**2);

	# output
	filefinal=open('final_tunes_and_growthrates'+opt.OUT+'.txt','w')
	#print >> filefinal, "Intensity[10^10]\tSigma_intensity[10^10]\tGrowthratex[s-1]\tGrowthratey[s-1]\tSigma_rate_x[s-1]\tSigma_rate_y[s-1]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]\tSigma_rate_fitx[s-1]\tSigma_rate_fity[s-1]\tTune_x\tTune_y\tSigma_tune_x\tSigma_tune_y"
	#print >> filefinal, inten_av, "\t", inten_sig, "\t", ratesx_av, "\t", ratesy_av,"\t", ratesx_sig, "\t", ratesy_sig,"\t", ratesfitx_av, "\t", ratesfity_av,"\t", ratesfitx_sig, "\t", ratesfity_sig,"\t", tunesx_av, "\t", tunesy_av,"\t", tunesx_sig, "\t", tunesy_sig;
	print >> filefinal, "Intensity[10^10]\tSigma_intensity[10^10]\tGrowthratefitx[s-1]\tGrowthratefity[s-1]\tSigma_rate_fitx[s-1]\tSigma_rate_fity[s-1]"
	print >> filefinal, inten_av, "\t", inten_sig, "\t", ratesfitx_av, "\t", ratesfity_av,"\t", ratesfitx_sig, "\t", ratesfity_sig;
	filefinal.close();

	if False:
	    # now compute and plot average tunes along the cycle
	    datax=np.array([el for el in tunebx if isfinite(el).any()])
	    datay=np.array([el for el in tuneby if isfinite(el).any()])
	    tunexav=np.average(datax,axis=0);
	    tuneyav=np.average(datay,axis=0);
	    tunexvar=np.sqrt(np.var(datax,axis=0));
	    tuneyvar=np.sqrt(np.var(datay,axis=0));

	    if not(isnan(tunexav).all()): axx.errorbar(newtimesAcq,tunexav,yerr=tunexvar,fmt='b',label='',lw=3);
	    if not(isnan(tuneyav).all()): axy.errorbar(newtimesAcq,tuneyav,yerr=tuneyvar,fmt='b',label='',lw=3);

	    axx.set_xlabel(timelab);
	    axx.set_ylabel('Average hor. tune');
	    axy.set_xlabel(timelab);
	    axy.set_ylabel('Average ver. tune');
	    end_figure(figx,axx,save='average_tunesx'+opt.OUT);
	    end_figure(figy,axy,save='average_tunesy'+opt.OUT);
	    
