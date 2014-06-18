#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,re,dateutil,random,pytz
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from optparse import OptionParser
from Timber import parseout
import math
import matplotlib
import matplotlib.dates
#import subprocess
import glob
from plot_lib import plot,set_fontsize,init_figure,end_figure
from string_lib import get_nice_string
from datetime_lib import tt
from io_lib import list_files
from read_Headtail_prt_fit import slideaver
from collimator_settings import concatenate_data_Timber

def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--add",action="store_true",
                      help="Specify if we add up all variables selected with a single -l or -r option (useful for RF voltage).",
                      metavar="ADD", default=False,dest="ADD")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the TIMBER .csv name of the file. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      #help="Specify the TIMBER .csv name of the file (several -f options possible -> several files); it can also be regular expressions encompassing several files at a time (with e.g. '*', etc.)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--legend",action="append",
                      help="Specify the legend for each variable (first those for left axis then those for right axis); by default take the variable name",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-l", "--leftaxis",action="append",
                      help="Specify the TIMBER variable name to plot on the left axis (several -l options possible -> several variables); only a part of the variable name (sufficient to identify it) is enough.",
                      metavar="LEFT", default=None,dest="LEFT")
    parser.add_option("-o", "--output",help="Specify filename for text file where variables are output. Default: no text file output.",
                      default=None,dest="OUT")
    parser.add_option("-r", "--rightaxis",action="append",
                      help="Specify the TIMBER variable name to plot on the right axis (several -r options possible -> several variables); only a part of the variable name (sufficient to identify it) is enough",
                      metavar="RIGHT", default=None,dest="RIGHT")
    parser.add_option("-s", "--slideaver",action="append",
                      help="Specify the period for a sliding average. Two numbers maximum (in two separate -s options) : one for the left axis + one for the right axis - if needed; by default no sliding average. 0 or 1 also mean no sliding average.",
                      metavar="SLIDE", default=None,dest="SLIDE")
    parser.add_option("-t", "--time",nargs=12,
                      help="Specify between which dates and times you want to plot and (with -o option) output the curves in a .txt file (year, month, day, hour, minutes, seconds - twice). Default: take all the times present in the initial file(s).",
                      type=int,metavar="TIME", default=None,dest="TIME")
    parser.add_option("-v", "--vs",action="store_true",
                      help="Specify if we do additional plot with right-axis variables as a function of the first left-axis one.",
                      metavar="VS", default=False,dest="VS")
    parser.add_option("-y", "--ylabel",action="append",
                      help="Specify the y axis name. Two labels maximum (in two separate -y options) : one for the left axis + one for the right axis - if needed; by default take the name of the first variable on the axis",
                      metavar="YLABEL", default=None,dest="YLABEL")
    parser.add_option("--sigl",action="append",
                      help="Specify the TIMBER variable name to use as an error bar for the left axis curves (several --sigl options possible -> one for each -l option); only a part of the variable name (sufficient to identify it) is enough. Default: no error bar. Use 0 to avoid giving error bars for a certain curve. Not compatible with -a or -s option.",
                      metavar="SIGLEFT", default=None,dest="SIGLEFT")
    parser.add_option("--sigr",action="append",
                      help="Specify the TIMBER variable name to use as an error bar for the right axis curves (several --sigr options possible -> one for each -r option); only a part of the variable name (sufficient to identify it) is enough. Default: no error bar. Use 0 to avoid giving error bars for a certain curve. Not compatible with -a or -s option.",
                      metavar="SIGRIGHT", default=None,dest="SIGRIGHT")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def extract_between_times(times,tarray,datarray,variable,dispflag=False,energy=4000.):

    # extract between two times given as an array [year1, month1, day1, hour1, minutes1, seconds1, year2, month2, day2, hour2, minutes2, seconds2];
    # variable is the TIMBER variable name
    
    tim1=datetime(times[0],times[1],times[2],times[3],times[4],times[5]);
    tim2=datetime(times[6],times[7],times[8],times[9],times[10],times[11]);
    dt1=pylab.date2num(tim1);
    dt2=pylab.date2num(tim2);
    if (dispflag):
	print "Time after which data is taken: ",pylab.num2date(dt1);
	print "Time before which data is taken: ",pylab.num2date(dt2);
    
    ind=np.where((tarray>(dt1*86400))*(tarray<(dt2*86400)));
    
    if (len(ind)>0)and(len(ind[0])>0):
    	ind=ind[0];
	t=tarray[ind];data=dataarray[ind];

	# compute average
	while (variable.endswith(" ")): variable=variable.rstrip(" ");
	variable=replace(variable,":","_");
	fileav=open('average_'+variable+'_'+tim1.strftime("%Y-%m-%d_%H-%M-%S")+'_'+tim2.strftime("%Y-%m-%d_%H-%M-%S")+'.txt','w');
	av=np.average(data);
	sigma=np.sqrt(np.average(data*data)-av**2);
	if (variable.find('ANA_GAIN')!=-1):
	    # particular case of damper gain - we use formula provided by D. Valuch to convert
	    # gain in DB in number of damping turns ntau.
	    # NOTE: coefficients 'a' should be changed sometimes 
	    # (below are values for after July 23rd, 2012)
	    if variable.find('B1')!=-1:
		if variable.startswith('ADTH'): a=10000.;
		else: a=9199.;
	    else:
		if variable.startswith('ADTH'): a=7942.;
		else: a=6197.;
	    print >> fileav, variable,"\tsigma_gain_dB\tntau";
	    ntau=2./(4095*10**(av/20.)/(a*energy/450.));
	    print >> fileav, av,"\t",sigma,"\t",ntau;
	else:
	    print >> fileav, variable,"\tsigma";
	    print  >> fileav, av,"\t",sigma;

	fileav.close();

    else:
	print "Pb in extract_between_times: no data found";
	t=[];data=[];

    return t, data;


if __name__ == "__main__":
    opt,args=parsse();

    if (opt.YLABEL!=None)and(len(opt.YLABEL)>2): print "Too many labels for the y axes ! Put only 2 !";

    gmt=pytz.timezone('Europe/Amsterdam');

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
	
    # concatenate datafiles in one single data set
    data=concatenate_data_Timber(listname);
    
    # initialize plot(s)
    fig,ax=init_figure();
    if (opt.VS): figvs,axvs=init_figure();
    
    color=['b','r','g','m','c','y','k'];

    if (opt.OUT!=None): dataout=[];
    
    # identify and plot the variables for the left axis
    for ivar,var in enumerate(opt.LEFT):
	
	if (opt.YLABEL==None): ylab=var;
	else: ylab=opt.YLABEL[0];
	
	v_count=0; # to count number of plotted variables
	for v in data['datavars']:
	
	    if (v.find(var)!=-1):
	    
		# legend
		leg='';
		if (v_count==0):
		    if (opt.LEG!=None)and(len(opt.LEG)>ivar): leg=opt.LEG[ivar];
		    else: leg=v;
		
		# make arrays with TIMBER data
		dataarray=np.array([float(j[0]) for j in data[v][1]]);
		tps=data[v][0];
	    	tarray=np.array([pylab.date2num(tt(j[:-4])) for j in tps])*86400.+np.array([float(j[-4:]) for j in tps]);
		
		if (opt.TIME!=None):
		    # extract between given times
		    ttmp,datatmp=extract_between_times(opt.TIME,tarray,dataarray,v,dispflag=(ivar==0),energy=opt.EN);
		    tarray=ttmp;dataarray=datatmp;

		if (opt.OUT!=None):
		    if ivar==0: tout=tarray;dataout.append(dataarray);
		    else: dataout.append(np.interp(tout,tarray,dataarray));
		
		if (opt.ADD):
		    if (v_count==0):
			# initialize summ arrays
			tsum=tarray;datasum=dataarray;
		    else:
		    	# add array to the total sum
			datasum += np.interp(tsum,tarray,dataarray);
		
		if (opt.SLIDE!=None)and(int(opt.SLIDE[0])>1):
		    # do a sliding average
		    dataarray=slideaver(dataarray,int(opt.SLIDE[0]));
		    tarray=slideaver(tarray,int(opt.SLIDE[0]));
		
		# extract date (first variable only)
		if (ivar==0):
		    dat=pylab.num2date(tarray[0]/86400.).date();
		    timeinterval=int(np.ceil((tarray[-1]-tarray[0])/8));
		
		if not(opt.ADD):
		    if (opt.SIGLEFT==None)or((len(opt.SIGLEFT)<=ivar)or(opt.SIGLEFT[ivar]==0)):
			# plot vs time
			plot(tarray/86400.,dataarray,leg,color[ivar],ylab,ax,0,lw=4,xlab="Time on "+dat.strftime("%Y-%m-%d"));
		    else:
		        # plot with error bars vs time
			for vsig in data['datavars']:
	    		    if (vsig.find(opt.SIGLEFT[ivar])!=-1): sigarray=np.array([float(j[0]) for j in data[vsig][1]]);
			ax.errorbar(tarray/86400,dataarray,yerr=sigarray,fmt=color[ivar],label=leg,lw=2.);
			ax.set_xlabel("Time on "+dat.strftime("%Y-%m-%d"));
			ax.set_ylabel(ylab);
			
		if (opt.VS)and(ivar==0):
		    ylabl=ylab;
		    tarrayl=tarray;
		    dataarrayl=dataarray;
			
		v_count += 1;

	if opt.ADD:
	    # plot vs date the total sum
	    plot(tsum/86400.,datasum,leg,color[ivar],ylab,ax,0,lw=4,xlab="Time on "+dat.strftime("%Y-%m-%d"));
		

		    

    ivar+=1;
    
    if (opt.RIGHT!=None):
	
	ylim=ax.get_ylim();
	ax2=pylab.twinx(ax=ax); # second y axis
	
	# identify and plot the variables for the left axis
	for jvar,var in enumerate(opt.RIGHT):

	    if (opt.YLABEL==None)or(len(opt.YLABEL)<2): ylab=var;
	    else: ylab=opt.YLABEL[1];

	    v_count=0; # to count number of plotted variables
	    for v in data['datavars']:
	    
		if (v.find(var)!=-1):

		    # legend
		    leg='';
		    if (v_count==0):
			if (opt.LEG!=None)and(len(opt.LEG)>ivar+jvar): leg=opt.LEG[ivar+jvar];
			else : leg=v;
		    
		    # make arrays with TIMBER data
		    dataarray=np.array([float(j[0]) for j in data[v][1]]);
		    tps=data[v][0];
	    	    tarray=np.array([pylab.date2num(tt(j[:-4])) for j in tps])*86400.+np.array([float(j[-4:]) for j in tps]);

		    if (opt.TIME!=None):
			# extract between given times
			ttmp,datatmp=extract_between_times(opt.TIME,tarray,dataarray,v,dispflag=False,energy=opt.EN);
			tarray=ttmp;dataarray=datatmp;

		    if (opt.OUT!=None):
			dataout.append(np.interp(tout,tarray,dataarray));
		
		    if (opt.ADD):
			if (v_count==0):
			    # initialize summ arrays
			    tsum=tarray;datasum=dataarray;
			else:
		    	    # add array to the total sum
			    datasum += np.interp(tsum,tarray,dataarray);
		
		    if (opt.SLIDE!=None)and(int(opt.SLIDE[1])>1):
			# do a sliding average
			dataarray=slideaver(dataarray,int(opt.SLIDE[1]));
			tarray=slideaver(tarray,int(opt.SLIDE[1]));

		    if not(opt.ADD):
			if (opt.SIGRIGHT==None)or((len(opt.SIGRIGHT)<=jvar)or(opt.SIGRIGHT[jvar]==0)):
			    # plot vs time
			    plot(tarray/86400.,dataarray,leg,color[ivar+jvar],ylab,ax2,0,lw=4,xlab="Time on "+dat.strftime("%Y-%m-%d"));
			else:
		            # plot with error bars vs time
			    for vsig in data['datavars']:
	    			if (vsig.find(opt.SIGRIGHT[jvar])!=-1): sigarray=np.array([float(j[0]) for j in data[vsig][1]]);
			    ax2.errorbar(tarray/86400,dataarray,yerr=sigarray,fmt=color[ivar+jvar],label=leg,lw=2.);
			    ax2.set_xlabel("Time on "+dat.strftime("%Y-%m-%d"));
			    ax2.set_ylabel(ylab);
	
		    if (opt.VS):
		    	# plot this right-axis variable vs. first left-axis one (on a separate figure)
			dataarrayr=np.interp(tarrayl,tarray,dataarray);
		    	plot(dataarrayl,dataarrayr,leg,'.'+color[ivar+jvar],ylab,axvs,0,lw=2.5,xlab=ylabl);
			
		    v_count += 1;
					    
	    if opt.ADD:
		# plot vs date the total sum
		plot(tsum/86400.,datasum,leg,color[ivar+jvar],ylab,ax,0,lw=4,xlab="Time on "+dat.strftime("%Y-%m-%d"));
		
    	
    	jvar+=1;
	ax.set_ylim(ylim);end_figure(fig,ax2,legpos=3);
    
    
    if (timeinterval<60):
	ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval));
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'));
    elif (timeinterval<3600):
	ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=int(np.floor(timeinterval/60.))));
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'));
    else:
	ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=int(np.floor(timeinterval/3600.))));
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'));
    
    end_figure(fig,ax,legpos=2);
    if (opt.VS): end_figure(figvs,axvs);
    
    if (opt.OUT!=None):
    	# write extracted data to file
	dataout=np.array(dataout).transpose();
    	filedata=open(opt.OUT,'w');
	if (opt.RIGHT==None):
	    print >> filedata, get_nice_string(opt.LEFT);
	else:
	    print >> filedata, get_nice_string(opt.LEFT), get_nice_string(opt.RIGHT);
	
	for it,tmp in enumerate(tout):
	    print >> filedata, get_nice_string(dataout[it]);
	
	filedata.close();
	
 	    
    pylab.show();
    sys.exit();
    
    
