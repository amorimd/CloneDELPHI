#!/usr/bin/python2

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,random
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from parser_lib import *
import math
import matplotlib
from read_cfg import read_cfg
from plot_lib import set_fontsize,init_figure,cmap,plot,end_figure
from io_lib import list_files
from read_Headtail_prt import extract_Headtail_param,check_data
from string_lib import takeout_common,takeout_spaces


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--analysis",action="store_true",
                      help="Specify if a data analysis has to be performed",
                      metavar="ANA", default=False,dest="ANA")
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the raw data plot (several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify boundaries for the fit of individual bunches (default=[1e-5, 1e50]): fit only when x & y are between those two values. Note: inverse the order to fit a damping time.",
                      metavar="BOUND", default=[1.e-5,1.e50],dest="BOUND")		      
    parser.add_option("-e", "--legend",action="append",
                      help="Specify the legend for the plot, for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-m", "--boundbeam",type=float,nargs=2,
                      help="Specify boundaries for the fit of the whole beam (default=[1e-6, 1e50]): fit only when x & y are between those two values. Note: inverse the order to fit a damping time.",
                      metavar="BOUNDBM", default=[1.e-6,1.e50],dest="BOUNDBM")		      
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the fit (default=0)",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-l", "--plot",type=int,
                      help="Specify 1 to plot also maxima, 2 to plot maxima and sliding average, 3 to plot sliding average (it also always plots raw data and fit)",
                      metavar="PLOT", default=0,dest="PLOT")
    parser.add_option("-o", "--output",help="Specify output suffix for rise times",
                      default="tau.txt",dest="OUT")
    parser.add_option("-p", "--plotlog",type=int,
                      help="Specify if we want to plot in linear scale (0), in semilogx (1), in semilogy (2) or in loglog (3)",
                      metavar="LOG", default=0,dest="LOG")
    parser.add_option("-r", "--period",type=int,nargs=2,
                      help="Specify periods used for 1) extracting the envelop and 2) the sliding average (default =[60,1])",
                      metavar="PERIOD", default=[60,1],dest="PERIOD")		      
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot (and fit with -a option) the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--plotraw",type=int,
                      help="Specify some option for the plot of the raw data: 1 (default)=plot all turns vs. number of turns, n=plot every n turns vs. time",
                      metavar="PLOTRAW", default=1,dest="PLOTRAW")
    # Note: options -c, -o, -l, -d, -r and -g used only when -a option activated	      
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    print "Legends:", opt.LEG
    return opt, args


def read_prt_file(filename,nbunch):

    # read prt file
    print filename

    #x=[];xp=[];y=[];yp=[];bunch_nos=[];
    for l,line in enumerate(open(filename)):
	#try:
	#    i=int(split(line)[20])
	#except IndexError:
	if (len(split(line))<21):
	    print 'Line ',l+1,': not enough columns!'
	    l-=1;
	    break

    # check that the data is not somehow corrupted (missing / added lines, etc.)
    print "length of data",l+1
    lmod=check_data(l+1,nbunch);
    # lmod is the number of line to take away at the end (when length of data is not a multiple of nbunch)

    # initialize arrays
    x=np.zeros(l+1-lmod)
    xp=np.zeros(l+1-lmod)
    y=np.zeros(l+1-lmod)
    yp=np.zeros(l+1-lmod)

    # read file
    for il,line in enumerate(open(filename)):
    	if (il<(l+1-lmod)):
	    x[il]=float(split(line)[1])
	    xp[il]=float(split(line)[2])
	    y[il]=float(split(line)[3])
	    yp[il]=float(split(line)[4])
	else: break
    
    #print len(x),len(y),len(xp),len(yp)
    
    return x,xp,y,yp


def envelop(turns,data,period):

    # extract the envelop of maxima of the absolute values of data,
    # (sliding maxima over nb points=period)
    env=[];turnsenv=[];
    for k,j in enumerate(data[:-period-1:period]):
	env.append(np.max(np.abs(data[k*period:(k+1)*period+1])));
	turnsenv.append(turns[np.argmax(np.abs(data[k*period:(k+1)*period+1]))+k*period]);
    
    return np.array(turnsenv),np.array(env);



def slideaver(data,period):
    # do a sliding average (average over nb points=period)
    
    # old version (very slow)
    #return np.array([np.average(data[k:min(k+period,len(data))]) for k,j in enumerate(data)]);
    
    # faster version
    result=np.zeros(len(data));
    for idat,dat in enumerate(data):
    	if (np.mod(idat,1000)==0)or((len(data)-idat)<=period):
	    # every 1000 steps (or close to the end) we compute from scratch the average (to avoid accumulation of errors)
	    result[idat]=np.average(data[idat:min(idat+period,len(data))]);
	else:
	    result[idat]=(float(period)*result[idat-1]-data[idat-1]+data[idat+period-1])/(float(period));

    return result;


def fit(turns,data,ax,flag,leg,ylab,Trev,bound,optplot,optperiod,optbeg,optlog,nevery=1,col='',datet=None):

    # routine to fit data by an exponential

    # Input: data is the data to fit (1D array, it's assumed to be vs. indices numbers,
    # i.e. turns), ax the axes on which to plot, flag is True if we want to plot, leg is the
    # legend to add to the plot, ylab the y axis label, Trev the revolution time 
    # (for conversion of the rise time into seconds), bound (2 elements array) the bounds 
    # in between which the data should be to be fitted
    
    # when datet is not None, turns is interpreted as times instead of turn numbers,
    # plot vs time instead of vs turns, and datet contains
    # the date (datetime format) (for the x-axis label)
    if (datet!=None):
    	timefactor=Trev/86400.;t0=turns[0]/86400.
	turns=(turns-turns[0])/Trev; # conversion from times to turn numbers
	labx="Time on "+datet.strftime("%Y-%m-%d")
    else:
    	timefactor=1.;t0=0.;
	labx="Number of turns";

    # take out nan and inf
    ind=~np.isnan(data);data=data[ind];turns=turns[ind];
    ind=~np.isinf(data);data=data[ind];turns=turns[ind];
    
    # Output: rise time in seconds
    
    # first step: construct "envelop" (i.e. curve of maxima) of the absolute value
    period=optperiod[0]; # number of turns on which we extract maxima
    turnsenv,env=envelop(turns,data,period);

    if flag:
	if (optplot==1)or(optplot==2):
	    # plot envelop
	    plot(turnsenv*timefactor+t0,env,leg+', envelop',':'+col,ylab,ax,optlog,plotevery=nevery,xlab=labx);


    # second step: do a "sliding average"
    periodaver=optperiod[1]; # number of points on which we average
    envaver=slideaver(env,periodaver);
    turnsenvaver=slideaver(turnsenv,periodaver);

    if flag:
	if (optplot==2)or(optplot==3):
	    plot(turnsenvaver*timefactor+t0,envaver,leg+', envelop average','-'+col,ylab,ax,optlog,plotevery=nevery,xlab=labx);

    # third step: find out bounds between which we should fit
    boundmin=bound[0]; # fit when signal becomes higher than this
    boundmax=bound[1]; # stop fit when we reach this value
    begt=np.where(turnsenvaver>=optbeg);begt=begt[0]; # for some reason begt itself is a single-element array of array
    
    if (boundmin<=boundmax):
    
	beg=np.where(envaver[begt[0]:-1]>=boundmin);
	beg=beg[0]+begt[0]; # for some reason beg itself is a single-element array of array
	if (len(beg)>0):
            if (boundmax>np.max(envaver[beg[0]:-1])): boundmax=np.max(envaver[beg[0]:-1]);
	    end=np.where(envaver[beg[0]:-1]>=boundmax);
	    if (len(end[0])==0): end=len(turnsenvaver)+beg[0];
	    else: end=end[0]+beg[0];end=end[0];

	    if (end-beg[0]>0):
		# fourth and final step: fit
		p=np.polyfit(turnsenvaver[beg[0]:end],np.log(envaver[beg[0]:end]),1);
		# note: the fit goes from beg[0] to end[0]-1
		tau=1./p[0];a=math.exp(p[1]); # tau is in number of turns here
		#print a,tau,p
		#sys.stdout.write('Rise time for file %s, bunch no %d : %lf s\n' % (fil, bunch, tau*Trev) );
		if flag:
		    astr="%0.3e" %(tau*Trev);bstr=astr.replace("e-0","\cdot 10^{-").replace("e+0","\cdot 10^{");
		    bstr=r"$"+bstr+"}$";
    		    plot(turnsenvaver[beg[0]:end]*timefactor+t0,np.exp(turnsenvaver[beg[0]:end]/tau+p[1]),leg+r', fit, $\tau=$'+bstr+' s','-'+col,ylab,ax,optlog,lw=4,plotevery=nevery,xlab=labx);
		    #print np.exp(turnsenvaver[beg[0]:end]/tau+p[1]),turnsenvaver[beg[0]:end]*timefactor+t0,np.exp(turnsenvaver[beg[0]:end]/tau+p[1])
	    
	    else: tau=float('inf');

        else: tau=float('inf');
	
    else:
    
	# the other way around to fit a damping rate
	boundmin=bound[1];boundmax=bound[0];
	beg=np.where(envaver[begt[0]:-1]<=boundmax);
	beg=beg[0]+begt[0]; # for some reason beg itself is a single-element array of array
	if (len(beg)>0):
            if (boundmin<np.min(envaver[beg[0]:])): boundmin=np.min(envaver[beg[0]:]);
	    end=np.where(envaver[beg[0]:]<=boundmin);
	    if (len(end[0])==0): end=len(turnsenvaver)+beg[0];
	    else: end=end[0]+beg[0];end=end[0];

	    if (end-beg[0]>1):
		# fourth and final step: fit
		p=np.polyfit(turnsenvaver[beg[0]:end],np.log(envaver[beg[0]:end]),1);
		# note: the fit goes from beg[0] to end[0]-1
		tau=1./p[0];a=math.exp(p[1]); # tau is in number of turns here
		#print a,tau,p,end-beg[0]
		#sys.stdout.write('Rise time for file %s, bunch no %d : %lf s\n' % (fil, bunch, tau*Trev) );
		if flag:
		    astr="%0.3e" %(tau*Trev);bstr=astr.replace("e-0","\cdot 10^{-").replace("e+0","\cdot 10^{");
		    bstr=r"$"+bstr+"}$";
    		    plot(turnsenvaver[beg[0]:end]*timefactor+t0,a*np.exp(turnsenvaver[beg[0]:end]/tau),leg+r', fit, $\tau=$'+bstr+' s','-'+col,ylab,ax,optlog,lw=4,plotevery=nevery,xlab=labx);

	    else: tau=float('-inf');

        else: tau=float('-inf');

    return tau*Trev;



if __name__ == "__main__":
    opt,args=parsse(); 

    if (opt.BNUM==None): bunches=[];
    else: 
        bunches=np.array(opt.BNUM);
        bunches.sort();
    print "Selected bunches:", bunches

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if (len(listname)>4):
    	axes=[0.12,0.13,0.55,0.8];legpos=(1.03,-0.1);
    else:
    	axes=None;legpos=0;
    
    # create list of associated legends (either with opt.LEG or 
    # the list of file names taking out all the common parameters
    # in the names)
    if (opt.LEG!=None):
        listleg=opt.LEG;
    else:
        listleg=takeout_common(listname);
	listleg=takeout_spaces(listleg);
    
    
    if (len(bunches)>0):
	# initialize plots
	figx,axx=init_figure(axes=axes);
	figy,axy=init_figure(axes=axes);
    else: axx=0;axy=0;

    if (opt.AVER):
        figbeamx,axbeamx=init_figure(axes=axes);
        figbeamy,axbeamy=init_figure(axes=axes);


    if (opt.ANA):
    	figtau,axtau=init_figure(axes=axes);
    
    # loop on filenames
    for ifile,filename in enumerate(listname):
    

	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    
	if (opt.ANA):
    	    taux=np.zeros(nbunch);
	    tauy=np.zeros(nbunch);
	    file=open(filename+opt.OUT,'w');
	    print >> file, "Bunch\ttaux[s]\ttauy[s]"

    	# for the legend
        #fil=filename.replace("_"," ").replace("prt.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read prt file
	x,xp,y,yp=read_prt_file(filename,nbunch)

        # extract data and plot each bunch chosen
	for bnum in range(nbunch,0,-1):
	
	    x1=x[bnum-1::nbunch];
    	    y1=y[bnum-1::nbunch];
	    turns=np.arange(0,len(x1),dtype=float);
	    
    	    if (bnum in bunches):
		# plot raw data (average x)
    		plot(turns,x1,fil+', bunch '+str(bnum),'.',"Average x position of the bunch [m]",axx,opt.LOG,plotevery=opt.PLOTRAW);
    		# plot raw data (average y)
    		plot(turns,y1,fil+', bunch '+str(bnum),'.',"Average y position of the bunch [m]",axy,opt.LOG,plotevery=opt.PLOTRAW);

	    if (opt.AVER):
	    	if (bnum==nbunch):
		    xave=x1;
		    yave=y1;
		else:
		    xave=x1+xave;	    
		    yave=y1+yave;	    

    	    if (opt.ANA):
		# data analysis (fit)
		flag=(bnum in bunches);
		taux[bnum-1]=fit(turns,x1,axx,flag,fil+', bunch '+str(bnum),"Average x position of the bunch [m]",Trev,opt.BOUND,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG,nevery=opt.PLOTRAW);
		tauy[bnum-1]=fit(turns,y1,axy,flag,fil+', bunch '+str(bnum),"Average y position of the bunch [m]",Trev,opt.BOUND,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG,nevery=opt.PLOTRAW);
		print >> file, bnum, "\t", taux[bnum-1], "\t", tauy[bnum-1];
		
	if (opt.AVER):
	    xave/=float(nbunch);	    
	    yave/=float(nbunch);
	    turns=np.arange(0,len(xave),dtype=float);
	    # plot raw data (average x)
    	    plot(turns,xave,fil,'.',"Average x position of the beam [m]",axbeamx,opt.LOG,plotevery=opt.PLOTRAW);
    	    # plot raw data (average y)
    	    plot(turns,yave,fil,'.',"Average y position of the beam [m]",axbeamy,opt.LOG,plotevery=opt.PLOTRAW);
	    
	    if (opt.ANA):
	    	# fit
		tauxbeam=fit(turns,xave,axbeamx,True,fil,"Average x position of the beam [m]",Trev,opt.BOUNDBM,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG,nevery=opt.PLOTRAW);
		tauybeam=fit(turns,yave,axbeamy,True,fil,"Average y position of the beam [m]",Trev,opt.BOUNDBM,opt.PLOT,opt.PERIOD,opt.BEG,opt.LOG,nevery=opt.PLOTRAW);
		print fil, " - Average rise times of the beam (x & y) in sec.: ", tauxbeam, "  ", tauybeam;
	    

	if (opt.ANA):
    	    file.close()
	    axtau.plot(range(1,nbunch+1),taux,'-',label=r'$\tau_x$'+', '+fil,lw=3);
	    axtau.plot(range(1,nbunch+1),tauy,'-',label=r'$\tau_y$'+', '+fil,lw=3);


    if (opt.ANA):
	axtau.set_xlabel("Bunch number");
	axtau.set_ylabel("Rise time [s]");
	end_figure(figtau,axtau,legpos=legpos);


    if (len(bunches)>0):
	end_figure(figx,axx,legpos=legpos,legfontsize=max(30-5*len(listname),10));
	end_figure(figy,axy,legpos=legpos,legfontsize=max(30-5*len(listname),10));
    if (opt.AVER):
	end_figure(figbeamx,axbeamx,legpos=legpos,legfontsize=max(30-5*len(listname),10));
	end_figure(figbeamy,axbeamy,legpos=legpos,legfontsize=max(30-5*len(listname),10));

    pylab.show();

    sys.exit()

