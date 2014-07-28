import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

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
from datetime_lib import tt
#import subprocess
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-p", "--plane",
                      help="Specify plane: h or v",
                      metavar="PLANE", default=None,dest="PLANE")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-n", "--nbfile",type=int,
                      help="Specify maximum number of BBQ data files to take",
                      metavar="NFILE", default=10000,dest="NFILE")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-m", "--multi",action="store_true",
                      help="Specify if there are several files saved with the multifile option in Timber",
                      metavar="MULT", default=False,dest="MULT")
    parser.add_option("-a", "--analysis",action="store_true",
                      help="Specify if a data analysis has to be performed",
                      metavar="ANA", default=False,dest="ANA")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify boundaries for the fit",
                      metavar="BOUND", default=[6.e7,1.5e9],dest="BOUND")		      
    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args

def plot(t,BBQ,leg,dat,col):
    gmt=pytz.timezone('Europe/Amsterdam');
    pylab.plot_date(t,BBQ,col+'-',label=leg,tz=gmt)
    pylab.xlabel("Time on "+dat.strftime("%Y-%m-%d"))
    pylab.ylabel("BBQ (arbitrary unit)");


if __name__ == "__main__":
    opt,args=parsse(); 
    gmt=pytz.timezone('Europe/Amsterdam');

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
	
    # LHC revolution frequency and period
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    LHCfrev=beta*299792458./26658.883; # the later is the LHC circumference
    LHCTrev=1./LHCfrev;
    # number of turns assumed to be present at eahc line
    nturn=8192;
    
    
    # concatenate datafiles in one single data set
    dataBBQ=[];timeBBQ=[];
    for j,name in enumerate(listname):
        sys.stdout.write('%s: '% name);
        data=parseout(name);
	for v in data['datavars']:
	    if (v.startswith('LHC.BQBBQ.')) and (v.endswith(beam+':ACQ_DATA_'+plane)): var=v;
	# note: in timber data column 0 is always time and column 1 is always the list 
	# of all the other columns
	#print data[var][1][0] # here 1 is column number, 0 is line number
	
	# connect lines between them
	for k,tps in enumerate(data[var][0]):
	    #if (k<=10):print k,dataBBQ[-1]
	    dataBBQ.extend(data[var][1][k]);
	    timebeg=pylab.date2num(tt(data[var][0][k][:-4]))*86400+float(data[var][0][k][-4:]);
	    timeBBQ.extend([timebeg+kk*LHCTrev for kk in range(0,nturn)]);
            #sys.stdout.write('%lf %lf %d\n'% (timebeg+(ind1+1)*LHCTrev,timeBBQ[-1],len(range(ind1+1,nturn))));


    # convert data to arrays
    tp=np.array(timeBBQ);
    BBQ=np.array([float(k) for k in dataBBQ])
    
    # initialize plot
    fig=pylab.figure();
    ax=fig.gca(); 

    # convert time to sth compatible with pylab.plot_date
    tplot=np.array([pylab.num2date(j/86400.) for j in tp]);
    # find data
    dat=pylab.num2date(tp[0]/86400.).date();
    
    # plot raw data
    plot(tplot,BBQ,'Measurement '+beam+' '+plane,dat,'b');
    #file=open('out','w');
    #for i,j in enumerate(BBQ): print >> file, tplot[i],j;
    #file.close();
	    
    
    if (opt.ANA):
	# data analysis
	# take out averagek,j in enumerate(BBQ)
	aver=np.average(BBQ);
	BBQ=BBQ-aver;
	
	# construct "envelop" (i.e. curve of maxima) of the absolute value
	period=60; # number of turns on which we extract maxima
	#BBQenv=np.array([np.max(np.abs(BBQ[k:k+period+1])) for k,j in enumerate(BBQ[:-period-1])]);
	#tenvplot=np.array([tplot[np.argmax(np.abs(BBQ[k:k+period+1]))+k] for k,j in enumerate(BBQ[:-period-1])]);
	## plot envelop
	##plot(tenvplot,BBQenv+aver,beam+plane+' envelop',dat,'r');
	##tenv=pylab.date2num(tenvplot);
	#tenv=np.array([tp[np.argmax(np.abs(BBQ[k:k+period+1]))+k]/86400. for k,j in enumerate(BBQ[:-period-1])]);
	## delete duplicates
	#indnotdup=np.where(np.diff(tenv)!=0);
	#tenv2=tenv[indnotdup];
	#BBQenv2=BBQenv[indnotdup];
	
	BBQenv=[];tenv=[];tenvplot=[];
	for k,j in enumerate(BBQ[:-period-1:period]):
	    BBQenv.append(np.max(np.abs(BBQ[k*period:(k+1)*period+1])));
	    tenv.append(tp[np.argmax(np.abs(BBQ[k*period:(k+1)*period+1]))+k*period]/86400.);
	    tenvplot.append(tplot[np.argmax(np.abs(BBQ[k*period:(k+1)*period+1]))+k*period]);
	BBQenv2=np.array(BBQenv);
	tenv2=np.array(tenv);
	tenvplot2=np.array(tenvplot);
	# plot envelop
	#plot(tenvplot2,BBQenv2+aver,beam+plane+' envelop',dat,'r');
	
			
	# do a "sliding average"
	periodaver=100; # number of points on which we average
	BBQenvaver=np.array([np.average(BBQenv2[k:k+periodaver+1]) for k,j in enumerate(BBQenv2[:-periodaver-1])]);
	tenvaver=np.array([np.average(tenv2[k:k+periodaver+1]) for k,j in enumerate(tenv2[:-periodaver-1])]);
	#plot(pylab.num2date(tenvaver),BBQenvaver,beam+plane+' envelop average',dat,'r');
	
	# find out bounds between which we should fit
	boundmin=opt.BOUND[0]; # fit when signal becomes higher than this
	boundmax=opt.BOUND[1]; # stop fit when we reach this value
	beg=np.where(BBQenvaver>boundmin);beg=beg[0]; # for some reason beg itself is an single-element array of array
	end=np.where(BBQenvaver[beg[0]:-1]>boundmax);
	end=end[0]+beg[0]; # now end[0] is the index looked for
	p=np.polyfit((tenvaver[beg[0]:end[0]]-tenvaver[beg[0]])*86400.,np.log(BBQenvaver[beg[0]:end[0]]),1);
	# note: the fit goes from beg[0] to end[0]-1
	#pylab.figure();
	#pylab.plot((tenvaver[beg[0]:end[0]]-tenvaver[beg[0]])*86400.,np.log(BBQenvaver[beg[0]:end[0]]));
	#pylab.plot((tenvaver[beg[0]:end[0]]-tenvaver[beg[0]])*86400.,p[0]*(tenvaver[beg[0]:end[0]]-tenvaver[beg[0]])*86400.+p[1]);
	tau=1./p[0];a=math.exp(p[1]);
	#print a,tau,p
	sys.stdout.write('Rise time: %lf s\n' % tau);
	#plot(pylab.num2date(tenvaver[beg[0]:end[0]]),BBQenvaver[beg[0]:end[0]],"Envelop average",dat,'g');
	plot(pylab.num2date(tenvaver[beg[0]:end[0]]),aver+a*np.exp((tenvaver[beg[0]:end[0]]-tenvaver[beg[0]])*86400./tau),r'Fit, $\tau$='+str(tau)+' s',dat,'r');

	# the following doesn't work
	#perioddiff=10 # step for the derivative
	#deriv=np.divide(np.diff(BBQenvaver[beg[0]:-1:perioddiff]),np.diff(tenvaver[beg[0]:-1:perioddiff]));
	#pylab.plot(deriv);
	#indices=np.where(np.diff(deriv)<0);indices=indices[0]; # points where derivative decreases


    ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=5))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
    pylab.legend(loc=0);
    pylab.show();

    sys.exit()

