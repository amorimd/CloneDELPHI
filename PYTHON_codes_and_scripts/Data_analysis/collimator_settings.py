#!/usr/bin/python

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
from plot_lib import init_figure,end_figure,build_colors
from datetime_lib import tt
import matplotlib
import matplotlib.dates
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--emittance",type=float,
                      help="Specify the normalized emittance (in mm.mrad) in x and y. Default=3.5 mm.mrad.",
                      metavar="EMIT", default=3.5,dest="EMIT")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify the energy in GeV. Default= 450 GeV",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-i", "--interval",type=int,
                      help="Specify the time interval (in seconds) for the x-axis ticks in the plot",
                      metavar="INTER", default=600,dest="INTER")
    parser.add_option("-m", "--multi",action="store_true",
                      help="Specify if there are several files saved with the multifile option in Timber",
                      metavar="MULT", default=False,dest="MULT")
    parser.add_option("-n", "--nbfile",type=int,
                      help="Specify maximum number of data files to take",
                      metavar="NFILE", default=100,dest="NFILE")
    parser.add_option("-o", "--output",help="Specify output filename",
                      default="coll.txt",dest="OUT")
    parser.add_option("-p", "--plot",
                      help="Specify a string: all the collimators containing this string will be plotted vs. time",
                      metavar="PLOT", default=None,dest="PLOT")
    parser.add_option("-r", "--ref",
                      help="Specify the reference file name (where coll. names can be found)",
                      metavar="RFILE", default=None,dest="RFILE")
    parser.add_option("-t", "--time",nargs=6,
                      help="Specify at which date and time you want to output the settings (year, month, day, hour, minutes, seconds)",
                      type=int,metavar="TIME", default=None,dest="TIME")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    print "Selected Reference file:", opt.RFILE
    return opt, args


def plot(t,hg,leg,dat,colr):
    gmt=pytz.timezone('Europe/Amsterdam');
    pylab.plot_date(t,hg,fmt='-',label=leg,tz=gmt,lw=3.,color=colr)
    pylab.xlabel("Time on "+dat.strftime("%Y-%m-%d"))
    pylab.ylabel("Half gap (mm)");


#def cmap(kys):
#    maps=[]
#    for u in range(kys):
#        maps.append((round(random.uniform(0,1),1),
#                     round(random.uniform(0,1),1),
#                     round(random.uniform(0,1),1)))    
#    return maps
                  
		  
def concatenate_data_Timber(listname):

    # concatenate datafiles in one single data set
    for j,name in enumerate(listname):
        sys.stdout.write('%s: '% name);
        data1=parseout(name);
	if (j>0):
	    for v in data1['datavars']:
	    	if (v in data['datavars']):
	    	    data[v][0].extend(data1[v][0]);
	    	    data[v][1].extend(data1[v][1]);
		else:
		    data['datavars'].append(v);
	    	    data[v]=data1[v];
		    
	else: data=data1;

    return data;


def compute_onesigma_halfgap(betax,betay,skewangle,gamma,emit):

    # compute a half-gap of one sigma for a general skew collimator, given
    # beta functions, skew angle, relativistic factor and emittance
    
    beta=np.sqrt(1.-1./(gamma*gamma));
    
    return np.sqrt(betax*emit*np.cos(skewangle)**2./(beta*gamma) + 
                    betay*emit*np.sin(skewangle)**2./(beta*gamma))
		    

def read_ref_coll_file(filename):

    # read reference collimator settings file
    collname=[];angle=[];betax=[];betay=[];
    fid=open(filename,"r")
    
    for j,l in enumerate(fid.readlines()):
	ll=l.strip().split();
	#sys.stdout.write(l);sys.stdout.flush();
    	if (j==0):
	    # find column name "name"
	    for k,col in enumerate(ll):
	    	if col.startswith('name'): namecol=k;
	    	if (col.startswith('Angle'))or(col.startswith('angle')): anglecol=k;
	    	if col.startswith('betax'): betaxcol=k;
	    	if col.startswith('betay'): betaycol=k;
	    # if first column is '#', do not take it into account in numbering
	    if ll[0].startswith('#'):
	    	namecol=namecol-1;
		anglecol-=1;betaxcol-=1;betaycol-=1;
    	else:
	    collname.append(ll[namecol]);
	    angle.append(float(ll[anglecol]));
	    betax.append(float(ll[betaxcol]));
	    betay.append(float(ll[betaycol]));
    fid.close()
    
    return collname,angle,betax,betay;
		    

if __name__ == "__main__":
    opt,args=parsse();

    gmt=pytz.timezone('Europe/Amsterdam');

    if opt.BEAM=="1":
    	beam="B1";tcsgrefname="TCSG.4R6.B1";
    elif opt.BEAM=="2":
    	beam="B2";tcsgrefname="TCSG.4L6.B2";
    else: print "Specify beam 1 or 2"; sys.exit()
    
    if (opt.RFILE==None): print "Specify reference file for the collimator settings!"; sys.exit();

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
	

    # concatenate datafiles in one single data set
    data=concatenate_data_Timber(listname);

    # read reference file
    collname,angle,betax,betay=read_ref_coll_file(opt.RFILE);
 
    
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    
    # open output file
    file=open(opt.OUT,'w');
    print >> file, "name\thalfgap[m]\tnsigma"

    if (opt.PLOT!=None): 
    	#fig=pylab.figure();
    	#pylab.axes([0.1,0.1,0.71,0.8])
    	#ax=fig.gca();
	fig,ax=init_figure(axes=[0.05,0.1,0.76,0.8]);
	ncol=0;
	for name in collname:
	    if (name.find(opt.PLOT)!=-1): ncol+=1;
	col=build_colors(ncol);

    
    kplot=0;
    for k,name in enumerate(collname):
    
        var=[];
	for v in data['datavars']:
	    if v.startswith(name.replace("."+beam,"")): var.append(v);
	
	# note: for a single collimator the measurement
	# time is assumed to be the same for GD, GU, etc.
	# (use advanced option of TIMBER in case -> Time scaled at fixed interval
	# and Interpolate)
	if (not name.startswith('TCDQA')):
	    for v in var:
	    	if v.endswith('GU'):
		    gu=np.array([float(j[0]) for j in data[v][1]]);tps=data[v][0];
	    	if v.endswith('GD'): 
		    #for j in data[v][1]: sys.stdout.write('toto --> %s'% j)
        	    #sys.stdout.flush();
	    	    gd=np.array([float(j[0]) for j in data[v][1]]);
	
	    halfgap=(gu+gd)/4.;
	    # time conversion
	    tp=np.array([pylab.date2num(tt(j[:-4])) for j in tps])*86400.+np.array([float(j[-4:]) for j in tps]);
	    #for j in tp: sys.stdout.write('toto --> %s\n'% j)
            #sys.stdout.flush();
	
	else:
	    
	    tcsgref=[];
	    for v in data['datavars']:
	    	if v.startswith(tcsgrefname): tcsgref.append(v);
	    
	    for v in var:
	    	if v.endswith('LU'): lu=np.array([float(j[0]) for j in data[v][1]]);tps=data[v][0];
	    	if v.endswith('LD'): ld=np.array([float(j[0]) for j in data[v][1]]);
	    
	    for v in tcsgref:
	    	if v.endswith('LU'): luref=np.array([float(j[0]) for j in data[v][1]]);tpsref=data[v][0];
	    	if v.endswith('LD'): ldref=np.array([float(j[0]) for j in data[v][1]]);
	    	if v.endswith('RU'): ruref=np.array([float(j[0]) for j in data[v][1]]);
	    	if v.endswith('RD'): rdref=np.array([float(j[0]) for j in data[v][1]]);
	    
	    # orbit offset at this position (from ref. collimator)
	    offsetref=(luref+ldref+ruref+rdref)/4.;
	    
	    # need to interpolate this offset at the times given for the TCDQ
	    # first time conversion
	    tp=np.array([pylab.date2num(tt(j[:-4])) for j in tps])*86400.+np.array([float(j[-4:]) for j in tps]);
	    tpr=np.array([pylab.date2num(tt(j[:-4])) for j in tpsref])*86400.+np.array([float(j[-4:]) for j in tpsref]);	    
	    offinterp=np.interp(tp,tpr,offsetref);
	    
	    halfgap=(lu+ld)/2.-offinterp;
	    # note: only 1 TCDQ in Timber; for the second one, lu and ld will actually be kept at the same value -> same half gap
	
	
	#tplot=np.array([pylab.num2date(j/86400.) for j in tp]);
	#print name,tp[0]
	dat=pylab.num2date(tp[0]/86400.).date();
	#for j in tplot: sys.stdout.write('toto --> %s\n'% j)
	#sys.stdout.write('toto --> %s\n'% dat.ctime())
        #sys.stdout.flush();
	if (opt.PLOT!=None):
	    if (name.find(opt.PLOT)!=-1):
	    	plot(tp/86400.,halfgap,name,dat,col[kplot]);
		kplot+=1;
	    
	if (opt.TIME!=None):
	    dt=pylab.date2num(datetime(opt.TIME[0],opt.TIME[1],opt.TIME[2],opt.TIME[3],opt.TIME[4],opt.TIME[5]));
	    #print dt,tp[0]/86400
	    if (k==0): print "Time at which data is taken: ",pylab.num2date(dt)
	else:
	    # by default takes first time
	    dt=pylab.date2num(pylab.num2date(tp[0]/86400.));
	    if (k==0): print "Time at which data is taken: ",pylab.num2date(dt)
	    
	hg=np.interp(dt*86400,tp,halfgap)*1.e-3;
	gap1sig=compute_onesigma_halfgap(betax[k],betay[k],angle[k],gamma,opt.EMIT*1.e-6)
	nsig=hg/gap1sig;
	print >> file, name, "\t", hg, "\t", nsig;
	#print name,dt,hg,halfgap,gu[0]
	
	
    file.close()
	
    if (opt.PLOT!=None): 
	if (opt.INTER<60):
    	    ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=opt.INTER))
    	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
	elif (opt.INTER<3600):
    	    ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=opt.INTER/60))
    	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
	else :
    	    ax.xaxis.set_major_locator(matplotlib.dates.HourLocator(interval=opt.INTER/3600))
    	    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
    	#pylab.legend(loc=(1.03,-0.1));
	end_figure(fig,ax,legpos=(1.03,-0.1));
	pylab.show();

	
    sys.exit()

