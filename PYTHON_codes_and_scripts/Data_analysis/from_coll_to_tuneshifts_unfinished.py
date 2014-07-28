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
from plot_lib import cmap
from datetime_lib import tt
import matplotlib
import matplotlib.dates


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-o", "--output",help="Specify output filename",
                      default="coll.txt",dest="OUT")
    parser.add_option("-p", "--plot",
                      help="Specify a string: all the collimators containing this string will be plotted vs. time",
                      metavar="PLOT", default=None,dest="PLOT")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv file name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-r", "--ref",
                      help="Specify the reference file name (where coll. names, beta function and so on can be found)",
                      metavar="RFILE", default=None,dest="RFILE")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    print "Selected Reference file:", opt.RFILE
    return opt, args


def plot(t,hg,leg,dat):
    gmt=pytz.timezone('Europe/Amsterdam');
    pylab.plot_date(t,hg,fmt='-',label=leg,tz=gmt)
    pylab.xlabel("Time on "+dat.strftime("%Y-%m-%d"))
    pylab.ylabel("Half gap (mm)");

                  
if __name__ == "__main__":
    opt,args=parsse(); data=parseout(opt.FILE);

    gmt=pytz.timezone('Europe/Amsterdam');

    if opt.BEAM=="1":
    	beam="B1";tcsgrefname="TCSG.4R6.B1";
    elif opt.BEAM=="2":
    	beam="B2";tcsgrefname="TCSG.4L6.B2";
    else: print "specify beam 1 or 2"; sys.exit()
    
    # read reference file
    collname=[];
    fid=open(opt.RFILE,"r")
    for j,l in enumerate(fid.readlines()):
	ll=l.strip().split();
    	if (j==0):
	    # find column names "name", "angle", "betax", "betay", "Material" and "Length"
	    for k,col in enumerate(ll):
	    	if col.startswith('name'): namecol=k;
		if col.startswith('angle'): anglecol=k;
		if col.startswith('betax'): betaxcol=k;
		if col.startswith('betay'): betaycol=k;
		if col.startswith('Material'): materialcol=k;
		if col.startswith('Length'): lengthcol=k;
	    # if first column is '#', do not take it into account in numbering
	    if ll[0].startswith('#'): 
	    	namecol=namecol-1;anglecol=anglecol-1;
		betaxcol=betaxcol-1;betaycol=betaycol-1;
		materialcol=materialcol-1;lengthcol=lengthcol-1;
    	else:
	    name.append(ll[namecol].replace(".B1","").replace(".B",""));
	    angle.append(ll[anglecol]);
	    betax.append(ll[betaxcol]);
	    betay.append(ll[betaycol]);
	    mat.append(ll[materialcol]);
	    length.append(ll[lengthcol]);
    fid.close()
 
    # open output file
    file=open(opt.OUT,'w');

    if (opt.PLOT!=None): 
    	fig=pylab.figure();
    	pylab.axes([0.1,0.1,0.71,0.8])
    	ax=fig.gca(); 

    
    for k,name in enumerate(collname):
    
        var=[];
	for v in data['datavars']:
	    if v.startswith(name.replace(".B1","")): var.append(v);
	
	# note: for a single collimator the measurement
	# time is assumed to be the same for GD, GU, etc.
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
	
	
	# tune shift from this collimator (at each time)
	
	
	tplot=np.array([pylab.num2date(j/86400.) for j in tp]);
	dat=pylab.num2date(tp[0]/86400.).date();
	#for j in tplot: sys.stdout.write('toto --> %s\n'% j)
	#sys.stdout.write('toto --> %s\n'% dat.ctime())
        #sys.stdout.flush();
	if (opt.PLOT!=None): 
	    if (name.find(opt.PLOT)!=-1): plot(tplot,halfgap,name,dat);
	    
	#print name,dt,hg,halfgap,gu[0]
	
	
    file.close()
	
    if (opt.PLOT!=None): 
    	ax.xaxis.set_major_locator(matplotlib.dates.MinuteLocator(interval=10))
    	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M'))
    	pylab.legend(loc=(1.03,-0.1));
	pylab.show();

	
    sys.exit()

