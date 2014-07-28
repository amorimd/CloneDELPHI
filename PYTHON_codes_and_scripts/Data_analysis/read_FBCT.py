#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from parser_lib import *
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
from Timber import parseout
import math
import matplotlib
import matplotlib.dates
from plot_lib import init_figure,end_figure,make_movie,plot_save_hist
from datetime_lib import tt
#import subprocess
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--average",type=int,nargs=2,
                      help="Specify if we average the bunch intensity over several bunches. 2 arguments: first and last bucket number (25ns) of the bunches on which to average",
                      metavar="AVER", default=[-1,-1],dest="AVER")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-d", "--boundstime",nargs=12,
                      help="With -k option: specify between which date and time you want to do the movie (year, month, day, hour, minutes, seconds - twice). Default = take the complete time window spanned by the Timber file(s). Without -k option: specify dates and times between which one computes the bunch-by-bunch intensity average (output to file). Default = no average computed. Also used by -p option.",
                      type=int,metavar="TIME", default=None,dest="TIME")
    parser.add_option("-e", "--every",type=int,
                      help="Specify every how many time steps we plot an histogram (default=10)",
                      metavar="EVERY", default=10,dest="EVERY")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--makemovie",action="store_true",
                      help="Specify if we make a movie with the FBCT histogram",
                      metavar="MAKE", default=False,dest="MAKE")
    parser.add_option("-l", "--losses",action="store_true",
                      help="Specify if we switch to the losses plot (instead of intensity plot) for the movie.",
                      metavar="LOSS", default=False,dest="LOSS")
    parser.add_option("-m", "--multi",action="store_true",
                      help="Specify if there are several files saved with the multifile option in Timber",
                      metavar="MULT", default=False,dest="MULT")
    parser.add_option("-n", "--nbfile",type=int,
                      help="Specify maximum number of FBCT data files to take",
                      metavar="NFILE", default=10000,dest="NFILE")
    parser.add_option("-o", "--output",
                      help="Specify an additional comment for the output movie file or output average file",
                      metavar="OUT", default='',dest="OUT")
    parser.add_option("-p", "--plotvslongrange",action="callback",callback=multistring_parse,
                      help="Specify if we plot also the losses vs number of long range, provided in an external file given in argument (11 columns: 25ns-slot number, then 5 unused columns, then number of long-range interaction in resp. IP1, IP2, IP5, IP8, and their total). Relative losses are calculated thanks to times given in option -d. Second argument (optional): threshold level for extraction of slots with highest losses per second (default =2.e6 p+/b/s).",
                      metavar="LONGRANGE", default=None,dest="LONGRANGE")
    parser.add_option("-r", "--plotraw",action="store_false",
                      help="Specify if we do not plot all bunches raw FBCT data (default=plot all of them).",
                      metavar="PLOTRAW", default=True,dest="PLOTRAW")
    parser.add_option("-s", "--sampling",type=int,
                      help="Specify sampling frequency, i.e. we take the FBCT data every ""sampling"" timesteps (default=1, i.e. keep all the data)",
                      metavar="SAMP", default=1,dest="SAMP")
    parser.add_option("-t", "--threshold",type=float,
                      help="Specify the threshold above which the number of particles of the bunch should be in order to be considered (default=0.)",
                      metavar="THR", default=0.,dest="THR")
    parser.add_option("-x", "--xlim",type=int,nargs=2,
                      help="Specify the bunch window on which we zoom the movie (2 arguments)",
                      metavar="XLIM", default=None,dest="XLIM")
    parser.add_option("-y", "--ylim",type=float,nargs=2,
                      help="Specify the intensity window on which we zoom the movie (2 arguments)",
                      metavar="YLIM", default=None,dest="YLIM")
    parser.add_option("-z", "--erase",action="store_true",
                      help="Specify if we erase all _tmp*.png files after making the movie",
                      metavar="RM", default=False,dest="RM")
    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args

def plot(t,data,leg,dat):
    gmt=pytz.timezone('Europe/Amsterdam');
    pylab.plot_date(t,data,label=leg,tz=gmt,fmt='-',lw=3.)
    pylab.xlabel("Time on "+dat.strftime("%Y-%m-%d"))
    pylab.ylabel("Bunch intensity (nb p+)");
    
def plot_save_hist_ref(bunches,FBCT0,FBCT,leg0,leg,ax,fig,i,xlim=None,ylim=None):

    # plot and save an FBCT histogram v. bunches, with the data in FBCT0 in red
    # and the data in FBCT  in black (FBCT0 is the reference)
    ax.cla()
    ax.cla()
    ax.bar(bunches,FBCT0,facecolor='r',label=leg0,edgecolor='r')
    ax.bar(bunches,FBCT,facecolor='k',label=leg,edgecolor='k')
    ax.set_xlabel("25ns-slot number");
    ax.set_ylabel("Bunch intensity (nb p+)");
    if (xlim!=None): ax.set_xlim(xlim);
    if (ylim!=None): ax.set_ylim(ylim);
    fname = '_tmp%d'%i
    end_figure(fig,ax,save=fname,legpos='lower right');

def read_longrange_file(filename):

    # read Xavier's kind of long-range file (11 columns)
    
    slot=[];lrIP1=[];lrIP2=[];lrIP5=[];lrIP8=[];lrtot=[];
    
    # read file
    for il,line in enumerate(open(filename)):
    	if (len(split(line))>1):
	    slot.append(int(split(line)[0]))
	    lrIP1.append(int(split(line)[6]))
	    lrIP2.append(int(split(line)[7]))
	    lrIP5.append(int(split(line)[8]))
	    lrIP8.append(int(split(line)[9]))
	    lrtot.append(int(split(line)[10]))
	
    return slot,np.array(lrIP1),np.array(lrIP2),np.array(lrIP5),np.array(lrIP8),np.array(lrtot);


def extract_beg_end(opttime,tp,optevery=0):

    # extract indices of table of times tp (date2num format*86400) such that dates-times provided 
    # in opttime (2 of them, expressed each as 6 integers, e.g. 2012,5,30,11,34,32) corresponds to 
    # tp[begt] and tp[endt+optevery]
    
    # initial time
    tim0=datetime(opttime[0],opttime[1],opttime[2],opttime[3],opttime[4],opttime[5]);
    t0=pylab.date2num(tim0);
    begt=np.where(tp/86400.>=t0);
    if (len(begt[0])==0):
	begt=0;
	print "Warning: -> first data point after first time needed !";
    else: begt=begt[0][0];

    # final time
    tim1=datetime(opttime[6],opttime[7],opttime[8],opttime[9],opttime[10],opttime[11]);
    t1=pylab.date2num(tim1);
    endt=np.where(tp/86400.>=t1);
    if (len(endt[0])==0):
        endt=len(tp)-optevery;
	print "Warning -> last data point before last time needed !";
    else: endt=endt[0][0];

    return begt,endt;


if __name__ == "__main__":
    opt,args=parsse(); 
    gmt=pytz.timezone('Europe/Amsterdam');

    if opt.BEAM=="1":
    	beam="B1";
    elif opt.BEAM=="2":
    	beam="B2";
    else: print "specify beam 1 or 2"; sys.exit()


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
    dataFBCT={};timeFBCT={};bunch=[];
    for j,name in enumerate(listname):
        sys.stdout.write('%s: '% name);
        data=parseout(name);
	for v in data['datavars']:
	    if (v.startswith('LHC.BCTFR.'))and(v.endswith(beam+':BUNCH_INTENSITY')): var=v;
	# note: in timber data column 0 is always time and column 1 is always the list 
	# of all the other columns
	#print data[var][1][0] # here 1 is column number, 0 is line number
	
	for i,line in enumerate(data[var][1][::opt.SAMP]):
	    for k,l in enumerate(line):
	    	#print k,l
	    	if (float(l)>opt.THR):
		    if (k not in bunch):
		    	bunch.append(k);
			dataFBCT[k]=[];
			timeFBCT[k]=[];
		    dataFBCT[k].append(float(l));
		    #print k,dataFBCT[k][-1]
		    timeFBCT[k].append(pylab.date2num(tt(data[var][0][i][:-4]))*86400+float(data[var][0][0][-4:]));
	

    print "Concatenation finished";
    
    # convert data to arrays
    FBCT=[];tp=[];tplot=[];
    bunches=np.array(bunch);bunches.sort();
    for bun in bunches:
        tp.append(np.array(timeFBCT[bun]));
    	FBCT.append(np.array([k for k in dataFBCT[bun]]));
    	# convert time to sth compatible with pylab.plot_date
    	tplot.append(np.array([pylab.num2date(j/86400.) for j in timeFBCT[bun]]));
    
    print "Conversion to array finished";
    
    # find date
    dat=pylab.num2date(tp[0][0]/86400.).date();

    # initialize plot
    #fig=pylab.figure();
    #pylab.axes([0.1,0.1,0.71,0.8])
    #ax=fig.gca();
    fig,ax=init_figure([0.1,0.1,0.71,0.8]);
    
    if (opt.PLOTRAW):
	# plot raw data
	for kk,bun in enumerate(bunches):
    	    plot(tplot[kk],FBCT[kk],beam+', bunch no '+str(bun),dat);

	print "25ns-bucket no: ",bunches
	print "number of bunches: ",len(bunches)

	timeinterval=np.ceil((timeFBCT[bunches[0]][-1]-timeFBCT[bunches[0]][0])/8.);
	print "Time interval for plot: ", timeinterval
	ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
	#pylab.legend(loc=(1.03,-0.1));
	end_figure(fig,ax,legpos=(1.03,-0.1));
    
    if (opt.AVER[0]>=0):
    	#fig=pylab.figure();
	#ax=fig.gca();
    	fig,ax=init_figure();
	beg=np.where(bunches>=opt.AVER[0]);beg=beg[0][0]; # for some reason beg itself is a single-element array of array
	end=np.where(bunches>=opt.AVER[1]);
	if (len(end[0])==0): end=len(bunches)-1;
	else: end=end[0][0];
	if (bunches[end]>opt.AVER[1]): end=end-1;
	#print beg,end
	tpav=tp[beg];FBCTav=FBCT[beg];tot=1;
	tpmax=tp[beg];FBCTmax=FBCT[beg];
	tpmin=tp[beg];FBCTmin=FBCT[beg];
	for kk in range(beg+1,end+1):
	    #print bunch[kk]
	    tpav=tpav+tp[kk];FBCTav=FBCTav+FBCT[kk];
	    tpmax=np.maximum(tpmax,tp[kk]);FBCTmax=np.maximum(FBCTmax,FBCT[kk]);
	    tpmin=np.minimum(tpmin,tp[kk]);FBCTmin=np.minimum(FBCTmin,FBCT[kk]);
	    tot+=1;
	
	tpav=tpav/tot;FBCTav=FBCTav/tot;
	tplotav=np.array([pylab.num2date(j/86400.) for j in tpav])
	tplotmax=np.array([pylab.num2date(j/86400.) for j in tpmax])
	tplotmin=np.array([pylab.num2date(j/86400.) for j in tpmin])
	plot(tplotav,FBCTav,beam+', average over bunches no '+str(bunches[beg])+' to '+str(bunches[end]),dat);
	plot(tplotmax,FBCTmax,beam+', maximum over bunches no '+str(bunches[beg])+' to '+str(bunches[end]),dat);
	plot(tplotmin,FBCTmin,beam+', minimum over bunches no '+str(bunches[beg])+' to '+str(bunches[end]),dat);
	#ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=20))
	
	timeinterval=np.ceil((tpav[-1]-tpav[0])/8.);
    	ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
    
    	end_figure(fig,ax);
    
    pylab.show();
    
    
    if (opt.MAKE):
    	# make a movie with the FBCT histogram
	# initialize histogram plot
	fighist,axhist=init_figure();
	if (opt.TIME==None):
	    # reference time
	    begt=0;endt=len(timeFBCT[bunches[0]])-opt.EVERY;
	else: 
	    begt,endt=extract_beg_end(opt.TIME,tp[0],optevery=opt.EVERY)

	# reference FBCT
	FBCT0=[np.average(dataFBCT[k][begt:begt+opt.EVERY]) for k in bunches];
	tim0=pylab.num2date(np.average(timeFBCT[bunches[0]][begt:begt+opt.EVERY])/86400.);
	
	if opt.LOSS:
	    if (opt.YLIM==None):
	    	optylim=[];
	    	optylim.append(0.);
	    	# calculate maximum losses (for scale -> it should not move)
	    	optylim.append(1.2*np.max(np.array([(dataFBCT[k][begt]-dataFBCT[k][endt]) for k in bunches])));
	    else: optylim=opt.YLIM;
	
	#print tim0.strftime("%Y-%m-%d %H:%M:%S")
	for i in range(begt+1,endt,opt.EVERY):
	    tim=pylab.num2date(np.average(timeFBCT[bunches[0]][i:i+opt.EVERY])/86400.);
	    #print i,tim.strftime("%Y-%m-%d %H:%M:%S")
	    #for k in bunches: print k,dataFBCT[k][i];
	    if opt.LOSS:
	    	FBCT1=[(FBCT0[ik]-np.average(dataFBCT[k][i:i+opt.EVERY])) for ik,k in enumerate(bunches)];
	    	plot_save_hist(bunches,FBCT1,beam+', '+tim.strftime("%Y-%m-%d %H:%M:%S"),axhist,fighist,
			i,xlim=opt.XLIM,ylim=optylim,tit='ref. intensities at '+tim0.strftime("%Y-%m-%d at %H:%M:%S"),
			ylab="Bunch losses (relative to the ref.)",legposition='upper left');
	    else:
	    	FBCT1=[np.average(dataFBCT[k][i:i+opt.EVERY]) for k in bunches];
	    	plot_save_hist_ref(bunches,FBCT0,FBCT1,beam+', '+tim0.strftime("%Y-%m-%d %H:%M:%S"),
			beam+', '+tim.strftime("%Y-%m-%d %H:%M:%S"),axhist,fighist,i,xlim=opt.XLIM,ylim=opt.YLIM);
	    
	os.system('cp _tmp%d'%i+'.png losses_'+beam+opt.OUT+'_'+tim.strftime("%Y%m%d_%H%M%S")+'.png');
	os.system('cp _tmp%d'%i+'.eps losses_'+beam+opt.OUT+'_'+tim.strftime("%Y%m%d_%H%M%S")+'.eps');
	make_movie(listname[0]+'_'+beam+opt.OUT+'.gif','_tmp',flagrm=opt.RM);
	pylab.close(fighist)
	
    elif (opt.TIME!=None):
    
	# initial time (=reference time)
	tim0=datetime(opt.TIME[0],opt.TIME[1],opt.TIME[2],opt.TIME[3],opt.TIME[4],opt.TIME[5]);
	# final time
	tim1=datetime(opt.TIME[6],opt.TIME[7],opt.TIME[8],opt.TIME[9],opt.TIME[10],opt.TIME[11]);
	print "Computes bunch-by-bunch intensity average between ",tim0.strftime("%Y-%m-%d %H:%M:%S")," and ",tim1.strftime("%Y-%m-%d %H:%M:%S")
	
	file=open('average_FBCT_'+tim0.strftime("%Y-%m-%d_%H-%M-%S")+'_'+tim1.strftime("%Y-%m-%d_%H-%M-%S")+'_'+beam+opt.OUT+'.txt','w');
	print >> file, "bunch\tintensity (10^11 p+)\tsigma (10^11 p+)"	
	for k,bnum in enumerate(bunches):

	    begt,endt=extract_beg_end(opt.TIME,tp[k]);
	    intav=np.average(FBCT[k][begt:endt])
	    sigma=np.sqrt(np.average(FBCT[k][begt:endt]*FBCT[k][begt:endt])-intav**2);
	    print  >> file, bnum,"\t",intav/1e11,"\t",sigma/1e11
	    
	file.close();
    
    
    if (opt.LONGRANGE!=None):
    
        slot,lrIP1,lrIP2,lrIP5,lrIP8,lrtot=read_longrange_file(opt.LONGRANGE[0]);
	lrIP15=lrIP1+lrIP5; # sum of IP1 and IP5 numbers of long-range interaction
	#lrtot=lrtot-lrIP2;
	
    	figlr,axlr=init_figure([0.15,0.1,0.8,0.7]);
    	#figlrtot,axlrtot=init_figure();

	# initial time (=reference time)
	tim0=datetime(opt.TIME[0],opt.TIME[1],opt.TIME[2],opt.TIME[3],opt.TIME[4],opt.TIME[5]);
	# final time
	tim1=datetime(opt.TIME[6],opt.TIME[7],opt.TIME[8],opt.TIME[9],opt.TIME[10],opt.TIME[11]);
	
	# file for slots with highest losses
	file=open('highest_losses_vs_longrange_'+beam+opt.OUT+'.txt','w');
	print >> file, "slot (25ns)\tinitial intensity (10^11 p+)\tlosses (10^9 p+/b)\tlosses per unit time (10^6 p+/b/s)\tnb longe-range (IP1)\tnb longe-range (IP1+IP5)\tnb long-range (IP8)\ttotal nb long-range"
	if (len(opt.LONGRANGE)>1): threshold_losses=float(opt.LONGRANGE[1]);
	else: threshold_losses=2.e6;

	for ibnum,bnum in enumerate(bunches):
	
	    if (opt.TIME==None):
		# reference time
		begt=0;endt=len(timeFBCT[bnum]);
	    else: 
		# initial time (=reference time)
		begt,endt=extract_beg_end(opt.TIME,tp[ibnum]);
		
	    if (bnum in slot):
	    
	        kbnum=slot.index(bnum);
		ref_int=dataFBCT[bnum][begt];
		loss=ref_int-dataFBCT[bnum][endt];
		losstime=loss/(tp[ibnum][endt]-tp[ibnum][begt]);
		
		if (kbnum==0):
		    axlr.plot(lrIP15[kbnum],losstime,'xb',label='Nb long-range in IP1 + IP5',lw=3.,ms=10.,mew=2.5);
		    axlr.plot(lrIP8[kbnum],losstime,'or',label='Nb long-range in IP8',lw=3.,ms=10.,mew=2.5);
		    axlr.plot(lrtot[kbnum],losstime,'+g',label='Total nb of long-range',lw=3.,ms=10.,mew=2.5);
		else:
		    axlr.plot(lrIP15[kbnum],losstime,'xb',lw=3.,ms=10.,mew=2.5);
		    axlr.plot(lrIP8[kbnum],losstime,'or',lw=3.,ms=10.,mew=2.5);
		    axlr.plot(lrtot[kbnum],losstime,'+g',lw=3.,ms=10.,mew=2.5);

		# save slots with highest losses
		if (losstime>threshold_losses): print  >> file, bnum,"\t",ref_int/1e11,"\t",loss/1e9,"\t",losstime/1e6,"\t",lrIP1[kbnum],"\t",lrIP15[kbnum],"\t",lrIP8[kbnum],"\t",lrtot[kbnum];
		

	axlr.set_title(beam+' relative losses at '+tim1.strftime("%H:%M:%S")+' w.r.t '+tim0.strftime("%H:%M:%S"));
	axlr.set_xlabel('Number of long-range interactions');
	axlr.set_ylabel('Losses per unit time (p+/bunch/s)');
	#axlrtot.set_title(beam+' relative losses at '+tim1.strftime("%H:%M:%S")+' w.r.t '+tim0.strftime("%H:%M:%S"));
	#axlrtot.set_xlabel('Number of long-range interactions');
	#axlrtot.set_ylabel('Losses (p+/bunch)');
     	end_figure(figlr,axlr,save='losses_vs_longrange_'+beam+opt.OUT,legpos=(0.7,1.1));
     	#end_figure(figlrtot,axlrtot);

	file.close();
   
    sys.exit()

