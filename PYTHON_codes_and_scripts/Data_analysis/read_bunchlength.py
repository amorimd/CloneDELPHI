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
from read_Headtail_prt_coupledbunch import plot_bunch
from plot_lib import plot,init_figure,end_figure
from datetime_lib import tt
#from read_FBCT import plot
#import subprocess
import glob


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--average",type=int,nargs=2,
                      help="Specify if we average the bunch lengths over several bunches. 2 arguments: first and last bunch number of the bunches on which to average",
                      metavar="AVER", default=[-1,-1],dest="AVER")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-c", "--bct",
                      help="Specify if we plot BCT data on the same plot and the TIMBER .csv file where to find the data",
                      metavar="BCT", default=None,dest="BCT")		      
    parser.add_option("-d", "--boundstime",nargs=12,
                      help="Specify dates and times between which one computes the bunch-by-bunch average (output to file). Default = no average computed.",
                      type=int,metavar="TIME", default=None,dest="TIME")
    parser.add_option("-f", "--file",
                      help="Specify the TIMBER .csv name of the first file",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-m", "--multi",action="store_true",
                      help="Specify if there are several files saved with the multifile option in Timber",
                      metavar="MULT", default=False,dest="MULT")
    parser.add_option("-n", "--nbfile",type=int,
                      help="Specify maximum number of bunch length data files to take",
                      metavar="NFILE", default=10000,dest="NFILE")
    parser.add_option("-o", "--output",
                      help="Specify an additional comment for the output average file",
                      metavar="OUT", default='',dest="OUT")
    (opt, args) = parser.parse_args()
    #print "Selected File:", opt.FILE
    return opt, args

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
    dataLen={};timeLen={};bunch=[];
    for j,name in enumerate(listname):
        sys.stdout.write('%s: '% name);
        data=parseout(name);
	for v in data['datavars']:
	    if (v.startswith('LHC.BQM.'+beam+':BUNCH_LENGTHS')): var=v;
	# note: in timber data column 0 is always time and column 1 is always the list 
	# of all the other columns
	#print data[var][1][0] # here 1 is column number, 0 is line number
	
	for i,line in enumerate(data[var][1]):
	    for k,l in enumerate(line):
	    	#print k,l
	    	if (l!='0'):
		    if (k not in bunch):
		    	bunch.append(k);
			dataLen[k]=[];
			timeLen[k]=[];
		    dataLen[k].append(l);
		    #print k,dataLen[k][-1]
		    timeLen[k].append(pylab.date2num(tt(data[var][0][i][:-4]))*86400+float(data[var][0][0][-4:]));
	

    # convert data to arrays
    BL=[];tp=[];tplot=[];
    for kbun,bun in enumerate(bunch):
    	#if (kbun==0)or(len(timeLen[bun])==len(tp[0])):
        tp.append(np.array(timeLen[bun]));
    	BL.append(np.array([float(k) for k in dataLen[bun]]));
    	# convert time to sth compatible with pylab.plot_date
    	tplot.append(np.array([pylab.num2date(j/86400.) for j in timeLen[bun]]));
    
    # find date
    dat=pylab.num2date(tp[0][0]/86400.).date();

    # initialize plot
    #fig=pylab.figure(facecolor='w', edgecolor=None);
    #pylab.axes([0.1,0.1,0.71,0.8])
    #ax=fig.gca(); 
    fig,ax=init_figure([0.1,0.1,0.71,0.8]);
    
    # plot raw data
    #col=['r','g','b'];
    for kk,bun in enumerate(bunch):
    	#plot(tplot[kk],BL[kk],beam+', bunch no '+str(bun),col[kk],'Bunch length (seconds)',ax,0,lw=3.,xlab="Time on "+dat.strftime("%Y-%m-%d"));
    	plot(tplot[kk],BL[kk],beam+', bunch no '+str(bun),'','Bunch length (seconds)',ax,0,lw=3.,xlab="Time on "+dat.strftime("%Y-%m-%d"));

    ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=20))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
    timeinterval=np.ceil((tp[0][-1]-tp[0][0])/10.);
    ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
    ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
    #pylab.legend(loc=(1.03,-0.1));
    end_figure(fig,ax,legpos=(1.03,-0.1));
    
    
    # plot vs bunches
    figb,axb=init_figure();
    interval=len(tp[0])/5-1;
    tab=np.arange(0,len(tp[0]),interval)
    for kt,t in enumerate(tp[0][::interval]):
    	tim=pylab.num2date(t/86400.).time().strftime("%H:%M:%S");
	indt=tab[kt];
	try:
	    a=[BL[kb][indt] for kb in range(len(bunch))];
	    plot_bunch(len(bunch),a,tim,'','Bunch length (seconds)','',axb,lw=3.,bunchoffset=1);
    	except IndexError:
	    print "Missing bunch lengths for time=", pylab.num2date(t/86400.).time();
	    #for kb in range(len(bunch)):
	    	#try:
		#    a=BL[kb][indt];
		#except IndexError:
		#    print "Missing bunch length for kb=", kb;
		    
    end_figure(figb,axb);
	

    # plot average, minimum and maximum
    if (opt.AVER[0]>=0):
        fig,ax=init_figure();
	tpav=tp[opt.AVER[0]];BLav=BL[opt.AVER[0]];tot=1;
	tpmax=tp[opt.AVER[0]];BLmax=BL[opt.AVER[0]];
	tpmin=tp[opt.AVER[0]];BLmin=BL[opt.AVER[0]];
	for kk in range(max(0,opt.AVER[0]+1),min(opt.AVER[1]+1,len(bunch))):
	    if kk in bunch:
		tpav=tpav+tp[kk];BLav=BLav+BL[kk];
		tpmax=np.maximum(tpmax,tp[kk]);BLmax=np.maximum(BLmax,BL[kk]);
		tpmin=np.minimum(tpmin,tp[kk]);BLmin=np.minimum(BLmin,BL[kk]);
		tot+=1;
	
	tpav=tpav/tot;BLav=BLav/tot;
	tplotav=np.array([pylab.num2date(j/86400.) for j in tpav])
	tplotmax=np.array([pylab.num2date(j/86400.) for j in tpmax])
	tplotmin=np.array([pylab.num2date(j/86400.) for j in tpmin])
	plot(tplotav,BLav,beam+', average over bunches no '+str(opt.AVER[0]+1)+' to '+str(min(opt.AVER[1]+1,len(bunch))),'b','Bunch length (seconds)',ax,0,lw=3.,xlab="Time on "+dat.strftime("%Y-%m-%d"));
	plot(tplotmax,BLmax,beam+', maximum over bunches no '+str(opt.AVER[0]+1)+' to '+str(min(opt.AVER[1]+1,len(bunch))),'r','Bunch length (seconds)',ax,0,lw=3.,xlab="Time on "+dat.strftime("%Y-%m-%d"));
	plot(tplotmin,BLmin,beam+', minimum over bunches no '+str(opt.AVER[0]+1)+' to '+str(min(opt.AVER[1]+1,len(bunch))),'g','Bunch length (seconds)',ax,0,lw=3.,xlab="Time on "+dat.strftime("%Y-%m-%d"));

	# add BCT data
	if (opt.BCT!=None):
            data2=parseout(opt.BCT);
	    for v in data2['datavars']:
		if (v.endswith(beam+':BEAM_INTENSITY')): var=v;
	    # note: in timber data column 0 is always time and column 1 is always the list 
	    # of all the other columns
	    #print data[var][1][0] # here 1 is column number, 0 is line number

	    ax2 = ax.twinx();
	    dataBCT=np.array([float(j[0]) for j in data2[var][1]]);
	    timeBCT=np.array([pylab.date2num(tt(i[:-4]))*86400+float(data2[var][0][0][-4:]) for i in data2[var][0]]);
	    # find beginning and end tiems of BBQ data
	    beg=np.where(timeBCT>=tp[0][0]);beg=beg[0];
	    if (len(beg)>0): beg=beg[0]-1;
	    else: beg=0;
	    end=np.where(timeBCT>=tp[0][-1]);end=end[0];
	    if (len(end)>0)and(end[0]<len(timeBCT)-1): end=end[0]+1;
	    else: end=-1;
	    # plot
	    plot(timeBCT[beg:end]/86400.,dataBCT[beg:end],"BCT "+beam,dat,'k',"BCT (nb part.)",ax2,lw=3.);
    	    ax2.legend(loc=0);

	# finalize plot
    	ax.xaxis.set_major_locator(matplotlib.dates.SecondLocator(interval=int(timeinterval)))
	ax.xaxis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))
    	pylab.legend(loc=2);
	end_figure(fig,ax);
    

    if (opt.TIME!=None):
    	# compute bunch-by-bunch average    
    
	# initial time (=reference time)
	tim0=datetime(opt.TIME[0],opt.TIME[1],opt.TIME[2],opt.TIME[3],opt.TIME[4],opt.TIME[5]);
	t0=pylab.date2num(tim0);
	# final time
	tim1=datetime(opt.TIME[6],opt.TIME[7],opt.TIME[8],opt.TIME[9],opt.TIME[10],opt.TIME[11]);
	t1=pylab.date2num(tim1);
	print "Computes bunch length average (bunch-by-bunch) between ",tim0.strftime("%Y-%m-%d %H:%M:%S")," and ",tim1.strftime("%Y-%m-%d %H:%M:%S")
    
	file=open('average_bunchlength_'+tim0.strftime("%Y-%m-%d_%H-%M-%S")+'_'+tim1.strftime("%Y-%m-%d_%H-%M-%S")+'_'+beam+opt.OUT+'.txt','w');
	print >> file, "bunch\tfull length (ns)\tsigma (ns)"	
	for k,bnum in enumerate(bunch):

	    begt=np.where(tp[k]/86400.>=t0);
	    if (len(begt[0])==0):
	        begt=0;
		print "Warning: computing average bunch ",bnum," -> first data point after first time needed !";
	    else: begt=begt[0][0]

	    endt=np.where(tp[k]/86400.>=t1);
	    if (len(endt[0])==0):
	    	endt=len(tp[k]);
		print "Warning: computing average bunch ",bnum," -> last data point before last time needed !";
	    else: endt=endt[0][0];

	    blav=np.average(BL[k][begt:endt])
	    sigma=np.sqrt(np.average(BL[k][begt:endt]*BL[k][begt:endt])-blav**2);
	    print  >> file, bnum,"\t",blav*1e9,"\t",sigma*1e9


    pylab.show();

    sys.exit()

