#!/usr/bin/python2.6

import pylab,dateutil,pytz
from datetime import time,datetime,date
from string import split, replace
import matplotlib
import matplotlib.dates


def tt(f1): return datetime.strptime(f1,"%Y-%m-%d %H:%M:%S")


def tb(f1): return time.strftime("%Y-%m-%d %H:%M:%S",time.gmtime(f1))


def tb_(f1):

    # return a string in the format %Y-%m-%d %H:%M:%S corresponding to the pylab datetime (numeric format) f1
    tim=str(pylab.num2date(f1/86400.)); # string of the format 2012-10-09 13:40:00+00:00
    return tim.split('+')[0];


def set_axisdate(axi,axisname,timeinterval,tz=None):
    # set date on the axis 'axisname' ('x' or 'y') of the axes axi, using the timeinterval given.
    timeinterval=int(timeinterval);
    eval('axi.'+axisname+'axis_date(tz=tz)');
    eval('axi.'+axisname+'axis.set_major_locator(matplotlib.dates.SecondLocator(interval=timeinterval))');
    eval('axi.'+axisname+"axis.set_major_formatter(matplotlib.dates.DateFormatter('%H:%M:%S'))");

    return;


def plotdate(time,data,leg,dat,patcol,ylab,ax,logflag,lw=3.,plotevery=1,colr=None):

    # function plotting w.r.t time. plot every "plotevery" turns.
    gmt=pytz.timezone('Europe/Amsterdam');

    if (colr==None):
       ax.plot_date(time[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,tz=gmt);
    else:
       ax.plot_date(time[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,tz=gmt,color=colr);
    
    if (logflag>=2): ax.set_yscale('log');
    if (logflag==1)or(logflag==3): ax.set_xscale('log');
	
    ax.set_xlabel("Time on "+dat.strftime("%Y-%m-%d"));
    ax.set_ylabel(ylab);
    

