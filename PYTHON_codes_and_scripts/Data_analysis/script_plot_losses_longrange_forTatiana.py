#!/usr/bin/python2.6


# Losses vs. long-range analysis

import pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
import numpy as np
import matplotlib.text as text
from string import split
import glob

def set_fontsize(fig,size):

    # set all fontsize to 'size' in the figure 'fig'
    for o in fig.findobj(text.Text): o.set_fontsize(size);


def init_figure(axes=None):

    pylab.rcParams['xtick.major.pad']='10';
    pylab.rcParams['ytick.major.pad']='10';
    fig=pylab.figure(figsize=(16, 12),facecolor='w', edgecolor=None)
    if (axes!=None): pylab.axes(axes)
    ax=fig.gca()
    
    return fig,ax
    

def end_figure(fig,ax,save=None,legpos=0):

    ax.legend(loc=legpos);
    ax.xaxis.labelpad=12;
    ax.yaxis.labelpad=12;
    set_fontsize(fig,'xx-large');
    
    if (save!=None):
	# save figure (.png and .eps files). "save" contains the root of the file names
	fig.savefig(save+".png")
	fig.savefig(save+".eps",format='eps')

    return;

def plot(x,data,leg,patcol,ylab,ax,logflag,lw=3.,xlab="Number of turns",plotevery=1,colr=None):

    # plotting function. plot every "plotevery" turns.

    if (colr==None):
	if (logflag==0): ax.plot(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5);
	elif (logflag==1): ax.semilogx(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5);
	elif (logflag==2): ax.semilogy(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5);
	elif (logflag==3): ax.loglog(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5);
    else:
	if (logflag==0): ax.plot(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,color=colr);
	elif (logflag==1): ax.semilogx(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,color=colr);
	elif (logflag==2): ax.semilogy(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,color=colr);
	elif (logflag==3): ax.loglog(x[::plotevery],data[::plotevery],patcol,label=leg,linewidth=lw,ms=10.,mew=2.5,color=colr);
	
    ax.set_xlabel(xlab);
    ax.set_ylabel(ylab);
    
def list_files(files):

    # list the files to be read. "files" is a list of names that can contain regular expressions
    # (like *, etc.) from which the list of files to be analysed is constructed
    
    listname=[];
    for f in files:
    	l=glob.glob(f);
	l.sort();
	listname.extend(l);
    #for j in listname: print j;
    
    return listname;
    
if __name__ == "__main__":

    pathmes='';

    # measurements subdirectories (part of the name only)
    fills=['3223','3229','3259','3266','3272','3285','3292','3298','3322','3323','3297','3378'];
    #timelosses=[100,150,258,110,70,360,420,300,43,70,211,40];
    fillgroups=[0,2,5,7,10,11,12]; # how to group fills together
    name='highest_losses_vs_longrange_B1_fill';suffix='_over_instab.txt';
    pat=['xb','or','+g','dk','vm','^c','.y'];

    fig15,ax15=init_figure();
    fig8,ax8=init_figure();
    fig158,ax158=init_figure();
    fig18,ax18=init_figure();

    for igroup in range(1,len(fillgroups)):

	lege='';legen='';

	for idir in range(fillgroups[igroup-1],fillgroups[igroup]):

            lege+=fills[idir];
	    if (idir!=fillgroups[igroup]-1): lege+=' - ';
	    else: legen='fills '+lege;

	    files=list_files([pathmes+name+fills[idir]+suffix]);

            filename=files[0];
            # read data of this file
	    slot=[];losses=[];lrIP1=[];lrIP15=[];lrIP8=[];lrtot=[];
	    # read file
	    for il,line in enumerate(open(filename)):
    		if (il>=1)and(len(split(line))>1):
		    slot.append(int(split(line)[0]))
		    losses.append(float(split(line)[2]))
		    lrIP1.append(int(split(line)[4]))
		    lrIP15.append(int(split(line)[5]))
		    lrIP8.append(int(split(line)[6]))
		    lrtot.append(int(split(line)[7]))


	    plot(lrIP15,np.array(losses),legen,pat[igroup-1],
    		'Losses during instability (10^9 p+/b)',ax15,2,
		xlab='Nb long-range interactions in IP1 and 5');

	    plot(lrIP8,np.array(losses),legen,pat[igroup-1],
    		'Losses during instability (10^9 p+/b)',ax8,2,
		xlab='Nb long-range interactions in IP8');

	    plot(np.array(lrIP15)+np.array(lrIP8),np.array(losses),legen,pat[igroup-1],
    		'Losses during instability (10^9 p+/b)',ax158,2,
		xlab='Nb long-range interactions in IP1, 5 and 8');

	    plot(np.array(lrIP1)+np.array(lrIP8),np.array(losses),legen,pat[igroup-1],
    		'Losses during instability (10^9 p+/b)',ax18,2,
		xlab='Nb long-range interactions in IP1 and 8');


    end_figure(fig15,ax15,legpos=(0.2,0.8));
    end_figure(fig8,ax8,legpos=(0.2,0.8));
    end_figure(fig158,ax158,legpos=(0.2,0.8));
    end_figure(fig18,ax18,legpos=(0.2,0.8));
    pylab.show();
    
