#!/usr/bin/python2.6


# Losses vs. long-range analysis

import pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
import numpy as np
from plot_lib import init_figure,end_figure,plot
from io_lib import list_files
from string import split

if __name__ == "__main__":

    pathmes='/home/nmounet/LHCTimberData/';

    dirgen=pathmes+"*/Fill_";
    FBCTname="*/FBCT*/";

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

	    files=list_files([dirgen+fills[idir]+FBCTname+name+fills[idir]+suffix]);

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


    end_figure(fig15,ax15,save=pathmes+'losses_vs_longrangeIP15_summary',legpos=(0.2,0.8));
    end_figure(fig8,ax8,save=pathmes+'losses_vs_longrangeIP8_summary',legpos=(0.2,0.8));
    end_figure(fig158,ax158,save=pathmes+'losses_vs_longrangeIP158_summary',legpos=(0.2,0.8));
    end_figure(fig18,ax18,save=pathmes+'losses_vs_longrangeIP18_summary',legpos=(0.2,0.8));
    #pylab.show();
    
