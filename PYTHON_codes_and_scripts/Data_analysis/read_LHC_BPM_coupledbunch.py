#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from optparse import OptionParser
from string import *
import pylab,os,re
from SussixNM import *
from parser_lib import bunch_parse
from read_cfg import read_cfg
from plot_lib import init_figure,plot,plot2D
from read_Headtail_prt_sussix import extract
from read_Headtail_prt_coupledbunch import gettune,find_max,plot_bunch,plot_spec,fft1d_bunch,fft1d_turn,fft2d,svd
from numpy.fft import fft,fft2
import numpy as np
from sddsdata_riccardodemaria import *


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--fftanalysis",action="store_true",
                      help="Specify if we perform 1D FFT analysis (both in bunches and in turns)",
                      metavar="FFT1", default=False,dest="FFT1")
    parser.add_option("-b", "--bunchoffset",type=int,
                      help="Specify the the offset to add to the bunch number",
                      metavar="BOFF", default=0,dest="BOFF")
    parser.add_option("-c", "--circum",type=float,
                      help="Specify the circumference in m (default=LHC)",
                      metavar="CIRC", default=26658.883,dest="CIRC")
    parser.add_option("-d", "--fft2danalysis",action="store_true",
                      help="Specify if we perform 2D FFT analysis",
                      metavar="FFT2", default=False,dest="FFT2")
    parser.add_option("-e", "--energy",type=float,
                      help="Specify energy in GeV (default=450 GeV)",
                      metavar="EN", default=450.,dest="EN")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot (several files possible)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-o", "--output",help="Specify output suffix for mode number and tune of highest mode",
                      default="mode2D.txt",dest="OUT")
    #parser.add_option("-t", "--turns",type=int,nargs=2,
    #                  help="Specify the first and last turn (+1) between which we plot and analyse",
    #                  metavar="TURNS", default=[0,10000000],dest="TURNS")
    parser.add_option("-q", "--tunes",type=float,nargs=2,
                      help="Specify the factional part of the tunes (x and y)",
                      metavar="TUNES", default=[0.28,0.31],dest="TUNES")
    parser.add_option("-r", "--raw",action="store_true",
                      help="Specify if we also plot raw data (1D for the last turn and last bunch, and 2D)",
                      metavar="RAW", default=False,dest="RAW")
    parser.add_option("-s", "--svd",action="store_true",
                      help="Specify if we perform Singular Value Decomposition (time domain analysis)",
                      metavar="SVD", default=False,dest="SVD")
    (opt, args) = parser.parse_args()
    return opt, args

    	    
def read_BPM(filename):

    # read BPM file
    bpm=sddsdata(filename)

    nturns=bpm.data[0]['nbOfCapTurns'][0]
    nbunch=bpm.data[0]['nbOfCapBunches'][0]
    nbpm=bpm.data[0]['bpmNames'].shape[0]

    x=bpm.data[0]['horPositionsConcentratedAndSorted']
    y=bpm.data[0]['verPositionsConcentratedAndSorted']
    
    # find bpm(s) with non zero data
    ix=[];iy=[];
    for i in range(0,nbpm):
        if (np.max(np.abs(x[i*nturns*nbunch:(i+1)*nturns*nbunch]))>0): ix.append(i)
        if (np.max(np.abs(y[i*nturns*nbunch:(i+1)*nturns*nbunch]))>0): iy.append(i)

    print filename,": ",len(ix)," BPMs in x, ",len(iy)," BPMs in y"
         
    for ii,i in enumerate(ix):
    	x1=x[i*nturns*nbunch:(i+1)*nturns*nbunch]
	# for each bunch, take out the average over all turns
	for j in range(nbunch):
	    x2=x1[j*nturns:(j+1)*nturns]
	    x1[j*nturns:(j+1)*nturns]=x2-np.average(x2)
	# make an average over all bpms
        if (ii==0): xave=x1
	else: xave=x1+xave
										    	
    for ii,i in enumerate(iy):
    	y1=y[i*nturns*nbunch:(i+1)*nturns*nbunch]
	# for each bunch, take out the average over all turns
	for j in range(nbunch):
	    y2=y1[j*nturns:(j+1)*nturns]
	    y1[j*nturns:(j+1)*nturns]=y2-np.average(y2)
	# make an average over all bpms
        if (ii==0): yave=y1
	else: yave=y1+yave
	
    xave/=float(len(ix));
    yave/=float(len(iy));
    											 
    return nturns,nbunch,xave,yave,bpm


if __name__ == "__main__":


    opt,args=parsse();
    
    print 'Energy: ',opt.EN,' GeV';
    gamma=opt.EN/0.938272 # the later is the rest mass of a proton
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./opt.CIRC; # the later is the circumference
    Trev=1./frev;
    
    Qx=opt.TUNES[0]
    Qy=opt.TUNES[1]
    Qx=Qx-floor(Qx)
    Qy=Qy-floor(Qy)
    
    
    # loop on filenames
    for filename in opt.FILE:
    
	# text for the legend
        fil=filename.replace("_"," ").replace(".data","").replace("/"," ").replace("-"," ");
	
    	# read BPM file
	nturns,nbunch,x,y,bpm=read_BPM(filename)
	
	# extract turns of interest
	#x=x[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(x))];
	#y=y[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(y))];
	
	# convert to 2D matrix (nturns*nbunch)
	x2d=np.transpose(np.reshape(x,(nbunch,-1)))
	y2d=np.transpose(np.reshape(y,(nbunch,-1)))
	
	
	if (opt.RAW):
	
	    # plot 1D data (along the bunch train, for the last turn)
	    fig1x,ax1x=init_figure();
	    fig1y,ax1y=init_figure();

	    plot_bunch(nbunch,x[nturns-1::nturns],'Turn '+str(nturns-1),'-b','Average x position of the bunch [a.u]','Raw position along the bunch train, '+fil+', horizontal',ax1x,bunchoffset=opt.BOFF);
	    plot_bunch(nbunch,y[nturns-1::nturns],'Turn '+str(nturns-1),'-b','Average y position of the bunch [a.u.]','Raw position along the bunch train, '+fil+', vertical',ax1y,bunchoffset=opt.BOFF);
	    
	    # plot 1D data (turn by turn, for the last bunch)
	    fig1x,ax1x=init_figure();
	    fig1y,ax1y=init_figure();

	    plot(np.arange(nturns),x[-nturns:],'Bunch '+str(nbunch+opt.BOFF-1),'-b','Average x position of the bunch [a.u.]',ax1x,0);
	    ax1x.set_title('Raw turn by turn position, '+fil+', horizontal')
	    ax1x.legend(loc=0);
	    plot(np.arange(nturns),y[-nturns:],'Bunch '+str(nbunch+opt.BOFF-1),'-b','Average y position of the bunch [a.u.]',ax1y,0);
	    ax1y.set_title('Raw turn by turn position, '+fil+', vertical')
	    ax1y.legend(loc=0);
    
	    # plot 2D raw data
	    fig2x,ax2x=init_figure();
	    fig2y,ax2y=init_figure();

	    plot2D(np.transpose(x2d),0,nturns-1,opt.BOFF,nbunch+opt.BOFF-1,'Turn number','Bunch number','Raw positions, '+fil+', horizontal',ax2x)		    
	    plot2D(np.transpose(y2d),0,nturns-1,opt.BOFF,nbunch+opt.BOFF-1,'Turn number','Bunch number','Raw positions, '+fil+', vertical',ax2y)		    


	if (opt.FFT1): 

	    x1=x2d.flatten()
	    y1=x2d.flatten()
	    
	    # 1D FFT for each bunch
	    fft1d_bunch(x1,nbunch,Qx,0.1,fil,'x');
	    fft1d_bunch(y1,nbunch,Qy,0.1,fil,'y');

	    # 1D FFT for each turn (along bunch train)
	    fft1d_turn(x1,nbunch,0,fil,'x')
	    fft1d_turn(y1,nbunch,0,fil,'y')
	
	
	if (opt.FFT2): 

	    # 2D FFT
	    fileout=open(filename[:-5]+'_'+opt.OUT,'w')

	    xmax=fft2d(x2d,Qx,0.1,fil,'x')
	    ymax=fft2d(y2d,Qy,0.1,fil,'y')
	    
	    print >> fileout, "2D FFT: for x & y: (mode,tune) of maximum: ", xmax, "\t", ymax;
	    print "2D FFT: for x & y: (mode,tune) of maximum: ", xmax, "\t", ymax;

	    fileout.close()
	
	
	
	if (opt.SVD):
	
	    # SVD (time domain analysis)
	    svd(x2d,Qx,0.1,0,nturns,fil,'x',bunchoffset=opt.BOFF)
	    svd(y2d,Qy,0.1,0,nturns,fil,'y',bunchoffset=opt.BOFF)
	    
    
			
    pylab.show();
