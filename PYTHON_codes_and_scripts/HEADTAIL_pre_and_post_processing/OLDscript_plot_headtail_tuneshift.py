#!/usr/bin/python


# obsolete

import sys,
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
import numpy as np
from parser_lib import *
from plot_lib import init_figure,end_figure,plot
from string_lib import takeout_common,split_and_takeout_spaces
from io_lib import list_files
from string import split

def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the name of the txt files to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
#    parser.add_option("-p", "--paramxaxis",type=int,
#                      help="Specify which changed parameter (in the order of apearance in the file names) is the one used for the x-axis. Default = 0 (first one).",
#                      metavar="PAR", default=0,dest="PAR")
#    parser.add_option("-m", "--maincolumn",type=int,
#                      help="Specify which column to use in the txt file for the y-axis data. Default = 0 (first one).",
#                      metavar="MAINCOL", default=0,dest="MAINCOL")
    parser.add_option("-l", "--legend",action="callback",callback=multistring_parse,
                      help="Specify legends when second parameter (default=parameter names from the filenames).",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-p", "--paramcol",type=int,nargs=2,
                      help="Specify the parameter used for the x-axis (first number), and the one used for the legend (2nd number). Default=[0,1].",
                      metavar="PARCOL", default=[0,1],dest="PARCOL")
    parser.add_option("-r", "--ratecolumns",type=int,nargs=2,
                      help="Specify which columns (beginning at 0) to use in the txt file for the growth rate data (2 cols: x and y). Default = [1,2].",
                      metavar="GRCOL", default=[1,2],dest="GRCOL")
    parser.add_option("-s", "--sigmacolumns",type=int,nargs=2,
                      help="Specify which columns (beginning at 0) to use in the txt file for the errror on the tuneshifts data (2 cols: x and y). Default = None (no errorbar).",
                      metavar="SIGCOL", default=None,dest="SIGCOL")
    parser.add_option("-t", "--tunecolumns",type=int,nargs=2,
                      help="Specify which columns (beginning at 0) to use in the txt file for the tuneshifts data (2 cols: x and y). Default = [3,4].",
                      metavar="TUCOL", default=[3,4],dest="TUCOL")
    (opt, args) = parser.parse_args()
    #print opt.FILE
    return opt, args


if __name__ == "__main__":

    opt,args=parsse();
    
    rxcol=opt.GRCOL[0];rycol=opt.GRCOL[1];
    tuxcol=opt.TUCOL[0];tuycol=opt.TUCOL[1];
    if (opt.SIGCOL!=None):
    	sigxcol=opt.GRCOL[0];sigycol=opt.GRCOL[1];
    
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    print listname

    figx,axx=init_figure();
    figy,axy=init_figure();
    figrx,axrx=init_figure();
    figry,axry=init_figure();

    # create list with the changed parameter(s) (taking out all the common parameters
    # in the names)
    listpar=takeout_common(listname);
    print listpar
    
    tmp=split_and_takeout_spaces(listpar[0]);
    print tmp,opt.PARCOL
    param=[];
    if (tmp[opt.PARCOL[0]].startswith('I')):
	parname='Intensity';
	parlabel="Intensity $ (10^{11} $ p+/b)";
	ipar=1;
    if (tmp[opt.PARCOL[0]].startswith('csi')):
	parname='Qprime';
	parlabel=" $ Q' $ ";
	ipar=3;

    if (len(tmp)==1):
	for il,line in enumerate(listpar):
	    tmp2=split_and_takeout_spaces(line);tmp2=tmp2[0];
	    param.append(float(tmp2[ipar:].replace('p','.')));

	tuneshiftx=[];tuneshifty=[];rx=[];ry=[];sigmax=[];sigmay=[];

	for iname,filename in enumerate(listname):

            # read data of this file (actually take only second line)
	    for il,line in enumerate(open(filename)):
    		if (il==1)and(len(split(line))>1):
		    rx.append(float(split(line)[rxcol]))
		    ry.append(float(split(line)[rycol]))
		    tuneshiftx.append(float(split(line)[tuxcol]))
		    tuneshifty.append(float(split(line)[tuycol]))
		    if (opt.SIGCOL!=None):
			sigmax.append(float(split(line)[sigxcol]))
			sigmay.append(float(split(line)[sigycol]))
		    else:
			sigmax.append(0.)
			sigmay.append(0.)


	axx.errorbar(param,tuneshiftx,yerr=sigmax,fmt='xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);
	axrx.plot(param,rx,'xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);

	axy.errorbar(param,tuneshifty,yerr=sigmay,fmt='xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);
	axry.plot(param,ry,'xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);

	if (parname=='Intensity'):
	    px=np.polyfit(param,tuneshiftx,1);
	    py=np.polyfit(param,tuneshifty,1);
	    axx.plot(param,px[1]+px[0]*np.array(param),'-b',label='Fit by a-x*%.6f'%(-px[0])+" $ /10^{11} $ ",lw=3.);
	    axy.plot(param,py[1]+py[0]*np.array(param),'-b',label='Fit by a-x*%.6f'%(-py[0])+" $ /10^{11} $ ",lw=3.);

    elif (len(tmp)>1):
	# identify a second parameter if it exists (-> different curves with legends)
	param2=[];param={};rx={};ry={};
	tuneshiftx={};tuneshifty={};sigmax={};sigmay={};
	if (opt.LEG==None): leg=[];
	else: leg=opt.LEG;
	
	for ip,par in enumerate(listpar):
	    tmp2=split_and_takeout_spaces(par);
	    if not(tmp2[opt.PARCOL[1]] in param2):
	    	param2.append(tmp2[opt.PARCOL[1]]);param[tmp2[opt.PARCOL[1]]]=[];
		if (opt.LEG==None): leg.append(tmp2[opt.PARCOL[1]]);
		rx[tmp2[opt.PARCOL[1]]]=[];ry[tmp2[opt.PARCOL[1]]]=[];
		tuneshiftx[tmp2[opt.PARCOL[1]]]=[];tuneshifty[tmp2[opt.PARCOL[1]]]=[];
		sigmax[tmp2[opt.PARCOL[1]]]=[];sigmay[tmp2[opt.PARCOL[1]]]=[];
		
	    param[tmp2[opt.PARCOL[1]]].append(float(tmp2[opt.PARCOL[0]][ipar:].replace('p','.')))

            # read data of this file (actually take only second line)
	    for il,line in enumerate(open(listname[ip])):
    		if (il==1)and(len(split(line))>1):
		    rx[tmp2[opt.PARCOL[1]]].append(float(split(line)[rxcol]))
		    ry[tmp2[opt.PARCOL[1]]].append(float(split(line)[rycol]))
		    tuneshiftx[tmp2[opt.PARCOL[1]]].append(float(split(line)[tuxcol]))
		    tuneshifty[tmp2[opt.PARCOL[1]]].append(float(split(line)[tuycol]))
		    if (opt.SIGCOL!=None):
			sigmax[tmp2[opt.PARCOL[1]]].append(float(split(line)[sigxcol]))
			sigmay[tmp2[opt.PARCOL[1]]].append(float(split(line)[sigycol]))
		    else:
			sigmax[tmp2[opt.PARCOL[1]]].append(0.)
			sigmay[tmp2[opt.PARCOL[1]]].append(0.)

	pat=['-xb','-or','-+g','-dk','-vm','-^c','-.y'];

        for ipar2,par2 in enumerate(param2):
	    
	    # conversion to array
	    param[par2]=np.array(param[par2]);
	    rx[par2]=np.array(rx[par2]);
	    ry[par2]=np.array(ry[par2]);
	    sigmax[par2]=np.array(sigmax[par2]);
	    sigmay[par2]=np.array(sigmay[par2]);
	    tuneshiftx[par2]=np.array(tuneshiftx[par2]);
	    tuneshifty[par2]=np.array(tuneshifty[par2]);
    	    # sort
	    ind=np.argsort(param[par2]);
	    param[par2]=param[par2][ind];
	    rx[par2]=rx[par2][ind];
	    ry[par2]=ry[par2][ind];
	    tuneshiftx[par2]=tuneshiftx[par2][ind];
	    tuneshifty[par2]=tuneshifty[par2][ind];
	    sigmax[par2]=sigmax[par2][ind];
	    sigmay[par2]=sigmay[par2][ind];
	
	    axx.errorbar(param[par2],tuneshiftx[par2],yerr=sigmax[par2],fmt=pat[ipar2],label=leg[ipar2],lw=3.,ms=10.,mew=2.5);
	    axrx.plot(param[par2],rx[par2],pat[ipar2],label=leg[ipar2],lw=3.,ms=10.,mew=2.5);

	    axy.errorbar(param[par2],tuneshifty[par2],yerr=sigmay[par2],fmt=pat[ipar2],label=leg[ipar2],lw=3.,ms=10.,mew=2.5);
	    axry.plot(param[par2],ry[par2],pat[ipar2],label=leg[ipar2],lw=3.,ms=10.,mew=2.5);


    axx.set_xlabel(parlabel);
    axx.set_ylabel("Horizontal tune shift");
    axrx.set_xlabel(parlabel);
    axrx.set_ylabel("Horizontal growth rate");

    axy.set_xlabel(parlabel);
    axy.set_ylabel("Vertical tune shift");
    axry.set_xlabel(parlabel);
    axry.set_ylabel("Vertical growth rate");

    end_figure(figx,axx);
    end_figure(figy,axy);
    end_figure(figrx,axrx);
    end_figure(figry,axry);
    pylab.show();
    
