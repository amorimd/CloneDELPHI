#!/usr/bin/python


# to plot tuneshifts and growth rates from several files, vs. a parameter changing
# in one of the columns (user defined)

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
from parser_lib import *
from plot_lib import init_figure,end_figure,plot
from io_lib import list_files
from string import split
from string_lib import split_and_takeout_spaces,takeout_common,takeout_spaces
from group_headtail_files import read_tuneshifts_headers

def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the name of the .txt files to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--legend",action="callback",callback=multistring_parse,
                      help="Specify legends (can be several, in one -l option). Default=parameter names from the filenames.",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-p", "--paramxaxis",
                      help="Specify the column name for the x-axis parameter. Default= first column. Note: if name given (e.g. Qp) applies for two columns, the column with its header containing 'x' is used for the x figures, idem with 'y'.",
                      metavar="PAR", default=None,dest="PAR")
    parser.add_option("-x", "--xlabel",
                      help="Specify the x-axis label. Default= determined from the -p option (parameter name)",
                      metavar="XLAB", default=None,dest="XLAB")
    (opt, args) = parser.parse_args()
    #print opt.FILE
    return opt, args


def read_tuneshifts_file(filename,param=None):

    # read tuneshift and growth rate data of the file (all lines)
    # param is the name of a varying parameter to indetify the corresponding header(s):
    # if there hare two headers beginning with param, one is for 'x' and the other one for 'y'

    file1=open(filename);
    parx=[];pary=[];tuneshiftx=[];tuneshifty=[];
    rx=[];ry=[];sigmax=[];sigmay=[];

    for il,line in enumerate(file1):

    	if (il==0):
	    # identify columns with tuneshifts, growth rates and sigmas
	    tuxcol,tuycol,rxcol,rycol,sigxcol,sigycol=read_tuneshifts_headers(line);
	    
	    # identify column with x-axis parameter
	    if (param!=None):
	    
		ll=line.strip().split();par=[];parcol=[];
		for k,col in enumerate(ll):
		    if (col.startswith(param)):
			par.append(col);parcol.append(k);

		if (len(par)==1):
	    	    par.append(par[0]);parcol.append(parcol[0]);
		elif (par[1].find('x')!=-1)or(par[0].find('y')!=-1):
		    par2=par[0];par[0]=par[1];par[1]=par2;
		    parcol2=parcol[0];parcol[0]=parcol[1];parcol[1]=parcol2;
		    
	    else:
	    	parcol[0]=0;parcol[1]=0;
		

	elif (len(split(line))>1):
	    parx.append(float(split(line)[parcol[0]]));
	    pary.append(float(split(line)[parcol[1]]));
	    rx.append(float(split(line)[rxcol]));
	    ry.append(float(split(line)[rycol]));
	    tuneshiftx.append(float(split(line)[tuxcol]));
	    tuneshifty.append(float(split(line)[tuycol]));
	    if (sigxcol==-1): sigmax.append(0.);
	    else: sigmax.append(float(split(line)[sigxcol]));
	    if (sigycol==-1): sigmay.append(0.);
	    else: sigmay.append(float(split(line)[sigycol]));
	
    file1.close();
    
    return par,(np.array(parx),np.array(pary)),(np.array(tuneshiftx),np.array(tuneshifty)),(np.array(rx),np.array(ry)),(np.array(sigmax),np.array(sigmay));	
		

if __name__ == "__main__":

    opt,args=parsse();
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE,flagprint=True);

    # initialize figures
    figtu=[];figr=[];axtu=[];axr=[];planes=[];
    for i in range(2):
    	fig,ax=init_figure();
	figtu.append(fig);axtu.append(ax);
    	fig,ax=init_figure();
	figr.append(fig);axr.append(ax);
	if (i==0): planes.append('x');
	else: planes.append('y');

    # create list with the changed parameter(s) (taking out all the common parameters
    # in the names)
    listpar=takeout_common(listname,nslash=1);

    # create legends
    if (opt.LEG==None)or(len(opt.LEG)<len(listname)): listleg=takeout_spaces(listpar);
    else: listleg=opt.LEG;
    
    print listleg;
    
    # patterns for the various curves
    pat=['-xb','-or','-+g','-dk','-vm','-^c','-.y'];
    patfit=['--b','--r','--g','--k','--m','--c','--y'];
    
    # main loop
    for iname,filename in enumerate(listname):

        # read data of the file
	parname,par,tuneshift,r,sigma=read_tuneshifts_file(filename,param=opt.PAR);
	
	for iplane,plane in enumerate(planes):

	    axtu[iplane].errorbar(par[iplane],tuneshift[iplane],yerr=sigma[iplane],fmt=pat[iname%7],label=listleg[iname],lw=3.,ms=10.,mew=2.5);
	    axr[iplane].plot(par[iplane],r[iplane],pat[iname%7],label=listleg[iname],lw=3.,ms=10.,mew=2.5);

	    if (parname[iplane].startswith('I')):
	        # case of intensity parameter: do a fit
		p=np.polyfit(par[iplane],tuneshift[iplane],1);
		axtu[iplane].plot(par[iplane],p[1]+p[0]*np.array(par[iplane]),patfit[iname%7],label=listleg[iname]+', fit by a-x*%.6f'%(-p[0])+" $ /10^{11} $ ",lw=3.);


    # finalize plots
    for iplane,plane in enumerate(planes):

	if (opt.XLAB==None):
	
	    if (parname[iplane].startswith('I')):
		parlabel="Intensity $ (10^{11} $ p+/b)";

	    elif (parname[iplane].startswith('csi'))or(parname[iplane].startswith('Qp')):
		parlabel=" $ Q^'_"+plane+" $ ";

	    else: parlabel=parname[iplane];
	    
	else: parlabel=opt.XLAB;

	axtu[iplane].set_xlabel(parlabel);
	axtu[iplane].set_ylabel(" $\Delta Q_"+plane+" $ ");
	axr[iplane].set_xlabel(parlabel);
	axr[iplane].set_ylabel("Growth rate in "+plane+" [s $ ^{-1} $]");

	end_figure(figtu[iplane],axtu[iplane]);
	end_figure(figr[iplane],axr[iplane]);

    
    pylab.show();
    
