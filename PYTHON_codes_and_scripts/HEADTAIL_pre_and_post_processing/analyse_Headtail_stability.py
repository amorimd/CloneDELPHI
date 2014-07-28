#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

import pylab,re,random
import numpy as np
from string import split, replace
from parser_lib import *
import math
import matplotlib
from plot_lib import set_fontsize,init_figure,end_figure
from io_lib import list_files
from string_lib import takeout_common
from read_cfg import read_cfg
from read_Headtail_prt_fit import read_prt_file
from read_Headtail_prt import extract_Headtail_param,check_data,read_allprt_file
from construct_bunchtable_from_scheme import check_nbunch


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-n", "--nbeg",type=int,
                      help="Specify the number of turns at the beginning on which we compute the initial maximum of the actions (default=100)",
                      metavar="NBEG", default=100,dest="NBEG")
    parser.add_option("-o", "--output",help="Specify output suffix for stability file(s) (default=stable.txt)",
                      metavar="OUT",default="stable.txt",dest="OUT")
    parser.add_option("-r", "--ratio",type=float,
                      help="Specify the minimum ratio between the maximum of the whole data and that of of the first turns, in order to consider the beam as unstable (default=1.1)",
                      metavar="RATIO", default=1.1,dest="RATIO")		      
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    return opt, args



if __name__ == "__main__":
    opt,args=parsse(); 

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    # tables with stability flag (True=stable, False=unstable)
    tabstabx=[];tabstaby=[];


    # loop on filenames
    for ifile,filename in enumerate(listname):
    

	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    
  	# read prt file
	x,y,epsx,epsy,invarx,invary,z,bl,zemit,nprleft_frac,bin_frac=read_allprt_file(filename,nbunch)
	#x,xp,y,yp=read_prt_file(filename,nbunch)
	
	# maxima on first turns
	maxinix=np.max(np.abs(invarx[:nbunch*opt.NBEG]));
	maxiniy=np.max(np.abs(invary[:nbunch*opt.NBEG]));

	# maxima on the whole data
	maxtotx=np.max(np.abs(invarx));
	maxtoty=np.max(np.abs(invary));
	
	if (maxtotx>maxinix*opt.RATIO): stabx="UNSTABLE";
	else: stabx="STABLE";
	
	if (maxtoty>maxiniy*opt.RATIO): staby="UNSTABLE";
	else: staby="STABLE";
	
	filestab=open(filename[:-8]+'_'+opt.OUT,'w');
	print >> filestab, "x\t", stabx;
	print >> filestab, "y\t", staby;
	filestab.close();
	
	tabstabx.append(maxtotx<=maxinix*opt.RATIO);
	tabstaby.append(maxtoty<=maxiniy*opt.RATIO);
	
    

    if (len(listname)>1)and(listname[0].find('oct')!=-1):
    
	# take out all the common parameters in the names
	listless=takeout_common(listname);
	
	listparam=[]; # list of different parameters scanned except octupole current and Qsec
	ilist=0;
	staboctx=[];stabocty=[]; # stabilizing octupole currents
	unstaboctx=[];unstabocty=[]; # error on stabilizing currents
	
	# loop on filenames
	for ifile,filename in enumerate(listname):
	
	    tmp=split(listless[ifile],' ');
	    n=len(tmp);
	    
	    for k,item in enumerate(reversed(tmp)):
	    	i=n-k-1;
	    	if (len(item)==0): tmp.pop(i);
	    	if (item.startswith('oct')):
		    octstr=tmp.pop(i);
		    curoct=int(octstr.replace('oct',''));
		if (item.startswith('qsecx')): 
		    qsecxstr=tmp.pop(i);
		    qsecx=int(qsecxstr.replace('qsecx',''));
		if (item.startswith('qsecy')): 
		    qsecystr=tmp.pop(i);
		    qsecy=int(qsecystr.replace('qsecy',''));
	    
	    if not(tmp in listparam):
	    	listparam.append(tmp);
		ilist+=1;
		
		if (tabstabx[ifile]): staboctx.append(curoct);
		else:
		    if curoct>=0: staboctx.append(float('inf'));
		    else: staboctx.append(float('-inf'));
		
	    	if (tabstaby[ifile]): stabocty.append(curoct);
		else:
		    if curoct>=0: stabocty.append(float('inf'));
		    else: stabocty.append(float('-inf'));
		
		unstaboctx.append(0);
		unstabocty.append(0);
		
	    else:
	    	j=listparam.index(tmp);
		if (tabstabx[ifile]):
		    if curoct>=0: staboctx[j]=min(staboctx[j],curoct);
		    else: staboctx[j]=max(staboctx[j],curoct);
		else:
		    if curoct>=0: unstaboctx[j]=max(unstaboctx[j],curoct);
		    else: unstaboctx[j]=min(unstaboctx[j],curoct);
		
		if (tabstaby[ifile]):
		    if curoct>=0: stabocty[j]=min(stabocty[j],curoct);
		    else: stabocty[j]=max(stabocty[j],curoct);
		else:
		    if curoct>=0: unstabocty[j]=max(unstabocty[j],curoct);
		    else: unstabocty[j]=min(unstabocty[j],curoct);
    
	
	
	# columns headers
	col="";
	for item in listparam[0]: col=col+''.join([c for c in item if not c.isdigit()])+"\t";
	
	# write in a file
	filestab=open(opt.OUT,'w');
	print >> filestab, col, "stab_cur_x\tstab_cur_y\terr_x\terr_y";
	for i,param in enumerate(listparam):
	    par="";
	    for item in param: par=par+''.join([c for c in item if c.isdigit()])+"\t";
	    errx=(staboctx[i]-unstaboctx[i])/2.;
	    erry=(stabocty[i]-unstabocty[i])/2.;
	    print staboctx[i],unstaboctx[i],stabocty[i],unstabocty[i]
	    print >> filestab, par, "\t", staboctx[i]-errx, "\t", stabocty[i]-erry, "\t", np.sign(staboctx[i])*errx, "\t", np.sign(stabocty[i])*erry;
	filestab.close();
	
    	

    sys.exit()

