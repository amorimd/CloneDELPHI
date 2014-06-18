#!/usr/bin/python2.6

import sys
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab
import numpy as np
import math,cmath
from string import split, replace
from optparse import OptionParser
from plot_lib import set_fontsize,init_figure,end_figure,build_colors
from io_lib import list_files
from tables_lib import create_list


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the file name (collimator settings). Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--legend",action="append",
                      help="Specify the legend for the plot, for each settings file (one -l option per file). Default = name of the file.",
                      metavar="LEG", default=None,dest="LEG")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def plot(halfgap,leg,col,widt,line):
    # plot coll. settings with labels, legend and color
    pylab.bar(np.arange(len(halfgap)),halfgap,facecolor='none',width=widt,label=leg,edgecolor=col,ls=line,lw=2.5)


def read(filename):
    # read collimator settings file
    halfgap=[];name=[]
    fh=open(filename,"r")
    for j,l in enumerate(fh.readlines()):
	ll=l.strip().split();
	#sys.stdout.write(l);sys.stdout.flush();
    	if (j==0):
	    # find column names "name" and "halfgap"
	    for k,col in enumerate(ll):
	    	if col.startswith('name'): namecol=k;
		if col.startswith('halfgap'): hgcol=k;
	    # if first column is '#', do not take it into account in numbering
	    if ll[0].startswith('#'): namecol=namecol-1;hgcol=hgcol-1;
    	else:
	    a=ll[namecol].replace(".B1","").replace(".B2","");
	    if (a.endswith(".B")): a=a[:-2];
	    if not( ( (a.startswith("TCSG.A4R7"))or(a.startswith("TCSG.A4L7")) )
	    	or( (a.startswith("TCSG.B5R7"))or(a.startswith("TCSG.B5L7"))) ):
	        a=a[:-2]+a[-1:];
	    if ((a.startswith("TDI"))or(a.startswith("TCLI"))):
	    	a=a[:-1];
	    name.append(a);
	    halfgap.append(ll[hgcol]);

    fh.close()
    return name,halfgap
    
                  
if __name__ == "__main__":
    opt,args=parsse();
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);

    if (len(listname)<5):
    	color=['b','r','g','m'];
    	line=['solid','dashed','dotted','dashdot'];
    else:
    	color=build_colors(len(listname),randomize=True);
	line=create_list('solid',n=len(listname));

    #fig=pylab.figure();
    #pylab.axes([0.03,0.05,0.96,0.9])
    fig,ax=init_figure(axes=[0.1,0.15,0.85,0.8])
    width=0.8; # width of the bars of the histogram

    for j,fil in enumerate(listname):

	# string for legend
	k=fil.rfind("/");
	if opt.LEG==None: leg=fil[k+1:].replace("_"," ").replace(".txt","");
	else: leg=opt.LEG[j];
	
	# read collimator settings
	name,halfgap=read(fil);
	
	if (j==0): name0=name;
	else:
	    halfgapnew=[];
	    for nam in name0:
	    	if (name.count(nam)>=1): halfgapnew.append(halfgap[name.index(nam)]);
		else: halfgapnew.append(0.0);
	    halfgap=halfgapnew;
	
	# plot
	plot(np.array([float(i)*1.e3 for i in halfgap]),leg,color[j],width,line[j]);

    
    pylab.xlabel("collimator name");
    pylab.ylabel("half-gap (mm)");
    pylab.legend(loc=0);
    set_fontsize(fig,'xx-large');
    pylab.xticks(np.arange(len(name0))+width/2., name0, rotation=75, size=10)
    ax.xaxis.labelpad=12;
    ax.yaxis.labelpad=12;
    #pylab.ylim([0,30]);
	
    pylab.show();sys.exit()

