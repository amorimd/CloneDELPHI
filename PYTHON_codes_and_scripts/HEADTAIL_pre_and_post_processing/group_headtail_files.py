#!/usr/bin/python


# grouping together different growthrate and tuneshift files (spanning several parameters)
# in order to plot them easily later.

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
from optparse import OptionParser
from plot_lib import init_figure,end_figure,plot
from io_lib import list_files
from string_lib import takeout_common,split_and_takeout_spaces
from string import split

def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",action="append",
                      help="Specify the name of the .txt files to group. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--line",type=int,
                      help="Specify which line to read and group from the data files. Choose arbitrary large number (e.g. 50) for last line. Default=second line (just after the headers).",
                      metavar="LINE", default=1,dest="LINE")
    parser.add_option("-o", "--output",help="Specify suffix for output file names (default='.txt')",
                      default=".txt",dest="OUT")
    parser.add_option("-p", "--paramxaxis",
                      help="Specify which parameter is the one used for the column of the file. Default=take first changing parameter in the file names.",
                      metavar="PAR", default=None,dest="PAR")
    parser.add_option("-t", "--tau",action="store_true",
                      help="Specify if data is coming from raw data fit giving rise times (tau) instead of growth rates. In that case it converts them to growth rates.",
                      metavar="TAU", default=False,dest="TAU")
    (opt, args) = parser.parse_args()
    #print opt.FILE
    return opt, args


def read_tuneshifts_headers(line):

    # read column headers and identify them
    # in output: column numbers of tunes (x & y), growth rates (x & y) and sigmas (x & y)
    
    tuxcol=-1;tuycol=-1;rxcol=-1;rycol=-1;sigxcol=-1;sigycol=-1;
    ll=line.strip().split();
    for k,col in enumerate(ll):
	if (col.startswith('Tune'))or(col.startswith('tune')):
	    if (col.find('x')!=-1): tuxcol=k;
	    if (col.find('y')!=-1): tuycol=k;
	if col.startswith('Qx'): tuxcol=k;
	if col.startswith('Qy'): tuycol=k;
	if ((col.startswith('Growth'))or(col.startswith('growth')))or((col.startswith('Tau'))or(col.startswith('tau'))):
	    if (col.find('x')!=-1): rxcol=k;
	    if (col.find('y')!=-1): rycol=k;
	if (col.startswith('Sigma'))or(col.startswith('sigma')):
	    if (col.find('x')!=-1): sigxcol=k;
	    if (col.find('y')!=-1): sigycol=k;
	    
    return tuxcol,tuycol,rxcol,rycol,sigxcol,sigycol;
	

def extract_param(name):

    # extract name of parameter and its length in a string of format parameter_name+value
    # example: if name="I5p5" -> lenpar=1, parname="I"
    # particular case: when parname should be csi, it is replaced by "Qp"
    # example: if name="csi5p5" -> lenpar=3, parname="Qp"

    numbers=set('0123456789-');
    
    lenpar=0;
    while not(name[lenpar] in numbers): lenpar+=1;
    parname=name[:lenpar];
    if (parname=='csi')or(parname=='xi'): parname='Qp';
	    
    return lenpar,parname;


def read_tuneshifts_file_singleline(filename,numline=1):

    # read tuneshift and growth rate data of the file (actually take 
    # only one line, the one of number numline - starting from zero)
    # The first line should be the headers.
    
    file1=open(filename);
    N=len([l for l in file1]);
    if (numline>N-1): numline=N-1; # take last line if numline is too large
    file1.close();
    
    if (numline>0):
	file1=open(filename);
	for il,line in enumerate(file1):
    	    if (il==0): tuxcol,tuycol,rxcol,rycol,sigxcol,sigycol=read_tuneshifts_headers(line);
	    if (il==numline)and(len(split(line))>1):
		rx=float(split(line)[rxcol]);
		ry=float(split(line)[rycol]);
		tuneshiftx=float(split(line)[tuxcol]);
		tuneshifty=float(split(line)[tuycol]);
		if (sigxcol==-1): sigmax=0.;
		else: sigmax=float(split(line)[sigxcol]);
		if (sigycol==-1): sigmay=0.;
		else: sigmay=float(split(line)[sigycol]);

	file1.close();
    
        return tuneshiftx,tuneshifty,rx,ry,sigmax,sigmay;
	
    else: return 0.,0.,0.,0.,0.,0.;
		

if __name__ == "__main__":

    opt,args=parsse();
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    #print listname

    # create list with the changed parameter(s) (taking out all the common parameters
    # in the names)
    listpar,com=takeout_common(listname,flagcommon=True);
    com=[nm+'_' for nm in com];
    common=''.join(com);common=common[:-1];
    #print listpar,common
    
    # identify the parameter to be put in a column of the final file
    # (we identify on the first file name)
    firstfilepar=split_and_takeout_spaces(listpar[0]);
    #print firstfilepar,opt.PAR
    
    if (opt.PAR==None): 
        if (len(firstfilepar)>0):
	    # take first changing parameter as the parameter to be put in a column
	    ipar=0;
	    lenpar,parname=extract_param(firstfilepar[0]);
	else:
	    print "Pb in file names";sys.exit();
    
    else:
    	param=opt.PAR;
	for iname,name in enumerate(firstfilepar):
    	    if name.startswith(param):
		# distinguish name of the parameter from numeric value associated with it
		ipar=iname;
		lenpar,parname=extract_param(name);
    
    
    # initialization of resulting lists
    outfiles=[];filenames=[];paramcol={};
    tuneshiftx={};tuneshifty={};rx={};ry={};sigmax={};sigmay={};
    # main loop on files to collect all data and open output files
    for iname,name in enumerate(listpar):

    	nameless=split_and_takeout_spaces(name);
	par=nameless.pop(ipar)
	paramvalue=float(par[lenpar:].replace('p','.'));
	nameless=['_'+nm for nm in nameless];
	newname=common+''.join(nameless);
	
	if not(newname in filenames):
    	    # create output file
	    filenames.append(newname);
    	    outfiles.append(open(newname+opt.OUT,'w'));
	    print >> outfiles[-1], parname, "\tTune_shiftx\tGrowthrate_x[s-1]\tSigma_tune_shiftx\tTune_shifty\tGrowthrate_y[s-1]\tSigma_tune_shifty"

	    # initialize all columns
	    paramcol[newname]=[];
	    tuneshiftx[newname]=[];
	    tuneshifty[newname]=[];
	    rx[newname]=[];
	    ry[newname]=[];
	    sigmax[newname]=[];
	    sigmay[newname]=[];

	# read the data file (actually take only line of number opt.LINE, starting from zero) and add to the columns
	tuneshiftx1,tuneshifty1,rx1,ry1,sigmax1,sigmay1=read_tuneshifts_file_singleline(listname[iname],numline=opt.LINE);
	
	if (abs(tuneshiftx1)+abs(tuneshifty1)+abs(rx1)+abs(ry1)+abs(sigmax1)+abs(sigmay1) !=0.):
	    paramcol[newname].append(paramvalue);
	    tuneshiftx[newname].append(tuneshiftx1);
	    tuneshifty[newname].append(tuneshifty1);
	    rx[newname].append(rx1);
	    ry[newname].append(ry1);
	    sigmax[newname].append(sigmax1);
	    sigmay[newname].append(sigmay1);
	

    # finally write data in files
    for ifile,file1 in enumerate(outfiles):
    
    	filename=filenames[ifile];
	print "Writing file ",filename+opt.OUT

	# conversion to array
	paramcol[filename]=np.array(paramcol[filename]);
	if opt.TAU:
	    rx[filename]=1./np.array(rx[filename]);
	    ry[filename]=1./np.array(ry[filename]);
	else:
	    rx[filename]=np.array(rx[filename]);
	    ry[filename]=np.array(ry[filename]);
	sigmax[filename]=np.array(sigmax[filename]);
	sigmay[filename]=np.array(sigmay[filename]);
	tuneshiftx[filename]=np.array(tuneshiftx[filename]);
	tuneshifty[filename]=np.array(tuneshifty[filename]);
	
    	# sort according to ascending parameter (in paramcol)
	ind=np.argsort(paramcol[filename]);
	paramcol[filename]=paramcol[filename][ind];
	rx[filename]=rx[filename][ind];
	ry[filename]=ry[filename][ind];
	tuneshiftx[filename]=tuneshiftx[filename][ind];
	tuneshifty[filename]=tuneshifty[filename][ind];
	sigmax[filename]=sigmax[filename][ind];
	sigmay[filename]=sigmay[filename][ind];
	
	# write to final file
	for i,par in enumerate(paramcol[filename]):
	    print >> file1, par, "\t", tuneshiftx[filename][i], "\t", rx[filename][i], "\t", sigmax[filename][i], "\t", tuneshifty[filename][i], "\t", ry[filename][i], "\t", sigmay[filename][i];


	file1.close();    


    sys.exit();
