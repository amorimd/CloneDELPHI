#!/usr/bin/python2.6


# grouping together different files (spanning several parameters)
# in order to plot them easily later.
# assume the columns are all the same in the files

import sys,pylab,re,dateutil,random,pytz,os
from datetime import time,datetime,date
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
import numpy as np
from optparse import OptionParser
from plot_lib import init_figure,end_figure,plot
from io_lib import list_files,read_file_singleline
from string_lib import takeout_common,split_and_takeout_spaces

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
    (opt, args) = parser.parse_args()
    #print opt.FILE
    return opt, args



def extract_param_simple(name):

    # extract name of parameter and its length in a string of format parameter_name+value
    # example: if name="I5p5" -> lenpar=1, parname="I"

    numbers=set('0123456789-');
    
    lenpar=0;
    while not(name[lenpar] in numbers): lenpar+=1;
    parname=name[:lenpar];
	    
    return lenpar,parname;



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
	    lenpar,parname=extract_param_simple(firstfilepar[0]);
	else:
	    print "Pb in file names";sys.exit();
    
    else:
    	param=opt.PAR;
	for iname,name in enumerate(firstfilepar):
    	    if name.startswith(param):
		# distinguish name of the parameter from numeric value associated with it
		ipar=iname;
		lenpar,parname=extract_param_simple(name);
        
    lines=[];paramcol=[];
    # main loop on files to collect all data and open output files
    for iname,name in enumerate(listpar):

    	nameless=split_and_takeout_spaces(name);
	par=nameless.pop(ipar)
	paramvalue=float(par[lenpar:].replace('p','.'));
	nameless=['_'+nm for nm in nameless];
	newname=common+''.join(nameless);
	
	if (iname==0):
    	    # create output file
	    filename=newname;
    	    outfile=open(newname+opt.OUT+'.txt','w');
	    headers=read_file_singleline(listname[iname],numline=0);
	    print >> outfile, parname, "\t", headers;

	# read the data file (actually take only line of number opt.LINE, starting from zero) and add to the columns
	lines.append(read_file_singleline(listname[iname],numline=opt.LINE));
	
	# parameter for first column
	paramcol.append(paramvalue);
	

    # finally write data in the output file (sorting the parameter first)
    print "Writing file ",filename+opt.OUT

    # conversion to arrays
    paramcol=np.array(paramcol);
    lines=np.array(lines);
	
    # sort according to ascending parameter (in paramcol)
    ind=np.argsort(paramcol);
    paramcol=paramcol[ind];
    lines=lines[ind];
	
    # write to final file
    for i,par in enumerate(paramcol):
	print >> outfile, par, "\t", lines[i];


    outfile.close();    


    sys.exit();
