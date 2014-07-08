#!/usr/bin/python

# library with in/out routines (read, write files)

import sys
import commands
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
from string import split,replace
import re,glob,os
import numpy as np

def list_files(files,flagprint=False):

    ''' list the files to be read. "files" is a list of names that can contain regular expressions
    (like *, etc.) from which the list of files to be analysed is constructed '''
    
    listname=[];
    for f in files:
    	l=glob.glob(f);
	l.sort();
	listname.extend(l);
    
    if flagprint:
    	for j in listname: print j;
    
    return listname;
    

def read_ncol_file(filename,ignored_rows=0,ignored_cols=0,delimiter=None):

    ''' kind of equivalent of MATLAB dlmread
    read a file in columns separated by spaces or tab (if 'delimiter' is None) or
    by any kind of delimiter specified in 'delimiter', ignoring a certain
    number of initial rows and columns
    it stops as soon as it cannot split a line into floats, or if a line is empty '''
    
    s=[];
    for il,line in enumerate(open(filename)):
    	if (il>=ignored_rows):
	    
	    try:
	    	s.append(np.array(split(line,delimiter)[ignored_cols:],dtype=float));
	    except ValueError:
		break;
	    
	    if (len(s[-1])==0): s.pop();break;

    return np.array(s);
    

def write_ncol_file(filename,data,header=None,delimiter="\t",format="%15.10e"):

    ''' kind of equivalent of MATLAB dlmwrite
    write from the 2D array "data" a file in columns separated by 'delimiter',
    with the first line given by "header"
    format is the data format to use when writing (e.g. %lf, %d, etc.) '''
    
    fileout=open(filename,'w');
    if (header!=None): print >> fileout, header;
    
    for irow in range(len(data[:,0])):
    	for icol in range(len(data[irow,:])-1): fileout.write(format % (data[irow,icol]));fileout.write(delimiter);
    	# last column separately
	fileout.write(format % (data[irow,-1])+"\n");

    fileout.close();
        

def write_Timber(t,data,filename,varname):

    ''' Write in "filename" the data in "data", vs time "t" (in s), in a Timber-like way
    The variable name is "varname".'''
    from io_lib import tb_

    fileTimber=open(filename,'a');
    print >> fileTimber, "VARIABLE: "+varname;
    print >> fileTimber, "";
    print >> fileTimber, "Timestamp (LOCAL_TIME),Value";
    for it,t in enumerate(t):
    	fracsec=round((t-np.floor(t))*1000)/1000.;
	un=int(fracsec>=1);
	fracsecstr=str(fracsec-un)
	for i in range(5-len(fracsecstr)): fracsecstr=fracsecstr+'0';
    	tprint=tb_(np.floor(t)+un)+fracsecstr.replace('0.','.')
	#print fracsec,tprint,t,np.floor(t),t-np.floor(t),round((t-np.floor(t))*1000)/1000.;
    	print >> fileTimber,tprint+","+str(data[it]);

    fileTimber.close();
	
	
def read_file_singleline(filename,numline=1):

    ''' read data of the file (actually take 
    only one line, the one of number numline - starting from zero)'''
    
    file1=open(filename);
    N=len([l for l in file1]);
    if (numline>N-1): numline=N-1; # take last line if numline is too large
    file1.close();
    
    file1=open(filename);
    for il,line in enumerate(file1):
	if (il==numline)and(len(split(line))>0): lineout=line;

    file1.close();
    if lineout[-1]=='\n': lineout=lineout[:-1];
    
    return lineout;


def find_string_in_file(filename,strin):
    ''' look for first line containing 'strin' in the file of name 'filename'
    return its line number (first line has number zero)'''
    lineres=-1;

    f=open(filename,'r')
    for il,line in enumerate(f.readlines()):
    	if (line.find(strin)!=-1)and(lineres==-1): lineres=il;

    f.close();
    return lineres;


def read_ncol_file_identify_header(filename,colheader,dispflag=True):
    ''' read file organized in columns with a header line, and extract the column
    which has [colheader] as header (or starts with it).
    if dispflag is True, print warning when there is no such column'''

    result=[];kcol=-1;
    flagstr=False; # true if final result contains alphanumeric character(s) (then it is returned as list of string)
    
    fid=open(filename,"r");

    for j,l in enumerate(fid.readlines()):
	ll=l.strip().split();
	#sys.stdout.write(l);sys.stdout.flush();
    	if (j==0):
	    # find proper column name
	    for k,col in enumerate(ll):
	    	if (re.compile(colheader).match(col)!=None): kcol=k;
	    # if first column is '#', do not take it into account in numbering
	    if ll[0].startswith('#'):
	    	kcol-=1;
    	else:
	    if (kcol>-1):
	    	try:
    		    res=float(ll[kcol]); # interpret as float
		except ValueError:
		    res=ll[kcol];flagstr=True;

		result.append(res);
	        #if (any(c.isalpha() for c in ll[kcol])): result.append(ll[kcol]);flagstr=True;
		#else: result.append(float(ll[kcol])); # if no letter present -> interpret as float

    fid.close()
    
    if (kcol<=-1)and(dispflag): print "In file "+filename+", column "+colheader+" not found; result will be empty";print result;
    
    if flagstr: return result;
    else: return np.array(result);
    

def test_and_create_dir(name):

    ''' test if a directory of name 'name' exists in the current directory,
    and if not, create it'''
    if not(os.path.exists(name)): os.makedirs(name);
    
    return;
   
    
