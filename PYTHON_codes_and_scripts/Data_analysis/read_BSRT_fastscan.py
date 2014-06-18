#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,os,re,dateutil,random,pytz
from datetime import time,datetime,date
import numpy as np
import numpy.ma as ma
from string import split, replace
from optparse import OptionParser
import math
import matplotlib
import matplotlib.dates
#import subprocess
from plot_lib import plot,init_figure,end_figure
from collimator_settings import tt
from io_lib import list_files,read_ncol_file,write_ncol_file,write_Timber,find_string_in_file
from tables_lib import complementary
from plot_correlations_multibunch import intersect

def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--average",action="store_true",
                      help="Specify if we average all the files chosen. Then sigma is also computed (default = no averaging). Everything is then put two 3 columns files (slot number, average emit, sigma, for either H or V).",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-b", "--beam",
                      help="Specify the beam (1 or 2)",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-o", "--output",
    		      help="Specify the end of the output filename suffix for TIMBER-like file of emittances vs time (default = nothing)",
                      default="",dest="OUT")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the BSRT fast scan name of the file. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-r", "--reffile",
                      help="Specify the reference file with the filling scheme in the first column (and a one line header). Typically from read_FBCT.py (average intensities per bunch). Used only with -a option.",
                      metavar="REF", default=None,dest="REF")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


if __name__ == "__main__":
    opt,args=parsse();

    if opt.BEAM=="1":
        beam="B1";
    elif opt.BEAM=="2":
        beam="B2";
    else: print "specify beam 1 or 2"; sys.exit()

    gmt=pytz.timezone('Europe/Amsterdam');

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if opt.AVER:
	# read reference file with filling scheme
	sref=read_ncol_file(opt.REF,ignored_rows=1);
	# 25ns slot numbers
	slots=sref[:,0];
	# initialize bunch-by-bunch emittances tables (masked arrays)
	emith=ma.zeros((len(slots),len(listname)),dtype=float);
	emitv=ma.zeros((len(slots),len(listname)),dtype=float);
    
    
    # initialization for averages vs bunches
    time=[];avEMITh=[];avEMITv=[];sigEMITh=[];sigEMITv=[];
    
    for ifile,filename in enumerate(listname):
    
    	print filename;
	
	# find first line of horizontal emittance
	l_emith=find_string_in_file(filename,"Horizontal Emit");
	l_emith+=1;
	# extract horizontal emittance part
	s=read_ncol_file(filename,ignored_rows=l_emith,delimiter=',');
	# take out too large or too small values
	ind=pylab.mlab.find((s[:,1]-6)*(s[:,1]-1)<=0.);sh=s[ind,:];
	# write it on a separate file
	write_ncol_file(filename.replace('.dat','_EMIT_H.dat'),sh,header="25ns_slot\thor_emit[mm.mrad]");
	avEMITh.append(np.average(sh[:,1]));sigEMITh.append(np.sqrt(np.var(sh[:,1])));
	
    	# find first line of horizontal emittance
	l_emitv=find_string_in_file(filename,"Vertical Emit");
	l_emitv+=1;
	# extract horizontal emittance part
	s=read_ncol_file(filename,ignored_rows=l_emitv,delimiter=',');
	# take out too large or too small values
	ind=pylab.mlab.find((s[:,1]-5)*(s[:,1]-1)<=0.);sv=s[ind,:];
	# write it on a separate file
	write_ncol_file(filename.replace('.dat','_EMIT_V.dat'),sv,header="25ns_slot\tver_emit[mm.mrad]");
	
	# extract date and time of fast scan (from filename)
	time.append(pylab.date2num(tt(filename[-28:-9].replace('_',' ').replace('.','-')))*86400.);
	# compute average emittances of all bunches, and sigma
	avEMITv.append(np.average(sv[:,1]));sigEMITv.append(np.sqrt(np.var(sv[:,1])));
	
	if opt.AVER:
	    # fill table of hor. emittance, for later computation of time average and sigma
	    s,ind1,ind2=intersect(slots,sh[:,0]);
	    emith[ind1,ifile]=sh[ind2,1];
	    # mask whenever the values are missing
	    ind3=complementary(ind1,np.arange(len(slots)));
	    if len(ind3)>0: emith[ind3,ifile]=ma.masked;
	    
	    # fill table of ver. emittance, for later computation of time average and sigma
	    s,ind1,ind2=intersect(slots,sv[:,0]);
	    emitv[ind1,ifile]=sv[ind2,1];
	    # mask whenever the values are missing
	    ind3=complementary(ind1,np.arange(len(slots)));
	    if len(ind3)>0: emitv[ind3,ifile]=ma.masked;
	
    # write average emittances (over bunches) and their sigma to a Timber-like file
    fileout='TIMBER_DATA_calculated_BSRT_fast_scans_average_emit_vs_time'+opt.OUT+'_'+beam+'.csv';
    os.system('rm -f '+fileout);
    rootvar='BSRT_fast_scan_'
    time=np.array(time);
    variables=[rootvar+'AVERAGE_EMIT_H',rootvar+'AVERAGE_EMIT_V',
    	rootvar+'SIGMA_EMIT_H',rootvar+'SIGMA_EMIT_V'];
    arrays=['avEMITh','avEMITv','sigEMITh','sigEMITv'];
    for ivar,var in enumerate(variables):
        # convert data to array
    	exec(arrays[ivar]+'=np.array('+arrays[ivar]+')');
	# write
	eval('write_Timber(time,'+arrays[ivar]+',fileout,var)');
    
    
    if opt.AVER:
    	# compute time average
	averh=ma.average(emith,axis=1);
	averv=ma.average(emitv,axis=1);
	sigh=np.sqrt(ma.var(emith,axis=1));
	sigv=np.sqrt(ma.var(emitv,axis=1));
	sh=np.hstack((slots[~averh.mask].reshape((-1,1)),averh[~averh.mask].data.reshape((-1,1)),sigh[~averh.mask].data.reshape((-1,1))));
	sv=np.hstack((slots[~averv.mask].reshape((-1,1)),averv[~averv.mask].data.reshape((-1,1)),sigv[~averv.mask].data.reshape((-1,1))));
	rootname='average_BSRT_'+listname[0][-28:-9].replace('.','-').replace(':','-')+'_'+listname[-1][-28:-9].replace('.','-').replace(':','-')+'_'+beam;
	write_ncol_file(rootname+'H.dat',sh,header="25ns_slot\taverage_hor_emit[mm.mrad]\tsigma_hor_emit[mm.mrad]");
	write_ncol_file(rootname+'V.dat',sv,header="25ns_slot\taverage_ver_emit[mm.mrad]\tsigma_ver_emit[mm.mrad]");
	

    sys.exit();
    
    
