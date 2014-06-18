#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
import pylab,re,random,pytz
from datetime import time,datetime,date
import numpy as np
from string import split, replace
from optparse import OptionParser
import math
from SussixNM import *
from plot_lib import make_movie,init_figure,cmap,plot,end_figure
from io_lib import list_files,read_ncol_file,write_ncol_file
from read_Headtail_prt_fit import slideaver
from read_Headtail_prt_sussix import extract
from read_Headtail_prt_coupledbunch import plot_spec,gettune


def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--fbct",
                      help="Specify the FBCT filename, to get filling scheme and find the correct slot numbers for ADT data. This should be a two column (at least) file, with a single line header, as obtained by read_FBCT.py (average intensities vs non-empty 25ns-slot numbers)",
                      metavar="FBCT", default=None,dest="FBCT")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the ADT (.csv) name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--nk",type=int,
                      help="Specify the numbers of highest lambdas to consider for the SVD (default=1)",
                      metavar="NK", default=1,dest="NK")
    parser.add_option("-l", "--limitsbad",type=int,nargs=2,
                      help="Specify limits (2 indices) of the bad data (determined visually - can be rather approximate). To be used when the algorithm to find the bad (dump) data fails. Default= determined automatically.",
                      metavar="LIM", default=None,dest="LIM")		      
    parser.add_option("-n", "--nsigoffset",type=float,
                      help="Specify opposite of offset (in number of sigma of the full signal) from the full signal average, below which the sliding average should be to define the 'bad' (i.e. dump) data. Can be negative (if bad data above the average) (default=2).",
                      metavar="OFF", default=2,dest="OFF")
    parser.add_option("-o", "--output",
    		      help="Specify some suffix for the output filenames (empty string by default)",
                      default="",dest="OUT")
    parser.add_option("-q", "--tune",type=float,
                      help="Specify the factional part of the tune (default=0.31)",
                      metavar="TUNE", default=0.31,dest="TUNE")
    parser.add_option("-s", "--slideaver",type=int,
                      help="Specify period used for the sliding average used to find the beginning of the data (default=20).",
                      metavar="SLIDE", default=20,dest="SLIDE")		      
    parser.add_option("-w", "--width",type=float,
                      help="Specify the half width of the window into which we should find the tune (default=0.05)",
                      metavar="WIDTH", default=0.05,dest="WIDTH")
    parser.add_option("-x", "--xlim",type=int,nargs=2,
                      help="Specify limits for the x axis, for the plots vs 25-ns slot numbers (default=[0,3654]).",
                      metavar="XLIM", default=[0,3564],dest="XLIM")		      
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    (opt, args) = parser.parse_args()
    #print "Selected files:", opt.FILE
    return opt, args


def find_beginning_postmortem(data,slide,nsigoffset,nslots=3564,FBCT=None,optlim=None):

    # find the correct beginning of the postmortem data.
    # If FBCT is not None, it should contain a 1D array of size 3564
    # giving the intensity (in p+) of each 25ns-slot. Then it tries to find
    # the best overlap with the FBCT.
    # if optlim is not None, it contains two numbers which are the indices
    # of resp. the beginning and the end of the 'bad' data. Then the part
    # of the algorithm that tries to find where is the bad data is bypassed.
    
    # compute global sigma
    sigma=np.sqrt(np.var(data));#print sigma
    
    # do a sliding average
    dataslideav=slideaver(data,slide);
    
    # compute average of this sliding average
    aver=np.average(dataslideav);#print aver
    
    indbeg=0;indfinal=len(data);
    
    if (optlim==None):
	# extraction of indices where sliding average is much lower or much higher
	# than the rest (by nsigoffset number of sigma)
	if (nsigoffset>0): ind=pylab.mlab.find(dataslideav<aver-nsigoffset*sigma);
	else: ind=pylab.mlab.find(dataslideav>aver-nsigoffset*sigma);
    else:
    	# 'bad' data is between the numbers in optlim
	ind=optlim;
    
    if (len(ind)>0):
    	# take the last of those indices and look for closest place where
	# sliding average is 0
	ind2=pylab.mlab.find(dataslideav[ind[-1]:]==0)+ind[-1];
	if (len(ind2)>0): 
	    indbeg=ind2[0];
	    if (FBCT!=None):
		# find best overlap with FBCT intensities
		# replace first data by array of 0 and 1
		data2=np.zeros(len(data));
		inddata=pylab.mlab.find(data>0.01);
		data2[inddata]=1;
		overlap=np.array([np.sum(FBCT*np.abs(data2[indbeg+i:indbeg+i+nslots])) for i in range(nslots)]);
		imax=np.argmax(overlap);
		indbeg+=imax;
		#pylab.plot(np.arange(nslots),overlap);
    	    
	    # take the last of those indices and look for closest place where
	    # sliding average is 0
	    ind3=pylab.mlab.find(dataslideav[:ind[0]+1]==0);
	    if (len(ind3)>0): indfinal=ind3[-1];
	    
    datanew=np.concatenate((data[indbeg:],data[:indfinal]))
	
    return datanew,dataslideav,indbeg,indfinal


def transform_FBCT_intensities(data,nslots=3564):
    # take a two columns FBCT array 'data' giving (slot number, intensity)
    # (ordered) and transform into an array of length 'nslot' such that
    # newdata[slot nb]=1 if there is a 'normal' intensity bunch, 0 otherwise
    # (newdata array has many zeros)
    
    newdata=np.zeros(nslots);
    
    for islot,slot in enumerate(data[:,0]):
    
    	newdata[slot]=(data[islot,1]>0.01);
    
    return newdata;
    

if __name__ == "__main__":
    opt,args=parsse(); 

    nslots=3564;
    
    if (opt.FBCT!=None):
        fbct=read_ncol_file(opt.FBCT,ignored_rows=1);
	dataFBCT=transform_FBCT_intensities(fbct,nslots=nslots);
	#for islot in range(nslots): print islot,dataFBCT[islot]
    else:
    	dataFBCT=None;
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    for ifile,filename in enumerate(listname):
    
    	# string for the title
        fil=filename.replace("_"," ").replace(".csv","").replace("/"," ").replace("-"," ");
    	print "Selected file:", filename
	
    	# read file
	s=read_ncol_file(filename,ignored_rows=2,delimiter=',');

	# plot raw data with and without pretreatment, showing the limits found
	fig,ax=init_figure();
	# plot initial signal
	plot(np.arange(len(s[:,0])),s[:,0],'raw data','b','ADT pickup signal [a.u.]',ax,0,xlab='');
	#pylab.show()
	
	# pre-treatment to find first turn
	data,dataslideav,indbeg,indfinal=find_beginning_postmortem(s[:,0],opt.SLIDE,opt.OFF,nslots=nslots,FBCT=dataFBCT,optlim=opt.LIM);
	
	# plot sliding average
	plot(np.arange(len(dataslideav)),dataslideav,'sliding average of raw data','r','ADT pickup signal [a.u.]',ax,0,xlab='');
	# plot limits of good signal
	smax=np.max(np.abs(s[:,0]));
	plot([indbeg-0.01,indbeg+0.01],[-smax,smax],'beginning of signal','--g','ADT pickup signal [a.u.]',ax,0,lw=5.,xlab='');
	plot([indfinal-0.01,indfinal+0.01],[-smax,smax],'end of signal','--m','ADT pickup signal [a.u.]',ax,0,lw=5.,xlab='');
	# plot good signal on top of the initial one
	plot(indbeg+np.arange(len(data)),data,'raw data after pretreatment','.k','ADT pickup signal [a.u.]',ax,0,xlab='');
	if (opt.XLIM!=[0,3564]): ax.set_xlim(indbeg+opt.XLIM);
	if (opt.FBCT!=None):
	    # plot FBCT on top of it
	    ax2=pylab.twinx(ax=ax); # second y axis
	    plot(indbeg+np.arange(nslots,dtype=float),dataFBCT,'Filling scheme','.c','Filling scheme: empty (0) or filled (1)',ax2,0,lw=3.)
	    if (opt.XLIM!=[0,3564]): ax2.set_xlim(indbeg+opt.XLIM);
	    ax2.legend(loc=2)
	
	end_figure(fig,ax,save=opt.SAVE*(filename.replace(".csv","")+"_raw"+opt.OUT),legpos=(opt.FBCT!=None));
	#pylab.show()
	
	# reshape data
	print len(data)
	data=data[:len(data)-mod(len(data),nslots)];
	x2d=data.reshape((-1,nslots));
	
	# make a movie with raw turn-by-turn and bunch-by-bunch (actually slot-by-slot) data
	print len(x2d[:,0])
	for turn in range(len(x2d[:,0])):
	    fig,ax=init_figure();
	    plot(np.arange(nslots),x2d[turn,:],'turn '+str(turn),'b','ADT pickup signal [a.u.]',ax,0,lw=2.5,xlab='25-ns slot number')
	    ax.set_xlim(opt.XLIM);
	    ax.set_ylim([-np.ceil(smax/10)*10,np.ceil(smax/10)*10]);
	    end_figure(fig,ax,legpos=1,save='_tmp'+str(turn));
	
	make_movie(filename.replace('.csv','')+'_raw_turnbyturn_movie'+opt.OUT+'.gif','_tmp',flagrm=True);


	# Reduced SVD (see "BEAM OBSERVATIONS WITH ELECTRON CLOUD IN THE CERN PS & SPS COMPLEX",
	# G. Arduini et al, Proceedings ECLOUD 2004, p. 31, CERN-2005-001)
	u,lambd,v=np.linalg.svd(x2d,full_matrices=False)

	# extract maximum lambdas
	#k=np.argmax(np.abs(lambd))
	indk=np.argsort(np.abs(lambd));

	tumost=np.zeros(opt.NK);
	
	filetune=open(filename.replace('.csv','')+'_SVD_time_tune'+opt.OUT+'.txt','w');
	print >> filetune, "SVD_mode_ik\tlambda\ttune";

	# analyse each of the opt.NK highest "mode"
	for ik,k in enumerate(indk[-1:-opt.NK-1:-1]):

	    # figures initialization
	    figSVDspatial,axSVDspatial=init_figure();
	    figSVDtime,axSVDtime=init_figure();

	    # plot spatial pattern of the maximum lambda
	    plot(np.arange(nslots),v[k,:],'Spatial pattern','-b','Spatial pattern [a.u]',axSVDspatial,0,lw=4.,xlab='25-ns slot number');
	    axSVDspatial.set_xlim(opt.XLIM);

	    # write SVD spatial pattern in a file
	    fileSVD=filename.replace(".csv","")+'_SVD_modek'+str(ik)+'_spatial'+opt.OUT+'.dat';
	    dataSVD=np.hstack((np.arange(nslots).reshape((-1,1)),v[k,:].reshape((-1,1))));
	    write_ncol_file(fileSVD,dataSVD,header="25ns_slot\tSVD_spatial");
	    
	    if (opt.FBCT!=None):
	        # plot FBCT filling scheme on top of SVD spatial pattern
	        axSVDspatial2=pylab.twinx(ax=axSVDspatial); # second y axis
	    	#plot(np.arange(nslots,dtype=float),dataFBCT,'FBCT','-r','FBCT [p+/bunch]',axSVDspatial2,0,lw=4.)
	    	plot(np.arange(nslots,dtype=float),dataFBCT,'Filling scheme','.c','Filling scheme: empty (0) or filled (1)',axSVDspatial2,0,lw=4.)
	    	axSVDspatial2.set_xlim(opt.XLIM);
	    	axSVDspatial2.legend(loc=2);
	    
	    #axSVDspatial.set_title(SVD, '+fil+', ik='+str(ik));
	    # plot time pattern of the maximum lambda
	    plot(np.arange(len(u[:,k]),dtype=float),u[:,k],'','-b','Time pattern [a.u]',axSVDtime,0,lw=4.)
	    #axSVDtime.set_title('SVD, '+fil+', ik='+str(ik))

	    end_figure(figSVDspatial,axSVDspatial,save=opt.SAVE*(filename.replace(".csv","")+'_SVD_modek'+str(ik)+'_spatial'+opt.OUT),legpos=(opt.FBCT!=None));
            end_figure(figSVDtime,axSVDtime,save=opt.SAVE*(filename.replace(".csv","")+'_SVD_modek'+str(ik)+'_time'+opt.OUT));

	    # figure initialization
	    figSVDtune,axSVDtune=init_figure();

	    # plot spectrum of the time pattern
	    suss=gettune(u[:,k],np.zeros(len(u[:,k])),opt.TUNE,0.,0.07)
	    # sort and extract only tunes at 'width' from nominal tune Q
	    tunes,amps=extract(suss.ox,suss.ax,opt.TUNE,opt.WIDTH)
	    plot_spec(tunes,amps,'','-',axSVDtune);
	    #axSVDtune.set_title('Time pattern spectrum, '+fil+', ik='+str(ik));
	    end_figure(figSVDtune,axSVDtune,save=opt.SAVE*(filename.replace(".csv","")+'_SVD_modek'+str(ik)+'_time_spectrum'+opt.OUT));

	    tumost[ik]=tunes[np.argmax(amps)];
	    print >> filetune, ik,lambd[k],tumost[ik];
		
    	filetune.close();
    
    if not(opt.SAVE): pylab.show();

    sys.exit()
