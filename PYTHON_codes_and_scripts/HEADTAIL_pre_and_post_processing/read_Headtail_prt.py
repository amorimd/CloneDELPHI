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
import math
import matplotlib
from parser_lib import *
from plot_lib import cmap,plot,init_figure,end_figure,set_legend_fontsize
from io_lib import list_files
from string_lib import takeout_common,takeout_spaces
from read_cfg import read_cfg
from construct_bunchtable_from_scheme import check_nbunch


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--actions",action="store_true",
                      help="Specify if we also plot transverse actions",
                      metavar="ACT", default=False,dest="ACT")
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the raw data plot (several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-e", "--legend",action="append",
                      help="Specify the legend for the plot, for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-l", "--longitudinal",action="store_true",
                      help="Specify if we also plot longitudinal quantities (bunch length, average z, longitudinal emittance)",
                      metavar="LONG", default=False,dest="LONG")
    parser.add_option("-p", "--plotlog",type=int,
                      help="Specify if we want to plot in linear scale (0 - by default), in semilogx (1), in semilogy (2) or in loglog (3)",
                      metavar="LOG", default=0,dest="LOG")
    parser.add_option("-s", "--losses",action="store_true",
                      help="Specify if we also plot losses and fraction of binned macroparticles ",
                      metavar="LOSS", default=False,dest="LOSS")
    parser.add_option("-t", "--transverseemit",action="store_true",
                      help="Specify if we also plot transverse emittances",
                      metavar="EMIT", default=False,dest="EMIT")
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot the average",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--plotraw",type=int,
                      help="Specify some option for the plot of the raw data: 1 (default)=plot all turns vs. number of turns, n=plot every n turns vs. time",
                      metavar="PLOTRAW", default=1,dest="PLOTRAW")
    (opt, args) = parser.parse_args()
    print "Selected files:", opt.FILE
    print "Legends:", opt.LEG
    return opt, args


def extract_Headtail_param(root):

    # extract number of non-empty bunch, gamma, circumference,
    # and calculate Trev
    # also check that longitudinal matchin number is correct (i.e. near 1)
    # 'root' is the cfg file name without its extension '.cfg'
    
    # read from cfg file gamma and circumference
    cfgfile=root+'.cfg';
    gamma=float(read_cfg(cfgfile,"Relativistic_gamma"));
    circ=float(read_cfg(cfgfile,"Ring_circumference"));
    # revolution frequency and period (assumes protons)
    beta=np.sqrt(1.-1./(gamma*gamma));
    frev=beta*299792458./circ;
    Trev=1./frev;
    
    # check matching and find nb of non-empty bunches
    fh=open(root+"_inph.dat","r")
    if (len(fh.readlines())>0):
	fh.close();fh=open(root+"_inph.dat","r")
	for l in fh.readlines():
    	    if l.startswith("Longitudinal matching number"):
		ll=l.strip().split();
		match=float(ll[3])
    	    #if l.startswith("total number of macroparticles"):
	    #	ll=l.strip().split();
	    #	nbunch=int(ll[9].replace(',',''))

	if (np.abs(match-1.))>0.02: print "Warning: longitudinal matching number=",match
	
    #else:
    #print "Number of bunches unknown (empty _inph.dat file)."
    #nbunch=int(raw_input("Please type the number of bunches here: "))
    ibunchtb=int(read_cfg(cfgfile,"Switch_for_bunch_table"));
    if (ibunchtb==0) : nbunch=int(read_cfg(cfgfile,"Number_of_bunches"));
    else :
	bunchfile=root+'.bunch';
	fb=open(bunchfile,'r');
	nbunch,tmp=check_nbunch(fb);
	fb.close();

    fh.close()
    print 'Gamma=',gamma,', Circumference=',circ,', Trev=',Trev,', nbunch=',nbunch;
    
    
    return gamma,circ,Trev,nbunch


def check_data(length,nbunch):

    # check data: length of data should be a multiple of nbunch
    if ((length%nbunch)!=0):
        print 'Length of data is not a multiple of ',nbunch,'!'
	print 'Skipping last ',length%nbunch,' lines...'

    return length%nbunch
    

def read_allprt_file(filename,nbunch):

    # read prt file

    #x=[];xp=[];y=[];yp=[];bunch_nos=[];
    for l,line in enumerate(open(filename)):
	#try:
	#    i=int(split(line)[20])
	#except IndexError:
	if (len(split(line))<21):
	    print 'Line ',l+1,': not enough columns!'
	    l-=1;
	    break

    # check that the data is not somehow corrupted (missing / added lines, etc.)
    print "length of data",l+1
    lmod=check_data(l+1,nbunch)
    # lmod is the number of line to take away at the end (when length of data is not a multiple of nbunch)

    # initialize arrays
    x=np.zeros(l+1-lmod)
    y=np.zeros(l+1-lmod)
    epsx=np.zeros(l+1-lmod)
    epsy=np.zeros(l+1-lmod)
    invarx=np.zeros(l+1-lmod)
    invary=np.zeros(l+1-lmod)
    z=np.zeros(l+1-lmod)
    bl=np.zeros(l+1-lmod)
    zemit=np.zeros(l+1-lmod)
    nprleft_frac=np.zeros(l+1-lmod)
    bin_frac=np.zeros(l+1-lmod)

    # read file
    for il,line in enumerate(open(filename)):
    	if (il<(l+1-lmod)):
	    x[il]=float(split(line)[1])
	    y[il]=float(split(line)[3])
	    epsx[il]=float(split(line)[11])
	    epsy[il]=float(split(line)[12])
	    invarx[il]=float(split(line)[14])
	    invary[il]=float(split(line)[15])
	    z[il]=float(split(line)[5])
	    bl[il]=float(split(line)[9])
	    zemit[il]=float(split(line)[16])
	    nprleft_frac[il]=float(split(line)[18])
	    bin_frac[il]=float(split(line)[19])
	else: break
    
    #print len(x),len(y)
    
    return x,y,epsx,epsy,invarx,invary,z,bl,zemit,nprleft_frac,bin_frac
    

if __name__ == "__main__":
    opt,args=parsse(); 

    if (opt.BNUM==None): bunches=[];
    else: 
        bunches=np.array(opt.BNUM);
        bunches.sort();
    print "Selected bunches:", bunches

    # create list of filenames to analyse
    listname=list_files(opt.FILE);

    # create list of associated legends (either with opt.LEG or 
    # the list of file names taking out all the common parameters
    # in the names)
    if (opt.LEG!=None):
        listleg=opt.LEG;
    else:
        listleg=takeout_common(listname);
	listleg=takeout_spaces(listleg);
    
    # find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
    gamma,circ,Trev,nbunch=extract_Headtail_param(listname[0][:-8])
    
    
    if (len(bunches)>0):
	# initialize plots
	figx,axx=init_figure();
	figy,axy=init_figure();
	if (opt.EMIT):
	    figepsx,axepsx=init_figure();
	    figepsy,axepsy=init_figure();
	if (opt.ACT):
	    figinvx,axinvx=init_figure();
	    figinvy,axinvy=init_figure();
	if (opt.LONG):
	    figz,axz=init_figure();
	    figbl,axbl=init_figure();
	    figzemit,axzemit=init_figure();
	if (opt.LOSS):
	    fignpr,axnpr=init_figure();
	    figbin,axbin=init_figure();
    else: axx=0;axy=0;axepsx=0;axepsy=0;axz=0;axbl=0;

    if (opt.AVER):
        figbeamx,axbeamx=init_figure();
        figbeamy,axbeamy=init_figure();
	if (opt.EMIT):
	    figepsbeamx,axepsbeamx=init_figure();
	    figepsbeamy,axepsbeamy=init_figure();
	if (opt.ACT):
	    figinvbeamx,axinvbeamx=init_figure();
	    figinvbeamy,axinvbeamy=init_figure();
	if (opt.LONG):
	    figbeamz,axbeamz=init_figure();
	    figbeambl,axbeambl=init_figure();
	    figbeamzemit,axbeamzemit=init_figure();
	if (opt.LOSS):
	    figbeamnpr,axbeamnpr=init_figure();
	    figbeambin,axbeambin=init_figure();

    # loop on filenames
    for ifile,filename in enumerate(listname):
    

	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    

    	# for the legend
    	#fil=filename.replace("_"," ").replace("prt.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read prt file
	x,y,epsx,epsy,invarx,invary,z,bl,zemit,nprleft_frac,bin_frac=read_allprt_file(filename,nbunch)

	# choose bunches to iterate on
	if (opt.AVER): bunches_iter=range(nbunch,0,-1);
	else: bunches_iter=bunches;
	
        # extract data and plot each bunch chosen	
	for bnum in bunches_iter:
	
	    x1=x[bnum-1::nbunch];
    	    y1=y[bnum-1::nbunch];
	    epsx1=epsx[bnum-1::nbunch];
    	    epsy1=epsy[bnum-1::nbunch];
	    invarx1=invarx[bnum-1::nbunch];
    	    invary1=invary[bnum-1::nbunch];
	    z1=z[bnum-1::nbunch];
    	    bl1=bl[bnum-1::nbunch];
	    zemit1=zemit[bnum-1::nbunch];
	    nprleft_frac1=nprleft_frac[bnum-1::nbunch];
    	    bin_frac1=bin_frac[bnum-1::nbunch];
	    turns=np.arange(0,len(x1),dtype=float);
	    
    	    if (bnum in bunches):
		# plot raw data (average positions)
    		plot(turns,x1,fil+', bunch '+str(bnum),'.',"Average x position of the bunch [m]",axx,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,y1,fil+', bunch '+str(bnum),'.',"Average y position of the bunch [m]",axy,opt.LOG,plotevery=opt.PLOTRAW);
		if (opt.EMIT):
		    # plot emittances
    		    plot(turns,epsx1,fil+', bunch '+str(bnum),'-',"Normalized horizontal emittance of the bunch [mm.mrad]",axepsx,opt.LOG,plotevery=opt.PLOTRAW);
    		    plot(turns,epsy1,fil+', bunch '+str(bnum),'-',"Normalized vertical emittance of the bunch [mm.mrad]",axepsy,opt.LOG,plotevery=opt.PLOTRAW);
		if (opt.ACT):
		    # plot invariants ("actions")
    		    plot(turns,invarx1,fil+', bunch '+str(bnum),'-',"Average horizontal invariant of the bunch [m]",axinvx,opt.LOG,plotevery=opt.PLOTRAW);
    		    plot(turns,invary1,fil+', bunch '+str(bnum),'-',"Average vertical invariant of the bunch [m]",axinvy,opt.LOG,plotevery=opt.PLOTRAW);
		if (opt.LONG):
		    # plot average z, bunch length and longitudinal emittance
    		    plot(turns,z1,fil+', bunch '+str(bnum),'-',"Average z position of the bunch [m]",axz,opt.LOG,plotevery=opt.PLOTRAW);
    		    plot(turns,bl1,fil+', bunch '+str(bnum),'-',"Bunch length [m]",axbl,opt.LOG,plotevery=opt.PLOTRAW);
    		    plot(turns,zemit1,fil+', bunch '+str(bnum),'-',"Longitudinal emittance of the bunch [eV.s]",axzemit,opt.LOG,plotevery=opt.PLOTRAW);
		if (opt.LOSS):
		    # plot fraction of macroparticles left and fraction of binned macroparticles
    		    plot(turns,nprleft_frac1,fil+', bunch '+str(bnum),'-',"Fraction of macroparticles still in the bunch",axnpr,opt.LOG,plotevery=opt.PLOTRAW);
    		    plot(turns,bin_frac1,fil+', bunch '+str(bnum),'-',"Fraction of binned macroparticles in the bunch",axbin,opt.LOG,plotevery=opt.PLOTRAW);

	    if (opt.AVER):
	    	if (bnum==nbunch):
		    xave=x1;yave=y1;
		    epsxave=epsx1;epsyave=epsy1;
		    invarxave=invarx1;invaryave=invary1;
		    zave=z1;blave=bl1;zemitave=zemit1;
		    nprleft_frac_ave=nprleft_frac1;
		    bin_frac_ave=bin_frac1;
		else:
		    xave+=x1;yave+=y1;	    
		    epsxave+=epsx1;epsyave+=epsy1;
		    invarxave+=invarx1;invaryave+=invary1;
		    zave+=z1;blave+=bl1;zemitave+=zemit1;
		    nprleft_frac_ave+=nprleft_frac1;
		    bin_frac_ave+=bin_frac1;

		
	if (opt.AVER):
	    xave/=float(nbunch);yave/=float(nbunch);
	    epsxave/=float(nbunch);epsyave/=float(nbunch);
	    invarxave/=float(nbunch);invaryave/=float(nbunch);
	    zave/=float(nbunch);blave/=float(nbunch);zemitave/=float(nbunch);
	    nprleft_frac_ave/=float(nbunch);
	    bin_frac_ave/=float(nbunch);
	    turns=np.arange(0,len(xave),dtype=float);
	    # plot raw data (average positions)
    	    plot(turns,xave,fil,'.',"Average x position of the beam [m]",axbeamx,opt.LOG,plotevery=opt.PLOTRAW);
    	    plot(turns,yave,fil,'.',"Average y position of the beam [m]",axbeamy,opt.LOG,plotevery=opt.PLOTRAW);
	    if (opt.EMIT):
		# plot emittances
    		plot(turns,epsxave,fil,'-',"Average normalized horizontal emittance of the beam [mm.mrad]",axepsbeamx,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,epsyave,fil,'-',"Average normalized vertical emittance of the beam [mm.mrad]",axepsbeamy,opt.LOG,plotevery=opt.PLOTRAW);
	    if (opt.ACT):
		# plot invariants ("actions")
    		plot(turns,invarxave,fil,'-',"Average horizontal invariant of the beam [mm.mrad]",axinvbeamx,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,invaryave,fil,'-',"Average vertical invariant of the beam [mm.mrad]",axinvbeamy,opt.LOG,plotevery=opt.PLOTRAW);
	    if (opt.LONG):
		# plot average z, bunch length and longitudinal emittance
    		plot(turns,zave,fil,'-',"Average z position of the beam [m]",axbeamz,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,blave,fil,'-',"Average bunch length of the beam [m]",axbeambl,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,zemitave,fil,'-',"Average longitudinal emittance of the beam [m]",axbeamzemit,opt.LOG,plotevery=opt.PLOTRAW);
 	    if (opt.LOSS):
		# plot fraction of macroparticles left and fraction of binned macroparticles
    		plot(turns,nprleft_frac_ave,fil,'-',"Fraction of macroparticles still in the beam",axbeamnpr,opt.LOG,plotevery=opt.PLOTRAW);
    		plot(turns,bin_frac_ave,fil,'-',"Fraction of binned macroparticles in the beam",axbeambin,opt.LOG,plotevery=opt.PLOTRAW);
	    
	    


    if (len(bunches)>0):
	set_legend_fontsize(figx,axx);
	set_legend_fontsize(figy,axy);
	if (opt.EMIT):
	    set_legend_fontsize(figepsx,axepsx);
	    set_legend_fontsize(figepsy,axepsy);
	if (opt.ACT):
	    set_legend_fontsize(figinvx,axinvx);
	    set_legend_fontsize(figinvy,axinvy);
	if (opt.LONG):
	    set_legend_fontsize(figz,axz);
	    set_legend_fontsize(figbl,axbl);
	    set_legend_fontsize(figzemit,axzemit);
	if (opt.LOSS):
	    set_legend_fontsize(fignpr,axnpr);
	    set_legend_fontsize(figbin,axbin);
	
    if (opt.AVER):
	set_legend_fontsize(figbeamx,axbeamx);
	set_legend_fontsize(figbeamy,axbeamy);
	if (opt.EMIT):
	    set_legend_fontsize(figepsbeamx,axepsbeamx);
	    set_legend_fontsize(figepsbeamy,axepsbeamy);
	if (opt.ACT):
	    set_legend_fontsize(figinvbeamx,axinvbeamx);
	    set_legend_fontsize(figinvbeamy,axinvbeamy);
	if (opt.LONG):
	    set_legend_fontsize(figbeamz,axbeamz);
	    set_legend_fontsize(figbeambl,axbeambl);
	    set_legend_fontsize(figbeamzemit,axbeamzemit);
	if (opt.LOSS):
	    set_legend_fontsize(figbeamnpr,axbeamnpr);
	    set_legend_fontsize(figbeambin,axbeambin);

    pylab.show();

    sys.exit()

