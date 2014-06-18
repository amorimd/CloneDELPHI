#!/usr/bin/python2.6

import sys
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX/Sussix")
sys.path.append("/home/nmounet/Documents/PYTHON/SUSSIX")
from parser import *
from string import *
import numpy as np
import pylab,os,re
from SussixNM import *
from read_cfg import read_cfg
from plot_lib import set_fontsize,init_figure,end_figure,plot,make_movie,plot2D,cmap
from io_lib import list_files
from string_lib import takeout_common
from read_Headtail_prt import extract_Headtail_param,read_allprt_file
from numpy.fft import fft,fft2,fftshift

def parsse():
    parser = OptionParser()
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the plot (several -b options possible)",
                      metavar="BNUM", default=[],dest="BNUM")
    parser.add_option("-c", "--coupledbunch",action="store_true",
                      help="Specify if we do an analysis of the full train motion (movie & turn-by-turn FFT)",
                      metavar="COUPL", default=False,dest="COUPL")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _hdtl.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--keeptmp",action="store_true",
                      help="Specify if temporary files are kept (not ereased) after making the coupled-bunch movie (with -c option)",
                      metavar="KEEP", default=False,dest="KEEP")
    parser.add_option("-m", "--maxtrace",type=int,
                      help="Specify the maximum number of traces to plot",
                      metavar="MAXT", default=10000,dest="MAXT")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    (opt, args) = parser.parse_args()
    return opt, args


def zeropad_middle(t,data):

    # takes irregularly spaced 1D-array 't' (monotonically increasing) and the corresponding
    # 1D-data 'data', and zero-pad where there are large "holes" in the middle (defined at points where diff(t)>10*(t[1]-t[0]) )
    
    # find places with holes
    deltatini=t[1]-t[0];
    if (deltatini<=0): print "Pb in zeropad: t not monotonically increasing !";sys.exit();
    
    # find where to pad with zeros in the middle
    indzero=pylab.mlab.find(np.diff(t)>10.*deltatini);
    # zero-padding in the middle
    for i in indzero:
    	tadd=np.arange(t[i],t[i+1],deltatini);lenadd=len(tadd);
        tnew=np.concatenate((t[:i+1],tadd,t[i+1:]));
        datanew=np.concatenate((data[:i+1],np.zeros(lenadd),data[i+1:]));
	
    return tnew,datanew;
    
def zeropad_end(t,data,nzero):

    # takes equidistant 1D-array 't' (monotonically increasing) and the corresponding
    # 1D-data 'data', and zero-pad at the end with 'nzero' zeros
    
    tnew=np.concatenate((t,np.arange(t[-1],t[-1]+nzero*(t[1]-t[0]),t[1]-t[0])));
    datanew=np.concatenate((data,np.zeros(nzero)));
    
    return tnew,datanew;
    

def zeropad_interp_equidistant(t,data,tbeg,tend,deltat,nzeroend=0):

    # takes irregularly spaced 1D-array 't' (monotonically increasing) and the corresponding
    # 1D-data 'data', zero-pad where there are large "holes" (defined at points where diff(t)>10*(t[1]-t[0]) )
    # and reinterpolate at equally spaced points, with spacing 'deltat'.
    # also add a number of zeros at the end (separated by 'deltat'), given by 'nzeroend'.
    
    # zero-padding in the middle (in the "holes")
    tpad,datapad=zeropad_middle(t,data)
    
    # interpolation on equidistant points seprated by deltat
    tint=np.arange(tbeg,tend,deltat);
    dataint=np.interp(tint,tpad,datapad);
    
    # zero-padding at the end
    tnew,datanew=zeropad_end(tint,dataint,nzeroend)
    
    return tnew,datanew;
    


if __name__ == "__main__":


    opt,args=parsse();
    
    bunches=np.array(opt.BNUM);
    bunches.sort();
    print "Selected bunches:", bunches

    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    # create list of associated legends (taking out from the list of file names 
    # all the common parameters in the names)
    listleg=takeout_common(listname);

    # loop on filenames
    for ifile,filename in enumerate(listname):
    
    	print filename;
    
	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-9])
	beta=np.sqrt(1.-1./gamma**2);
    
	# read from cfg file the number of slices, of kicks per turn, turn nubmers of first and last traces acquisition, and number of turns 
	# between each trace
	cfgfile=filename[:-9]+'.cfg';
	nsl=int(read_cfg(cfgfile,"Number_of_slices"));
	nini=int(read_cfg(cfgfile,"Start_turn_for_bunch_shape_acquisitions"));
	nfin=int(read_cfg(cfgfile,"Last_turn_for_bunch_shape_acquisitions"));
	nbet=int(read_cfg(cfgfile,"Number_of_turns_between_two_bunch_shape_acquisitions"));
	nsigmaz=float(read_cfg(cfgfile,"Longitud_extension_of_the_bunch"));
	
	print nsl,'slices, first trace at ',nini,'turns, last trace at ',nfin,'turns, ',nbet,'turns between traces, ',nsigmaz,'sigmas';

        #fil=filename.replace("_"," ").replace("_hdtl.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read file
	tau=[];xs=[];ys=[];
	for l,line in enumerate(open(filename)):
	    if (len(split(line))>0)and(len(tau)<opt.MAXT*nbunch*nsl):
        	tau.append(float(split(line)[0]))
        	xs.append(float(split(line)[1]))
		ys.append(float(split(line)[2]))
	
	# convert data to arrays
	tau=np.array(tau);
	xs=np.array(xs);
	ys=np.array(ys);

	for bnum in bunches:
	
	    print "bunch no ",bnum
	    
	    # initialize headtail plots (x and y)
	    figx,axx=init_figure();
	    figy,axy=init_figure();
		
	    
	    for iturn in range(len(xs[(bnum-1)*nsl::nbunch*nsl])):
	    
	        tau1=tau[iturn*nsl*nbunch+(bnum-1)*nsl:iturn*nsl*nbunch+(bnum-1)*nsl+nsl];
	        xs1=xs[iturn*nsl*nbunch+(bnum-1)*nsl:iturn*nsl*nbunch+(bnum-1)*nsl+nsl];
	        ys1=ys[iturn*nsl*nbunch+(bnum-1)*nsl:iturn*nsl*nbunch+(bnum-1)*nsl+nsl];

		# plot
		plot(tau1,xs1,'','-r',r'$n_S \cdot x_S$ [m] in the slices',axx,0,xlab=r'$\tau$ (ns) of the slices');
		plot(tau1,ys1,'','-r',r'$n_S \cdot y_S$ [m] in the slices',axy,0,xlab=r'$\tau$ (ns) of the slices');
		
		axx.set_title(fil+', bunch number '+str(bnum)+' , horizontal');
		axy.set_title(fil+', bunch number '+str(bnum)+' , vertical');
		
	    # finalize figures
	    if opt.SAVE:
		end_figure(figx,axx,save=filename[:-4]+'b'+str(bnum)+'_x');
		end_figure(figy,axy,save=filename[:-4]+'b'+str(bnum)+'_y');
	    else: 
		end_figure(figx,axx);
		end_figure(figy,axy);
	
	
	if opt.COUPL:
	
	    figcx,axcx=init_figure();
	    figcy,axcy=init_figure();
	    
	    deltat=0.01; # time interval for the FFT in ns
	    
	    # read prt file
	    x,y,epsx,epsy,invarx,invary,z,bl,zemit,nprleft_frac,bin_frac=read_allprt_file(filename[:-9]+'_prt.dat',nbunch);
	    
	    xsmax=np.max(np.abs(xs));ysmax=np.max(np.abs(ys));
	    
	    for iturn in range(len(xs[::nbunch*nsl])):
	    
	    	# find turn number where to look at in _prt.dat file, to find some additional info (z average and size of bunches)
		turnprt=(iturn*nbet+nini-1);
		print turnprt
		zave=z[nbunch*turnprt:nbunch*(turnprt+1)];
		szsize=bl[nbunch*turnprt:nbunch*(turnprt+1)];
		
		# time (ns) along the full train (adding the longitudinal offset of each bunch, in ns)
		tau2=tau[iturn*nsl*nbunch:(iturn+1)*nsl*nbunch]+np.reshape(np.transpose(np.ones((nsl,1))*zave/(beta*0.299792458)),(1,-1));
		tau2=np.squeeze(tau2)
	        xs2=xs[iturn*nsl*nbunch:(iturn+1)*nsl*nbunch];
	        ys2=ys[iturn*nsl*nbunch:(iturn+1)*nsl*nbunch];
		
		axcx.cla();axcy.cla();
		plot(tau2,xs2,'','-r',r'$n_S \cdot x_S$ [m] in the slices',axcx,0,xlab=r'$\tau$ (ns) of the slices');
		plot(tau2,ys2,'','-r',r'$n_S \cdot y_S$ [m] in the slices',axcy,0,xlab=r'$\tau$ (ns) of the slices');
		
		axcx.set_title(fil);axcy.set_title(fil);
		axcx.set_ylim([-xsmax,xsmax]);axcy.set_ylim([-ysmax,ysmax]);
	    
	    	fnamex=filename[:-4]+'_tmpx'+'%d'%iturn;
	    	fnamey=filename[:-4]+'_tmpy'+'%d'%iturn;
		
		end_figure(figcx,axcx,save=fnamex);
		end_figure(figcy,axcy,save=fnamey);
		
		# fft analysis
		if iturn==0: tbeg=tau2[-1];tend=tau2[0]; # beginning and end of times considered (fixed once for all)
		
		# zero-pad in the middle and interpolate on regular mesh
		taunew,xsnew=zeropad_interp_equidistant(tau2[::-1],xs2[::-1],tbeg,tend,deltat,nzeroend=0);
		taunew,ysnew=zeropad_interp_equidistant(tau2[::-1],ys2[::-1],tbeg,tend,deltat,nzeroend=0);
		# fft
		if iturn==0:
		    xfft1d=np.zeros((len(xs[::nbunch*nsl]),len(taunew)));
	    	    yfft1d=np.zeros((len(xs[::nbunch*nsl]),len(taunew)));

		xfft1d[iturn,:]=np.abs(fft(xsnew));
    		yfft1d[iturn,:]=np.abs(fft(ysnew));
		
	    # 2D plot of the turn by turn fft
	    fig2Dx,ax2Dx=init_figure();
	    fig2Dy,ax2Dy=init_figure();
	    plot2D(xfft1d,0,1e9/(tau2[0]-tau2[1]),nini,nfin,'Frequency [MHz]','Turn','1D FFT for each turn, '+fil+', horizontal',ax2Dx)		    
	    plot2D(yfft1d,0,1e9/(tau2[0]-tau2[1]),nini,nfin,'Frequency [MHz]','Turn','1D FFT for each turn, '+fil+', vertical',ax2Dy)		    
	    # finalize figures
	    if opt.SAVE:
		end_figure(fig2Dx,ax2Dx,save=filename[:-4]+'_turn_by_turn_fft_x');
		end_figure(fig2Dy,ax2Dy,save=filename[:-4]+'_turn_by_turn_fft_y');
	    else: 
		end_figure(fig2Dx,ax2Dx);
		end_figure(fig2Dy,ax2Dy);	    
	    
	    make_movie(filename[:-4]+'_x.gif',filename[:-4]+'_tmpx',flagrm=not(opt.KEEP));
	    make_movie(filename[:-4]+'_y.gif',filename[:-4]+'_tmpy',flagrm=not(opt.KEEP));
	    

    if (not(opt.SAVE))and(len(bunches)>0): pylab.show();
