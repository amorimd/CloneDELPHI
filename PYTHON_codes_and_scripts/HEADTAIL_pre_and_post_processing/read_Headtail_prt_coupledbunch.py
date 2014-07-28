#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from optparse import OptionParser
from string import *
import pylab,os,re
from SussixNM import *
from read_cfg import read_cfg
from plot_lib import init_figure,end_figure,plot,plot2D
from io_lib import list_files
from string_lib import takeout_common
from read_Headtail_prt import extract_Headtail_param
from read_Headtail_prt_fit import read_prt_file
from read_Headtail_prt_laplace import find_peak_lapl,Laplace
from numpy.fft import fft,fft2,fftshift
import numpy as np


def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--fftanalysis",action="store_true",
                      help="Specify if we perform 1D FFT analysis (both in bunches and in turns)",
                      metavar="FFT1", default=False,dest="FFT1")
    parser.add_option("-d", "--fft2danalysis",action="store_true",
                      help="Specify if we perform 2D FFT analysis",
                      metavar="FFT2", default=False,dest="FFT2")
    parser.add_option("-e", "--legend",action="append",
                      help="Specify the legend for the plot, for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--nk",type=int,
                      help="Specify the numbers of highest lambdas to consider for the SVD (default=1)",
                      metavar="NK", default=1,dest="NK")
    parser.add_option("-n", "--nmodes",type=int,
                      help="Specify number of modes to extract with Laplace (default=20)",
                      metavar="NMOD", default=20,dest="NMOD")
    parser.add_option("-o", "--output",help="Specify output suffix for mode number and tune of highest mode (default=\"mode2D.txt\")",
                      default="mode2D.txt",dest="OUT")
    parser.add_option("-p", "--spectrumsvd",action="store_true",
                      help="Specify if we perform spectral analysis (Sussix and Laplace) on the SVD time pattern",
                      metavar="SPECSVD", default=False,dest="SPECSVD")
    parser.add_option("-r", "--raw",action="store_true",
                      help="Specify if we also plot raw data (1D for the last turn, and 2D)",
                      metavar="RAW", default=False,dest="RAW")
    parser.add_option("-s", "--svd",action="store_true",
                      help="Specify if we perform Singular Value Decomposition (time domain analysis)",
                      metavar="SVD", default=False,dest="SVD")
    parser.add_option("-t", "--turns",type=int,nargs=2,
                      help="Specify the first and last turn (+1) between which we plot and analyse",
                      metavar="TURNS", default=[0,10000000],dest="TURNS")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    (opt, args) = parser.parse_args()
    return opt, args


def gettune(x,y,Qx,Qy,Qs,nharm=300):

    # get tunes and amplitudes from Sussix (case without px and py information)
    nturns=len(x);
    px=np.zeros(nturns)
    py=np.zeros(nturns)
    sussix_inp(ir=1, turns=nturns, tunex=Qx, tuney=Qy, tunez=Qs, istun=0.02, idam=2, narm=nharm)
    a=sussixBS(x,y,px,py)
    
    return a


def find_max(data2d,Q,width):

    # find max in 2d data, in a window between Q-width and Q+width for the second dimension
    tunes=np.arange(0,data2d.shape[1],dtype=float)/float(data2d.shape[1]);
    m=np.searchsorted(tunes,Q-width)
    p=np.searchsorted(tunes,Q+width,side='right')
    newtunes=tunes[m:p]
    newdata2d=data2d[:,m:p]
    datamax=np.unravel_index(np.argmax(newdata2d),newdata2d.shape)
    datamax1 = (datamax[0], newtunes[datamax[1]])
    #datamax1 = (datamax[0], datamax[1])
    return datamax1


def plot_bunch(nbunch,data,leg,col,ylab,tit,ax,lw=2.5,bunchoffset=0):
    
    # plot vs bunch numbers
    ax.plot(np.arange(bunchoffset,bunchoffset+nbunch),data,col,label=leg,linewidth=lw);

    #ax.set_xlabel("Bunch number (0 is the head of the bunch train)")
    ax.set_xlabel("Bunch number")
    ax.set_ylabel(ylab)
    ax.legend(loc=0)
    #ax.set_title(tit)


def plot_spec(x,data,title,col,ax,lw=1.):

    #lab='x' or 'y'
    # plotting spectrum
    ax.semilogy(x,data,col,linewidth=lw);
    ax.set_xlabel("Tune")
    ax.set_ylabel("Spectrum amplitude [a.u.]")
    ax.set_title(title)
    #ax.legend(loc=0)


def fft1d_bunch(x,nbunch,Q,width,fil,plane,save=None):

    # 1D FFT for each bunch
    # plane='x' or 'y'

    # initialization
    fig,ax=init_figure();
    xfft1d=np.zeros((nbunch,len(x[::nbunch])));

    # extract and do the fft for each bunch
    for bnum in range(0,nbunch):
	x1d=x[bnum::nbunch];
	xfft1d[bnum,:]=np.abs(fft(x1d))

    if (plane=='x'): plan='Horizontal'
    else: plan='Vertical'

    # plot
    plot2D(xfft1d,0,1-1./float(len(x[::nbunch])),0,nbunch-1,plan+' tune','Bunch number (0 is the head of the bunch train)','1D FFT for each bunch, '+fil+', '+plane,ax)		    
    ax.set_xlim([Q-width,Q+width])

    if (save==None):
        end_figure(fig,ax);
    else:
    	end_figure(fig,ax,save=save+'_bunch_by_bunch_fft');
  

def fft1d_turn(x,nbunch,firstturn,fil,plane,save=None):

    # 1D FFT for each turn (along bunch train)
    # initialization
    fig,ax=init_figure();
    xfft1d=np.zeros((nbunch,len(x[::nbunch])));

    # extract and do the fft for each turn
    for turn in range(0,len(x[::nbunch])):
	x1d=x[turn*nbunch:(turn+1)*nbunch];
	xfft1d[:,turn]=np.abs(fft(x1d))

    if (plane=='x'): plan='Horizontal'
    else: plan='Vertical'

    xmax=find_max(xfft1d,0,1);
    print plan," plane, highest mode: ", xmax[0];
    
    # plot
    plot2D(xfft1d,firstturn,firstturn+len(x[::nbunch])-1,0,nbunch-1,'Turn','Mode number','1D FFT for each turn, '+fil+', '+plane,ax)		    

    if (save==None):
    	end_figure(fig,ax);
    else:
    	end_figure(fig,ax,save=save+'_turn_by_turn_fft');


def fft2d(x2d,Q,width,fil,plane,save=None):

    # 2D FFT
    # initialization
    fig,ax=init_figure()
    nbunch=x2d.shape[1]

    xfft2d=np.transpose(np.abs(fft2(x2d)))

    # find maxima at 0.2 from initial tunes
    xmax=find_max(xfft2d,Q,0.2)
    #print xfft2d[:,xmax[1]]

    if (plane=='x'): plan='Horizontal'
    else: plan='Vertical'

    # plot
    plot2D(xfft2d[:,:],0,1-1./float(x2d.shape[0]),0,nbunch-1,plan+' tune','Mode number','2D FFT, '+fil+', '+plane,ax)		    
    ax.set_xlim([Q-width,Q+width])

    if (save==None):
    	end_figure(fig,ax);
    else:
    	end_figure(fig,ax,save=save+'_2Dfft');
    
    return xmax


def svd(x2d,Q,width,firstturn,lastturn,fil,plane,Qs=0.,nk=1,flagspec=1,figSVDspatial=None,axSVDspatial=None,bunchoffset=0,nmodes=20,save=None):

    # SVD (time domain analysis)
    from read_Headtail_prt_sussix import extract
    
    nbunch=x2d.shape[1]

    if (plane=='x'):
        plan='Horizontal';col='b';
    else:
        plan='Vertical';col='r';

    if (x2d.shape[0]<nbunch):
	print "SVD error: numbers of turns needs to be higher than number of bunches !"
	sys.exit()

    # Reduced SVD (see "BEAM OBSERVATIONS WITH ELECTRON CLOUD IN THE CERN PS & SPS COMPLEX",
    # G. Arduini et al, Proceedings ECLOUD 2004, p. 31, CERN-2005-001)
    u,lambd,v=np.linalg.svd(x2d[firstturn:lastturn,:],full_matrices=False)

    # extract maximum lambdas
    #k=np.argmax(np.abs(lambd))
    indk=np.argsort(np.abs(lambd));
    
    tumost=np.zeros((nk,3,2));
    if (Qs==0.):
    	mmin=0;mmax=0;
    	tuneshift=np.zeros((nk,1));tau=np.zeros((nk,1));amp=np.zeros((nk,1));
    else:
	mmin=-1;mmax=1;
    	tuneshift=np.zeros((nk,3));tau=np.zeros((nk,3));amp=np.zeros((nk,3));

    # analyse each of the nk highest "mode"
    for ik,k in enumerate(indk[-1:-nk-1:-1]):
	    
	# figure initialization
	if (axSVDspatial==None)or(figSVDspatial==None): figSVDspatial,axSVDspatial=init_figure();
	figSVDtime,axSVDtime=init_figure();

	# plot spatial pattern of the maximum lambda
	plot_bunch(nbunch,v[k,:],plan,col,'Spatial pattern [a.u]','SVD, '+fil+', '+plane+', ik='+str(ik),axSVDspatial,bunchoffset=bunchoffset)

	# plot time pattern of the maximum lambda
	plot(np.arange(firstturn,lastturn,dtype=float),u[:,k],'','-b',plan+' time pattern [a.u]',axSVDtime,0,lw=2.5)
	#axSVDtime.set_title('SVD, '+fil+', '+plane+', ik='+str(ik))
        
	if (save==None):
	    end_figure(figSVDspatial,axSVDspatial);
            end_figure(figSVDtime,axSVDtime);
	else:
	    end_figure(figSVDspatial,axSVDspatial,save=save+'_SVD_modek'+str(ik)+'_spatial');
            end_figure(figSVDtime,axSVDtime,save=save+'_SVD_modek'+str(ik)+'_time');
	    

	if (flagspec):
	
	    # figure initialization
	    figSVDtune,axSVDtune=init_figure();
	    figSVDlapl,axSVDlapl=init_figure();
	    
	    # plot spectrum of the time pattern
	    suss=gettune(u[:,k],np.zeros(len(u[:,k])),Q,0.,0.07)
	    # sort and extract only tunes at 'width' from nominal tune Q
	    tunes,amps=extract(suss.ox,suss.ax,Q,width)
	    plot_spec(tunes,amps,'Time pattern spectrum, '+fil+', '+plane+', ik='+str(ik),'-',axSVDtune);

	    # complex frequency spectral analysis of the time pattern
	    tunes,rates,amps=Laplace(u[:,k],nmodes);

	    # find spectral lines associated with headtail modes
	    for j,m in enumerate(range(mmin,mmax+1)):
		tuneshift[ik,j],ra,amp[ik,j]=find_peak_lapl(tunes,rates,amps,Q,Qs,m);
		if (ra!=0): tau[ik,j]=1./ra;
	    
	    # find spectral line of highest amplitude overall (within width/2 from the tune)
	    tumost[ik,0,0],tumost[ik,1,0],tumost[ik,2,0]=find_peak_lapl(tunes,rates,amps,Q,width,0)
	    # find spectral line of highest growth rate overall (within width/2 from the tune)
	    # (to do this we just need to invert 'amps' and 'rates' in find_peak_lapl)
	    tumost[ik,0,1],tumost[ik,2,1],tumost[ik,1,1]=find_peak_lapl(tunes,amps,rates,Q,width,0)

	    # plot Laplace spectra
	    axSVDlapl.plot(tunes,np.abs(np.array(amps)),'x',lw=2.5,ms=10.);
    	    axSVDlapl.set_xlabel("Tune")
    	    axSVDlapl.set_ylabel(plan+" spectrum amplitude [a.u.]")
    	    axSVDlapl.set_title('Time pattern spectrum (Laplace), '+fil+', '+plane+', ik='+str(ik))

            if (save==None):
		end_figure(figSVDtune,axSVDtune);
		end_figure(figSVDlapl,axSVDlapl);
	    else:
	        end_figure(figSVDtune,axSVDtune,save=save+'_SVD_modek'+str(ik)+'_time_tune_Sussix')
		end_figure(figSVDlapl,axSVDlapl,save=save+'_SVD_modek'+str(ik)+'_time_Laplace');
	    
	    
    
    return tuneshift,tau,amp,tumost;
     

if __name__ == "__main__":


    opt,args=parsse();
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    # create list of associated legends (either with opt.LEG or 
    # the list of file names taking out all the common parameters
    # in the names)
    if (opt.LEG!=None):
        listleg=opt.LEG;
    else:
        listleg=takeout_common(listname);
    
    
    # loop on filenames
    for ifile,filename in enumerate(listname):
    
    	print filename;
    
	if not(opt.SAVE): txtsave=None;
	else: txtsave=filename[:-4];
	
	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    
	cfgfile=filename[:-8]+".cfg";
	Qx=float(read_cfg(cfgfile,"Horizontal_tune"));
	Qy=float(read_cfg(cfgfile,"Vertical_tune"));  
	Qs=float(read_cfg(cfgfile,"Synchrotron_tune"));
	Qx=Qx-floor(Qx);Qy=Qy-floor(Qy); 
	betax=float(read_cfg(cfgfile,"Horizontal_beta_function"));
	betay=float(read_cfg(cfgfile,"Vertical_beta_function"));  

	# text for the legend
        #fil=filename.replace("_"," ").replace("_prt.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read prt file
	x,xp,y,yp=read_prt_file(filename,nbunch)
	# convert x and y to 2D matrix (nturns*nbunch)
	x2d=np.reshape(x,(-1,nbunch))
	y2d=np.reshape(y,(-1,nbunch))
	
	# extract turns of interest
	x=x[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(x))];
	xp=xp[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(xp))];
	y=y[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(y))];
	yp=yp[nbunch*opt.TURNS[0]:min(nbunch*opt.TURNS[1],len(yp))];
	# complex x & y
	xcomp=x-1j*betax*xp;
	ycomp=y-1j*betay*yp

	# convert compelx x and y to 2D matrix (nturns*nbunch)
	x2dcomp=np.reshape(xcomp,(-1,nbunch))
	y2dcomp=np.reshape(ycomp,(-1,nbunch))
	
	if (opt.RAW):
	
	    # plot 1D data (along the bunch train, for the last turn)
	    fig1x,ax1x=init_figure();
	    fig1y,ax1y=init_figure();

	    plot_bunch(nbunch,x[-nbunch:],'Turn '+str(opt.TURNS[1]-1),'-b','Average x position of the bunch [m]','Raw position along the bunch train, '+fil+', horizontal',ax1x);
	    plot_bunch(nbunch,y[-nbunch:],'Turn '+str(opt.TURNS[1]-1),'-b','Average y position of the bunch [m]','Raw position along the bunch train, '+fil+', horizontal',ax1y);
	    
	    # plot 2D raw data
	    fig2x,ax2x=init_figure();
	    fig2y,ax2y=init_figure();

	    plot2D(np.transpose(x2d),opt.TURNS[0],opt.TURNS[1]-1,0,nbunch-1,'Turn number','Bunch number (0 is the head of the bunch train)','Raw positions, '+fil+', horizontal',ax2x)		    
	    plot2D(np.transpose(y2d),opt.TURNS[0],opt.TURNS[1]-1,0,nbunch-1,'Turn number','Bunch number (0 is the head of the bunch train)','Raw positions, '+fil+', vertical',ax2y)		    

            if not(opt.SAVE):
		end_figure(fig1x,ax1x);
        	end_figure(fig1y,ax1y);
		end_figure(fig2x,ax2x);
        	end_figure(fig2y,ax2y);
	    else:
		end_figure(fig1x,ax1x,save=txtsave+'_raw_x_vs_bunches_last_turn');
        	end_figure(fig1y,ax1y,save=txtsave+'_raw_y_vs_bunches_last_turn');
		end_figure(fig2x,ax2x,save=txtsave+'_raw_x_2D');
        	end_figure(fig2y,ax2y,save=txtsave+'_raw_y_2D');


	if (opt.FFT1): 

	    # 1D FFT for each bunch
	    fft1d_bunch(xcomp,nbunch,Qx,0.02,fil,'x',save=txtsave);
	    fft1d_bunch(ycomp,nbunch,Qy,0.02,fil,'y',save=txtsave);

	    # 1D FFT for each turn (along bunch train)
	    fft1d_turn(x,nbunch,opt.TURNS[0],fil,'x',save=txtsave)
	    fft1d_turn(y,nbunch,opt.TURNS[0],fil,'y',save=txtsave)
	
	
	if (opt.FFT2): 

	    # 2D FFT
	    fileout=open(filename[:-4]+'_'+opt.OUT,'w')

	    xmax=fft2d(x2dcomp,Qx,0.02,fil,'x',save=txtsave)
	    ymax=fft2d(y2dcomp,Qy,0.02,fil,'y',save=txtsave)
	    
	    print >> fileout, "2D FFT: for x & y: (mode,tune) of maximum: ", xmax, "\t", ymax;
	    print "2D FFT: for x & y: (mode,tune) of maximum: ", xmax, "\t", ymax;

	    fileout.close()
	
	
	
	if (opt.SVD):
	
	    figSVDspat,axSVDspat=init_figure();
	    
	    # SVD (time domain analysis)
	    tuneshiftx,taux,ampx,tumostx=svd(x2d,Qx,0.1,max(opt.TURNS[0],0),min(opt.TURNS[1],max(opt.TURNS[0],0)+x2d.shape[0]),fil,'x',Qs=Qs,nk=opt.NK,flagspec=opt.SPECSVD,figSVDspatial=figSVDspat,axSVDspatial=axSVDspat,nmodes=opt.NMOD,save=txtsave)
	    tuneshifty,tauy,ampy,tumosty=svd(y2d,Qy,0.1,max(opt.TURNS[0],0),min(opt.TURNS[1],max(opt.TURNS[0],0)+y2d.shape[0]),fil,'y',Qs=Qs,nk=opt.NK,flagspec=opt.SPECSVD,figSVDspatial=figSVDspat,axSVDspatial=axSVDspat,nmodes=opt.NMOD,save=txtsave)
	    
	    if (opt.SPECSVD):
		# write Laplace tune shifts and tau to files
		mmin=-1;mmax=1;
		for ik in range(opt.NK):
		    filetau=open(filename[:-4]+'_SVD_modek'+str(ik)+'_Laplace_tau.txt','w')
		    print >> filetau, "Headtail_mode\ttaux[s]\ttauy[s]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"

		    for j,m in enumerate(range(mmin,mmax+1)):
	    		print >> filetau, m, "\t", taux[ik,j]*Trev, "\t", tauy[ik,j]*Trev,"\t", tuneshiftx[ik,j], "\t", tuneshifty[ik,j], "\t", abs(ampx[ik,j]), "\t", abs(ampy[ik,j]);

	            filetau.close();

		    # highest amplitude
		    filemostamp=open(filename[:-4]+'_SVD_modek'+str(ik)+'_Laplace_most_amp_tau.txt','w')
		    print >> filemostamp, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
		    print >> filemostamp, tumostx[ik,1,0]/Trev, "\t", tumosty[ik,1,0]/Trev, "\t", tumostx[ik,0,0], "\t", tumosty[ik,0,0], "\t", abs(tumostx[ik,2,0]), "\t", abs(tumosty[ik,2,0]);
	            filemostamp.close();

		    # highest growth rate
		    filemostrate=open(filename[:-4]+'_SVD_modek'+str(ik)+'_Laplace_most_rate_tau.txt','w')
		    print >> filemostrate, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tAmp_x\tAmp_y"
		    print >> filemostrate, tumostx[ik,1,1]/Trev, "\t", tumosty[ik,1,1]/Trev, "\t", tumostx[ik,0,1], "\t", tumosty[ik,0,1], "\t", abs(tumostx[ik,2,1]), "\t", abs(tumosty[ik,2,1]);
	            filemostrate.close();
			
    if not(opt.SAVE): pylab.show();
