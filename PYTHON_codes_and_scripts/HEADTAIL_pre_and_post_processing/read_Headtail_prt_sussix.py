#!/usr/bin/python

import sys
import commands
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

from parser_lib import *
from string import *
import numpy as np
import pylab,os,re
from scipy import optimize
from SussixNM import *
from read_cfg import read_cfg
from plot_lib import plot,init_figure,end_figure,make_movie,plot2D
from io_lib import list_files
from read_Headtail_prt import extract_Headtail_param
from string_lib import takeout_common,takeout_spaces
from read_Headtail_prt_laplace import find_Laplace
from read_Headtail_prt_fit import read_prt_file

def parsse():
    parser = OptionParser()
    parser.add_option("-a", "--analysis",action="store_true",
                      help="Specify if we do an analysis of headtail modes (i.e. fit for rise time + tune shifts determination for several m from -m option)",
                      metavar="ANA", default=False,dest="ANA")
    parser.add_option("-b", "--bunch",action="callback",callback=bunch_parse,
                      help="Specify bunch number (beginning at 1) or list of bunches (i.e. 1:5:2 = [1 3 5]) for the spectra animated plot and 2D plot(several -b options possible)",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--clim",type=float,nargs=2,
                      help="Specify the limits of the color scale for the 2D plot (default: autoscale, which can be different for each plot)",
                      metavar="CLIM", default=None,dest="CLIM")
    parser.add_option("-d", "--bound",type=float,nargs=2,
                      help="Specify bounds for the vertical scale for the spectrum plot (default=[1e-8,1e-1])",
                      metavar="BOUND", default=[1.e-8,1.e-1],dest="BOUND")		      
    parser.add_option("-e", "--every",type=int,
                      help="Specify the number of bunches to skip in the analysis- i.e. we analyse only one every 'this number' bunch (default=1)",
                      metavar="EVERY", default=1,dest="EVERY")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the Headtail _prt.dat name of the file to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-g", "--beginning",type=int,
                      help="Specify the number of turns from the beginning after which we want to begin the fit (default=0)",
                      metavar="BEG", default=0,dest="BEG")
    parser.add_option("-i", "--istun",type=float,
                      help="Specify istun option (default=0.02)",
                      metavar="ISTUN", default=0.02,dest="ISTUN")
    parser.add_option("-j", "--legend",action="append",
                      help="Specify the legend for the plot, for each file",
                      metavar="LEG", default=None,dest="LEG")
    parser.add_option("-k", "--makemovie",action="store_true",
                      help="Specify if we plot sucessive spectra and make movie with them",
                      metavar="MAKE", default=False,dest="MAKE")
    parser.add_option("-l", "--laplace",action="store_true",
                      help="Specify if we also analyse with Laplace method",
                      metavar="LAPL", default=False,dest="LAPL")
    parser.add_option("-m", "--mminmmax",type=int,nargs=2,
                      help="Specify minimum and maximum headtail mode to find (with -a option)",
                      metavar="MMAX", default=[0,1],dest="MMAX")
    parser.add_option("-n", "--bunchmin",type=int,
                      help="Specify the minimum bunch number to analyse (we analyse the bunches from nbunch - last bunch - to this bunch number) (default=1)",
                      metavar="NBMIN", default=1,dest="NBMIN")
    parser.add_option("-o", "--output",help="Specify output suffix names of files with general (multifiles) plots (default=empty string)",
                      default="",dest="OUT")
    parser.add_option("-t", "--turns",type=int,nargs=2,
                      help="Specify number of turns to calculate the spectra on, and number of turns between two successive windows (default=[1024,100])",
                      metavar="TURNS", default=[1024,100],dest="TURNS")
    parser.add_option("-s", "--stopturn",type=int,
                      help="Specify the turn at which we stop analysing",
                      metavar="STOP", default=10000000,dest="STOP")
    parser.add_option("-u", "--boundmin",type=float,
                      help="Specify minimum amplitude of spectral lines to begin the fit (default=1e-7)",
                      metavar="BMIN", default=1.e-7,dest="BMIN")
    parser.add_option("-v", "--average",action="store_true",
                      help="Specify if we also plot (and analyse) the average beam position",
                      metavar="AVER", default=False,dest="AVER")
    parser.add_option("-w", "--width",type=float,
                      help="Specify the spectrum width around the unperturbed tune (default=0.02)",
                      metavar="WIDTH", default=0.02,dest="WIDTH")
    parser.add_option("--save",action="store_true",
                      help="Specify if all the plots are saved right away (.png and .eps files). Then they are not printed on the screen.",
                      metavar="SAVE", default=False,dest="SAVE")
    (opt, args) = parser.parse_args()
    return opt, args


    

def gettune(x,xp,y,yp,betax,betay,Qx,Qy,Qs,istun=0.02,narm=300):

    px=betax*xp;
    py=betay*yp;
    nturns=len(x);
    #xfft=[];yfft=[];

    #alpha=optimize.fminbound(maxfftshifted,0,1./float(len(x)),args=(x,0,len(x)),xtol=1.e-12);
    #beta=optimize.fminbound(maxfftshifted,0,1./float(len(y)),args=(y,0,len(y)),xtol=1.e-12);
    #tunex=float(argmax(np.abs(fftshifted(alpha,x))))/float(len(x))-alpha;
    #tuney=float(argmax(np.abs(fftshifted(beta,y))))/float(len(y))-beta;
    sussix_inp(ir=0, turns=nturns, tunex=Qx, tuney=Qy, tunez=Qs, istun=istun, idam=2, narm=narm)
    a=sussixBS(x,y,px,py)
    #print 'sussix tune shifts: ', a.tunexy[0]-Qx,a.tunexy[1]-Qy
    #print 'fftshifted tune shifts: ', tunex-Qx,tuney-Qy
    # those two are exactly the same
    #print a.amplitude
    #for i in range(len(a.ox)):
    #    print a.ox[i], a.ax[i],a.oy[i], a.ay[i]

    
    return a#,xfft,yfft


def plot_spec(x,data,title,col,lab,ax,fig,Q,i,bounds=None,lw=1.,width=0.01,flagleg=True,flaglog=True,flagclear=True,leg=None):

    #lab='x' or 'y'
    # plotting spectrum
    if (leg==None): lege='turn '+str(i);
    else: lege=leg;
    
    if flagclear: ax.cla()
    
    if flaglog:
	if (flagleg): ax.semilogy(x,data,col,label=lege,linewidth=lw);
	else: ax.semilogy(x,data,col,linewidth=lw);	
    else:
	if (flagleg): ax.plot(x,data,col,label=lege,linewidth=lw);
	else: ax.plot(x,data,col,linewidth=lw);	
    
    ax.set_xlabel("Tune")
    ax.set_ylabel("Spectrum amplitude [a.u.]")
    ax.set_title(title)
    ax.set_xlim([Q-width,Q+width])
    if (bounds!=None): ax.set_ylim(bounds)
    ax.legend(loc=0)
    if flagclear:
	fname = '_tmp'+lab+'%d.png'%i
	fnameeps = '_tmp'+lab+'%d.eps'%i
	fig.savefig(fname)
	fig.savefig(fnameeps,format='eps')
    

def extract(tunes,amps,Q,width):

    # sort and extract only tunes at width from the tune Q
    dic=dict(transpose([tunes,amps]))
    tunessort=np.sort(tunes)
    m=np.searchsorted(tunessort,Q-width)
    p=np.searchsorted(tunessort,Q+width,side='right')
    newtunes=tunessort[m:p]
    newamps=np.array([dic[i] for i in newtunes]);
    
    return newtunes,newamps;
    

def find_peak(tunes,amps,Q,Qs,m):

    # find the headtail mode m from the spectrum of amplitude 'amp' vs. the tunes 'tunes'
    tu,am=extract(tunes,amps,Q+m*Qs,Qs/2.)
    if (len(tu)>0):
	# tune shift of the peak from the zero intensity value Q+mQs
	tunem=tu[np.argmax(am)]-(Q+m*Qs)
	# peak amplitude
	ampm=np.max(am)
    else:
        tunem=0;ampm=0;

    return tunem,ampm
	

def read_bunch_table(filename):

    # read bunch table file
    # return intensities for non empty bunches
    # and flag that is True if intensities are different
    inte=[];flag=True;

    for il,line in enumerate(open(filename)):
        if (len(split(line))>0)and(float(split(line)[0])>0):
            inte.append(float(split(line)[0]))

    inte=np.array(inte);
    flag=(np.max(np.abs(np.diff(inte)))!=0.);
    print inte,len(inte),flag;
    
    return inte,flag;
    

def collect_and_fit(tunes,amps,turns,lab,plane,begin,flagplot,axtu,axamp,Trev,tuneaxis,ampaxis,firstturn=0,lastturn=None,flagaver=False,xlab="Number of turns",xscaleplot=1.,xoffplot=0.):

    # collect tune shifts and amplitudes, then fit amplitudes by an exponential
    # and plot (label for legend='lab', axes for tune shift plot='axtu',
    # axes for amplitudes plot='axamp')
    # tuneaxis is the y label for the tuneshift plot
    # ampaxis is the y label for the amplitude plot
    # plane is 'x' or 'y'
    # xscaleplot and xoffplot are resp. a scaling factor and an offset for the xaxis, in case of plot 
    # -> we multiply turns by xscaleplot, then add xoffplot (for the plot only)
    tunes=np.array(tunes)
    amps=np.array(amps)
    turns=np.array(turns)

    if (lastturn!=None): last=lastturn;
    else: last=turns[-1];
	
    #fit amplitude for tau
    begt=np.where(turns>=firstturn);begt=begt[0];
    if (len(begt)==0): begt=[0];print "  warning: plane ",plane,", first turn for fit too large ",firstturn
    beg=np.where(amps[begt[0]:]>=begin);beg=beg[0];
    if (len(beg)==0): beg=[0];print "  warning: plane ",plane,", boundary for fit too large ",begin	
    beg=beg+begt[0];
    endt=np.where(turns<=last);endt=endt[0];
    if (len(endt)==0)or(endt[-1]<=beg[0]): endt=[-2];
    poly=np.polyfit(turns[beg[0]:endt[-1]+1],np.log(amps[beg[0]:endt[-1]+1]),1)
    tau=1./poly[0];ct=math.exp(poly[1]);

    if flagaver:
	# compute tuneshift average
	tuneshift=np.average(tunes[beg[0]:endt[-1]+1]);
	sigma=np.sqrt(np.average(tunes[beg[0]:endt[-1]+1]*tunes[beg[0]:endt[-1]+1])-tuneshift**2);
    else:
	#find final tune shift
	tuneshift=tunes[endt[-1]];
	
    if (flagplot):
	#plot vs turns
	plot(turns*xscaleplot+xoffplot,tunes,lab,'-',tuneaxis,axtu,0,xlab=xlab)
	plot(turns*xscaleplot+xoffplot,amps,lab,'-',ampaxis,axamp,2,xlab=xlab)
	plot(turns[beg[0]:endt[-1]+1]*xscaleplot+xoffplot,ct*np.exp(turns[beg[0]:endt[-1]+1]/tau),lab+', fit, '+r'$\tau_'+plane+'$='+"%5.2e s" % (tau*Trev),'-',ampaxis,axamp,2,lw=5.,xlab=xlab)
	axtu.legend(loc=0)
	axamp.legend(loc=0)
	
    if flagaver:
    	return tuneshift,tau,sigma
    else:
    	return tuneshift,tau
    

def intercalate(xdata,ydata,value,limit,xmin=None,xmax=None):
    # Intercalate in 'ydata' the value 'value' when the difference
    # between two successive points in 'xdata' is higher than 'limit'.
    # New points with the value 'value' are put regularly between the two initial
    # (far away) points in xdata.
    # Also add two points at 'xmin' and 'xmax' of value 'value', if no points 
    # of xdata are close to them (at less than 'limit').
    
    if (xmin!=None)and(len(pylab.mlab.find(np.abs(xdata-xmin)<limit))==0):
    	xdata=np.concatenate(([xmin],xdata));
    	ydata=np.concatenate(([value],ydata));

    if (xmax!=None)and(len(pylab.mlab.find(np.abs(xdata-xmax)<limit))==0):
    	xdata=np.concatenate((xdata,[xmax]));
    	ydata=np.concatenate((ydata,[value]));
    
    xdiff=np.diff(xdata);
    indx=pylab.mlab.find(xdiff>limit);
    
    xnew=[];ynew=[];ind0=0;ind=ind0-1;
    for i,ind in enumerate(indx):
    	if (i>0): ind0=indx[i-1]+1;
	xnew=np.concatenate((xnew,xdata[ind0:ind],np.arange(xdata[ind],xdata[ind+1],limit)));
    	ynew=np.concatenate((ynew,ydata[ind0:ind+1],np.ones(len(np.arange(xdata[ind],xdata[ind+1],limit))-1)*value));

    xnew=np.concatenate((xnew,xdata[ind+1:]));
    ynew=np.concatenate((ynew,ydata[ind+1:]));
    
    return xnew,ynew;


if __name__ == "__main__":


    opt,args=parsse();
    
    mmax=opt.MMAX[1]
    mmin=opt.MMAX[0]

    bunches=np.array(opt.BNUM);
    bunches.sort();
    print "Selected bunches:", bunches
    
    # create list of filenames to analyse
    listname=list_files(opt.FILE);
    
    if (len(listname)>4):
    	axes=[0.12,0.13,0.55,0.8];legpos=(1.03,-0.1);
    else:
    	axes=None;legpos=0;
    
    # extract common parameters from the list of file names 
    # (this is for the names of global pictures)
    # and create list of associated legends (taking out all the 
    # common parameters in the names).
    # If opt.LEG is not None, replace these legends by opt.LEG.
    listleg,com=takeout_common(listname,flagcommon=True);
    com=[nm+'_' for nm in com];
    common=''.join(com);common=common[:-1];
    if (opt.LEG!=None):
    	# replace these legends by opt.LEG.
        listleg=opt.LEG;
    else:
	listleg=takeout_spaces(listleg)
    
    # find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
    gamma,circ,Trev,nbunch=extract_Headtail_param(listname[0][:-8])
    

    if (opt.ANA):
        # initialize plots and variables for tune shifts and rise times of headtail modes
    	figtau,axtau=init_figure(axes=axes);
    	figtuneshift,axtuneshift=init_figure(axes=axes);
    	taux=np.zeros((nbunch,-mmin+mmax+1));
	tauy=np.zeros((nbunch,-mmin+mmax+1));
    	tuneshiftx=np.zeros((nbunch,-mmin+mmax+1));
	tuneshifty=np.zeros((nbunch,-mmin+mmax+1));
    	sigmax=np.zeros((nbunch,-mmin+mmax+1));
	sigmay=np.zeros((nbunch,-mmin+mmax+1));
	if (opt.AVER):
    	    taubeamx=np.zeros(-mmin+mmax+1);
	    taubeamy=np.zeros(-mmin+mmax+1);
    	    tuneshiftbeamx=np.zeros(-mmin+mmax+1);
	    tuneshiftbeamy=np.zeros(-mmin+mmax+1);
    	    sigmabeamx=np.zeros(-mmin+mmax+1);
	    sigmabeamy=np.zeros(-mmin+mmax+1);
	    
    width=opt.WIDTH; # width chosen for the spectrum 

    if (opt.AVER):
        # initialize plots for tune shifts and amplitudes of most unstable mode
	figtumostx,axtumostx=init_figure(axes=axes);
	figampmostx,axampmostx=init_figure(axes=axes);	
	figtumosty,axtumosty=init_figure(axes=axes);
	figampmosty,axampmosty=init_figure(axes=axes);	


    # loop on filenames
    for ifile,filename in enumerate(listname):
    
    	print filename;
    
	# read from cfg file gamma, circumference, beta functions and tunes
	
	# find gamma, circumference, Trev and nbunch from .cfg and _inph.dat files
	gamma,circ,Trev,nbunch=extract_Headtail_param(filename[:-8])
    
	cfgfile=filename[:-8]+'.cfg';
	betax=float(read_cfg(cfgfile,"Horizontal_beta_function"));
	betay=float(read_cfg(cfgfile,"Vertical_beta_function"));  
	Qx=float(read_cfg(cfgfile,"Horizontal_tune"));
	Qy=float(read_cfg(cfgfile,"Vertical_tune"));  
	Qs=float(read_cfg(cfgfile,"Synchrotron_tune"));
	ibunchtb=int(read_cfg(cfgfile,"Switch_for_bunch_table"));
	Qx=Qx-floor(Qx);Qy=Qy-floor(Qy); 
	print 'Betax=',betax,', Betay=',betay;
	print 'Qx=',Qx,', Qy=',Qy,', Qs=',Qs;

	if (opt.AVER)and(opt.MAKE):
            figbeamx,axbeamx=init_figure();
            figbeamy,axbeamy=init_figure();

	if (opt.ANA):
	    file=[];
	    for j,m in enumerate(range(mmin,mmax+1)):
		file.append(open(filename[:-4]+'m'+str(m)+'tau.txt','w'));
		print >> file[j], "Bunch\ttaux[s]\ttauy[s]\tTune_shiftx\tTune_shifty\tSigma_tune_shiftx\tSigma_tune_shifty"
	    if (opt.AVER):
	    	filebeam=open(filename[:-4]+'_aver_tau.txt','w')
		print >> filebeam, "Headtail_mode\ttaux[s]\ttauy[s]\tTune_shiftx\tTune_shifty\tSigma_tune_shiftx\tSigma_tune_shifty"
	
	if (opt.AVER):
	    filebeammost=open(filename[:-4]+'_aver_most_tau.txt','w')
	    print >> filebeammost, "Growthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tSigma_tune_shiftx\tSigma_tune_shifty"

	if (opt.BNUM!=None):
	    filebunchmost=open(filename[:-4]+'_most_tau.txt','w')
	    print >> filebunchmost, "Bunch\tGrowthratex[s-1]\tGrowthratey[s-1]\tTune_shiftx\tTune_shifty\tSigma_tune_shiftx\tSigma_tune_shifty"

	if (opt.LAPL):
	    file2=open(filename[:-4]+'_Laplacetau.txt','w');
	    print >> file2, "Bunch\tmx\ttaux[s]\tTune_shiftx\tmy\ttauy[s]\tTune_shifty"

    	# for the legend
        #fil=filename.replace("_"," ").replace("_prt.dat","").replace("/"," ");
	fil=listleg[ifile];
	
    	# read prt file
	x,xp,y,yp=read_prt_file(filename,nbunch)
	
	stopturn=min(len(x[::nbunch]),opt.STOP);
	turns=np.arange(0,stopturn-opt.TURNS[0],opt.TURNS[1]);
	
	if (opt.BNUM!=None):
	    # initialize 2D plots (amps vs tunes and turns) ("hump buster" like)
	    # and initialize 2D data
	    fig2Dx,ax2Dx=init_figure();
	    fig2Dy,ax2Dy=init_figure();
	    ntunes=ceil(2.*width/2.e-6);
	    data2Dx=np.zeros((len(turns),ntunes),dtype=float);
	    data2Dy=np.zeros((len(turns),ntunes),dtype=float);


	if (opt.AVER):

	    # compute average first
	    for bnum in range(nbunch,0,-1):
	
		x1=x[bnum-1::nbunch];
		xp1=xp[bnum-1::nbunch];
    		y1=y[bnum-1::nbunch];
    		yp1=yp[bnum-1::nbunch];

	    	if (bnum==nbunch):
		    xave=x1;
		    yave=y1;
		    xpave=xp1;
		    ypave=yp1;
		else:
		    xave=x1+xave;
		    yave=y1+yave;	    
		    xpave=xp1+xpave;	    
		    ypave=yp1+ypave;	    


	for bnum in range(nbunch,opt.NBMIN-1,-opt.EVERY):

	    print "bunch no ",bnum
	    
	    if (bnum in bunches)and(opt.MAKE):
		# initialize spectra plots
		figx,axx=init_figure();
		figy,axy=init_figure();
	
	    if (bnum in bunches)and(opt.ANA):	
		#initialize figures for tune shift and amplitude vs turns plots
    		figampx,axampx=init_figure();
    		figampy,axampy=init_figure();
    		figtux,axtux=init_figure();
    		figtuy,axtuy=init_figure();

	    else:
	       axampx=0;axampy=0;axtux=0;axtuy=0;

	    
	    if (opt.ANA):
		# initialize tune shifts and amplitudes vs turns (headtail modes)
		tunemx=[];tunemy=[];ampmx=[];ampmy=[];
		for m in range(mmin,mmax+1):
		    tunemx.append([]);tunemy.append([])
		    ampmx.append([]);ampmy.append([])
		    		    
		    
	    if (opt.ANA) or (bnum in bunches):
	    
		x1=x[bnum-1::nbunch];
		xp1=xp[bnum-1::nbunch];
    		y1=y[bnum-1::nbunch];
    		yp1=yp[bnum-1::nbunch];

		if (bnum in bunches):
		    tunebx=[];tuneby=[];ampbx=[];ampby=[];
		
		for iturn,turn in enumerate(turns):

		    x2=x1[turn:turn+opt.TURNS[0]]
		    xp2=xp1[turn:turn+opt.TURNS[0]]
		    y2=y1[turn:turn+opt.TURNS[0]]
		    yp2=yp1[turn:turn+opt.TURNS[0]]
		    #print turn," ",len(x2)

		    suss=gettune(x2,xp2,y2,yp2,betax,betay,Qx,Qy,Qs,istun=opt.ISTUN)

		    # sort and extract only tunes at "width" from nominal tunes
		    tunesx,ampx=extract(suss.ox,suss.ax,Qx,width)
		    tunesy,ampy=extract(suss.oy,suss.ay,Qy,width)

		    if (bnum in bunches):
			# intercalate some zero values spaced by opt.WIDTH/2. if tunes & amps are too "empty"
			# (too avoid strange colored bands in 2D plots)
			tunesx2D,ampx2D=intercalate(tunesx,ampx,np.min(suss.ax),width/2.,xmin=Qx-width,xmax=Qx+width);
			tunesy2D,ampy2D=intercalate(tunesy,ampy,np.min(suss.ay),width/2.,xmin=Qy-width,xmax=Qy+width);
			# construct data for 2D plot
			data2Dx[iturn,:]=np.interp(Qx-width+2.*width*np.arange(ntunes)/float(ntunes),tunesx2D,log(ampx2D));
			data2Dy[iturn,:]=np.interp(Qy-width+2.*width*np.arange(ntunes)/float(ntunes),tunesy2D,log(ampy2D));
			
			if (opt.MAKE):
			    # Plot amplitude of the lines vs tune (spectrum)
			    # in x
			    plot_spec(tunesx,ampx,fil+', x, bunch '+str(bnum),'-',"x",axx,figx,Qx,turn,opt.BOUND);
			    # in y
			    plot_spec(tunesy,ampy,fil+', y, bunch '+str(bnum),'-',"y",axy,figy,Qy,turn,opt.BOUND);
			
		    if (opt.ANA):
		    	# find spectral lines associated with headtail modes
			for j,m in enumerate(range(mmin,mmax+1)):
			    tu,am=find_peak(tunesx,ampx,Qx,Qs,m)
			    tunemx[j].append(tu);ampmx[j].append(am);
			    tu,am=find_peak(tunesy,ampy,Qy,Qs,m)
			    tunemy[j].append(tu);ampmy[j].append(am);
		    
		    if (bnum in bunches):
		    	# find spectral line of highest amplitude
			tu,am=find_peak(tunesx,ampx,Qx,width,0)
			tunebx.append(tu);ampbx.append(am);
			tu,am=find_peak(tunesy,ampy,Qy,width,0)
			tuneby.append(tu);ampby.append(am);
			

		if (bnum in bunches):
		    # make the 2D plot ("hump buster" like: turns vs tunes with color=amplitudes)
    		    plot2D(data2Dx,Qx-width,Qx+width-2.*width/float(ntunes),turns[0],turns[-1],'Tune',
		    	'Turn number',fil+', x, bunch '+str(bnum),ax2Dx,
			colorlabel="Log(spectrum amplitude)",colorlim=opt.CLIM,fig=fig2Dx);
    		    plot2D(data2Dy,Qy-width,Qy+width-2.*width/float(ntunes),turns[0],turns[-1],'Tune',
		    	'Turn number',fil+', y, bunch '+str(bnum),ax2Dy,
			colorlabel="Log(spectrum amplitude)",colorlim=opt.CLIM,fig=fig2Dy);
		    
		    if (opt.MAKE):
			# make a movie with the spectra
			make_movie(filename[:-4]+'_b'+str(bnum)+'_x.gif','_tmpx');
			make_movie(filename[:-4]+'_b'+str(bnum)+'_y.gif','_tmpy');
			#pylab.close(figx)
			#pylab.close(figy)
		
		    # collect tune shifts and amplitudes of most unstable mode
		    tuneshiftxmost,tauxmost,sigmaxmost=collect_and_fit(tunebx,ampbx,turns,
			    '','x',opt.BMIN,False,0,0,Trev,'','',opt.BEG,flagaver=True)
		    tuneshiftymost,tauymost,sigmaymost=collect_and_fit(tuneby,ampby,turns,
			    '','y',opt.BMIN,False,0,0,Trev,'','',opt.BEG,flagaver=True)
		    print >> filebunchmost, bnum, "\t", 1./(tauxmost*Trev), "\t", 1./(tauymost*Trev),"\t", tuneshiftxmost, "\t", tuneshiftymost,"\t", sigmaxmost, "\t", sigmaymost;


		if (opt.ANA):
		    # collect tune shifts and amplitudes of headtail modes
		    for j,m in enumerate(range(mmin,mmax+1)):
			tuneshiftx[bnum-1,j],taux[bnum-1,j],sigmax[bnum-1,j]=collect_and_fit(tunemx[j],ampmx[j],turns,
				'mode m='+str(m)+', '+fil+', x, bunch '+str(bnum),'x',opt.BMIN,
				(bnum in bunches),axtux,axampx,Trev,'Tune shift w.r.t '+r'$Q_x-mQ_s$',
				'Amplitude of the mode',opt.BEG,flagaver=True)
			tuneshifty[bnum-1,j],tauy[bnum-1,j],sigmay[bnum-1,j]=collect_and_fit(tunemy[j],ampmy[j],turns,
				'mode m='+str(m)+', '+fil+', y, bunch '+str(bnum),'y',opt.BMIN,
				(bnum in bunches),axtuy,axampy,Trev,'Tune shift w.r.t '+r'$Q_y-mQ_s$',
				'Amplitude of the mode',opt.BEG,flagaver=True)
			print >> file[j], bnum, "\t", taux[bnum-1,j]*Trev, "\t", tauy[bnum-1,j]*Trev,"\t", tuneshiftx[bnum-1,j], "\t", tuneshifty[bnum-1,j],"\t", sigmax[bnum-1,j], "\t", sigmay[bnum-1,j];
							
	
	
	    if (opt.LAPL):
		# try the same thing with complex frequency spectral analysis
		x2=x1[opt.BEG::];y2=y1[opt.BEG::];
		for j,m in enumerate(range(mmin,mmax+1)):
		    tauxlapl,tunex,x3=find_Laplace(x2);
		    tauylapl,tuney,y3=find_Laplace(y2);
		    print tunex,tuney
		    mx=round((tunex-Qx)/Qs);my=round((tuney-Qy)/Qs);
		    tuneshiftxlapl=tunex-Qx-mx*Qs;
		    tuneshiftylapl=tuney-Qy-my*Qs;
		    print >> file2, bnum, "\t", mx, "\t", tauxlapl*Trev, "\t", tuneshiftxlapl, "\t", my, "\t", tauylapl*Trev,"\t", tuneshiftylapl;
		    x2=x3;y2=y3;
		
		
	if (opt.LAPL):
	    file2.close();
	    

	if (opt.BNUM!=None):
	    filebunchmost.close();
	    end_figure(fig2Dx,ax2Dx,save=opt.SAVE*(filename[:-4]+'_x_Sussix_spectrogram'));
	    end_figure(fig2Dy,ax2Dy,save=opt.SAVE*(filename[:-4]+'_y_Sussix_spectrogram'));
	    
	
	if (opt.ANA):
    	    for j,m in enumerate(range(mmin,mmax+1)): file[j].close();
	    pat=['-','--',':','.-','x','o','+']
	    for j,m in enumerate(range(mmin,mmax+1),):
		axtau.plot(range(nbunch,opt.NBMIN-1,-opt.EVERY),taux[range(nbunch-1,opt.NBMIN-2,-opt.EVERY),j]*Trev,pat[j%8],label='mode m='+str(m)+', '+fil+' '+r'$\tau_x$',lw=2.5);
		axtau.plot(range(nbunch,opt.NBMIN-1,-opt.EVERY),tauy[range(nbunch-1,opt.NBMIN-2,-opt.EVERY),j]*Trev,pat[j%8],label='mode m='+str(m)+', '+fil+' '+r'$\tau_y$',lw=2.5);
		axtuneshift.plot(range(nbunch,opt.NBMIN-1,-opt.EVERY),tuneshiftx[range(nbunch-1,opt.NBMIN-2,-opt.EVERY),j],pat[j%8],label='mode m='+str(m)+', '+fil+' '+', x',lw=2.5);
		axtuneshift.plot(range(nbunch,opt.NBMIN-1,-opt.EVERY),tuneshifty[range(nbunch-1,opt.NBMIN-2,-opt.EVERY),j],pat[j%8],label='mode m='+str(m)+', '+fil+' '+', y',lw=2.5);


	if (opt.AVER):
	    
	    print "Beam average"
	    
	    xave/=float(nbunch);
	    yave/=float(nbunch);
	    xpave/=float(nbunch);
	    ypave/=float(nbunch);
	    
	    if (opt.ANA):
		# initialize figures for the plots of tune shift and amplitude vs turns (headtail modes)
    		figampbeamx,axampbeamx=init_figure();
    		figampbeamy,axampbeamy=init_figure();
    		figtubeamx,axtubeamx=init_figure();
    		figtubeamy,axtubeamy=init_figure();
		# initialize tune shifts and amplitudes vs turns (headtail modes)
		tunemx=[];tunemy=[];ampmx=[];ampmy=[];
		for m in range(mmin,mmax+1):
		    tunemx.append([]);tunemy.append([])
		    ampmx.append([]);ampmy.append([])
		    
	    tunebx=[];tuneby=[];ampbx=[];ampby=[];
		    
	    for turn in range(0,min(len(xave),opt.STOP)-opt.TURNS[0],opt.TURNS[1]):

		xave2=xave[turn:turn+opt.TURNS[0]]
		xpave2=xpave[turn:turn+opt.TURNS[0]]
		yave2=yave[turn:turn+opt.TURNS[0]]
		ypave2=ypave[turn:turn+opt.TURNS[0]]
		#print turn," ",len(xave2)

		suss=gettune(xave2,xpave2,yave2,ypave2,betax,betay,Qx,Qy,Qs,opt.ISTUN)

		# sort and extract only tunes at "width" from nominal tunes
		tunesx,ampx=extract(suss.ox,suss.ax,Qx,width)
		tunesy,ampy=extract(suss.oy,suss.ay,Qy,width)

		# find spectral line of highest amplitude
		tu,am=find_peak(tunesx,ampx,Qx,width,0)
		tunebx.append(tu);ampbx.append(am);
		tu,am=find_peak(tunesy,ampy,Qy,width,0)
		tuneby.append(tu);ampby.append(am);

		if (opt.MAKE):
		    # Plot amplitude of the lines vs tune (spectrum)
		    # in x
		    plot_spec(tunesx,ampx,fil+", x, average over all bunches",'-',"x",axbeamx,figbeamx,Qx,turn,opt.BOUND);
		    # in y
		    plot_spec(tunesy,ampy,fil+", y, average over all bunches",'-',"y",axbeamy,figbeamy,Qy,turn,opt.BOUND);

		if (opt.ANA):
		    # find spectral lines associated with headtail modes
		    for j,m in enumerate(range(mmin,mmax+1)):
			tu,am=find_peak(tunesx,ampx,Qx,Qs,m)
			tunemx[j].append(tu);ampmx[j].append(am);
			tu,am=find_peak(tunesy,ampy,Qy,Qs,m)
			tunemy[j].append(tu);ampmy[j].append(am);

	    if (opt.MAKE):
		make_movie(filename[:-4]+'_aver_x.gif','_tmpx');
		make_movie(filename[:-4]+'_aver_y.gif','_tmpy');


	    # collect tune shifts and amplitudes of most unstable mode
	    turns=np.array(range(0,min(len(x1),opt.STOP)-opt.TURNS[0],opt.TURNS[1]))
	    tuneshiftxmost,tauxmost,sigmaxmost=collect_and_fit(tunebx,ampbx,turns,
		    fil+', aver. x','x',opt.BMIN,True,axtumostx,axampmostx,Trev,'Tune shift w.r.t '+r'$Q_x$',
		    'Amplitude of the mode',opt.BEG,flagaver=True)
	    tuneshiftymost,tauymost,sigmaymost=collect_and_fit(tuneby,ampby,turns,
		    fil+', aver. y','y',opt.BMIN,True,axtumosty,axampmosty,Trev,'Tune shift w.r.t '+r'$Q_y$',
		    'Amplitude of the mode',opt.BEG,flagaver=True)
	    print >> filebeammost, 1./(tauxmost*Trev), "\t", 1./(tauymost*Trev),"\t", tuneshiftxmost, "\t", tuneshiftymost,"\t", sigmaxmost, "\t", sigmaymost;
	    filebeammost.close()


	    if (opt.ANA):
		# collect tune shifts and amplitudes of headtail modes
		for j,m in enumerate(range(mmin,mmax+1)):
		    tuneshiftbeamx[j],taubeamx[j],sigmabeamx[j]=collect_and_fit(tunemx[j],ampmx[j],turns,
		    	'Average, mode m='+str(m)+', '+fil+', x','x',opt.BMIN,True,axtubeamx,axampbeamx,
			Trev,'Tune shift w.r.t '+r'$Q_x-mQ_s$',
			'Amplitude of the mode',opt.BEG,flagaver=True)
		    tuneshiftbeamy[j],taubeamy[j],sigmabeamy[j]=collect_and_fit(tunemy[j],ampmy[j],turns,
		    	'Average, mode m='+str(m)+', '+fil+', y','y',opt.BMIN,True,axtubeamy,axampbeamy,
			Trev,'Tune shift w.r.t '+r'$Q_y-mQ_s$',
			'Amplitude of the mode',opt.BEG,flagaver=True)
		    print fil, ", mode m=",m," - Average rise times of the beam (x & y) in sec.: ", taubeamx[j]*Trev, "  ", taubeamy[j]*Trev;
		    print fil, ", mode m=",m," - Average tune shifts of the beam (x & y): ", tuneshiftbeamx[j], "  ", tuneshiftbeamy[j];
		    print >> filebeam, m, "\t", taubeamx[j]*Trev, "\t", taubeamy[j]*Trev,"\t", tuneshiftbeamx[j], "\t", tuneshiftbeamy[j],"\t", sigmabeamx[j], "\t", sigmabeamy[j];
   	        
		filebeam.close()



	if (opt.ANA)and(ibunchtb==1):
	    inte,flagint=read_bunch_table(filename[:-8]+'.bunch');
	    if (flagint):
		# plot also tuneshifts of headtail modes vs intensities, and fit
		for j,m in enumerate(range(mmin,mmax+1)):
    		    figtuintx,axtuintx=init_figure(axes=axes);
    		    figtuinty,axtuinty=init_figure(axes=axes);
		    
		    px=np.polyfit(inte,tuneshiftx[:,j],1);
		    py=np.polyfit(inte,tuneshifty[:,j],1);

		    axtuintx.errorbar(inte,tuneshiftx[:,j],yerr=sigmax[:,j],fmt='xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);
		    axtuintx.plot(inte,px[1]+px[0]*inte,'-b',label='Fit by a-x*%.6f'%(-px[0])+" $ /10^{11} $ ",lw=3.);
		    axtuintx.set_xlabel("Intensity $ (10^{11} $ p+/b)");
		    axtuintx.set_ylabel("Horizontal tune shift");
		    axtuintx.set_title('mode m='+str(m)+', '+fil+', x');

		    axtuinty.errorbar(inte,tuneshifty[:,j],yerr=sigmay[:,j],fmt='xb',label='Headtail data',lw=3.,ms=10.,mew=2.5);
		    axtuinty.plot(inte,py[1]+py[0]*inte,'-b',label='Fit by a-x*%.6f'%(-py[0])+" $ /10^{11} $ ",lw=3.);
		    axtuinty.set_xlabel("Intensity $ (10^{11} $ p+/b)");
		    axtuinty.set_ylabel("Vertical tune shift");
		    axtuinty.set_title('mode m='+str(m)+', '+fil+', y');
		    
		    end_figure(figtuintx,axtuintx,save=opt.SAVE*(common+'_x_Sussix_most_tau_tuneshift_vs_bunch_int'),legpos=legpos);
		    end_figure(figtuinty,axtuinty,save=opt.SAVE*(common+'_y_Sussix_most_tau_tuneshift_vs_bunch_int'),legpos=legpos);



    if (opt.ANA):
    	# finalize figures (headtail modes)
	axtau.set_xlabel("Bunch number");
	axtau.set_ylabel("Rise time [s]");
	end_figure(figtau,axtau,save=opt.SAVE*(common+'_Sussix_most_tau_tuneshift_vs_bunch'+opt.OUT),legpos=legpos);
	axtuneshift.set_xlabel("Bunch number");
	axtuneshift.set_ylabel("Tune shift");
	end_figure(figtuneshift,axtuneshift,save=opt.SAVE*(common+'_Sussix_most_tau_tuneshift_vs_bunch'+opt.OUT),legpos=legpos);
	
    if (opt.AVER):

	end_figure(figtumostx,axtumostx,save=opt.SAVE*(common+'_x_Sussix_aver_most_tau_tuneshift_vs_turns'+opt.OUT),legpos=legpos,legfontsize=max(30-5*len(listname),10));
	end_figure(figampmostx,axampmostx,save=opt.SAVE*(common+'_x_Sussix_aver_most_tau_amp_vs_turns'+opt.OUT),legpos=legpos,legfontsize=max(30-5*len(listname),10));
	end_figure(figtumosty,axtumosty,save=opt.SAVE*(common+'_y_Sussix_aver_most_tau_tuneshift_vs_turns'+opt.OUT),legpos=legpos,legfontsize=max(30-5*len(listname),10));
	end_figure(figampmosty,axampmosty,save=opt.SAVE*(common+'_y_Sussix_aver_most_tau_amp_vs_turns'+opt.OUT),legpos=legpos,legfontsize=max(30-5*len(listname),10));

    #pylab.close(figbeamx)
    #pylab.close(figbeamy)

    if not opt.SAVE: pylab.show();

    sys.exit();

