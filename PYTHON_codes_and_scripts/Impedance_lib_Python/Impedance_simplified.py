#!/usr/bin/python2.6

import sys
#from optparse import OptionParser
from string import *
import numpy as np
import pylab,re,os
from plot_lib import *
from io_lib import *
from tables_lib import *
from string_lib import invert_selection
from copy import deepcopy

# library for impedance models

class freq_param(object):
    '''class for definition of frequency sampling (all in Hz)'''
    	
    def __init__(self,fmin=10,fmax=1e13,ftypescan=2,fsamplin=1e8,nflog=20,
    	fminrefine=1e11,fmaxrefine=5e12,nrefine=5000,fadded=[1e-2,0.1,1,1e15]):
    
        self.fmin=fmin;
	self.fmax=fmax;
	self.ftypescan=ftypescan;
	self.fsamplin=fsamplin;
	self.nflog=nflog;
	self.fminrefine=fminrefine;
	self.fmaxrefine=fmaxrefine;
	self.nrefine=nrefine;
	self.fadded=fadded;


class z_param(object):
    '''class for definition of z (distance) sampling (all in m)'''
    	
    def __init__(self,zmin=0.01,zmax=1e6,ztypescan=2,nzlog=100,
        zminrefine=0.5e-5,zmaxrefine=0.01,zsamplin=0.5e-5,zadded=[]):
    
        self.zmin=zmin;
	self.zmax=zmax;
	self.ztypescan=ztypescan;
	self.nzlog=nzlog;
	self.zminrefine=zminrefine;
	self.zmaxrefine=zmaxrefine;
	self.zsamplin=zsamplin;
	self.zadded=zadded;


def resonator_impedance(R,fr,Q,freq,save=None):

    '''resonator model impedance (transverse)
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor
    computes impedance at the frequencies freq (2 columns array)
    if save different from None and contains a string, write imp. to this file name (3 columns)'''
    
    Z=np.zeros((len(freq),2));
    Zt=(R*fr/freq)/(1.-1j*Q*(fr/freq-freq/fr));
    Z[:,0]=np.real(Zt);Z[:,1]=np.imag(Zt);
    
    if (save!=None):
    	data=np.hstack((freq.reshape((-1,1)),Z));
    	write_ncol_file(save,data,header="Frequency[Hz]\tRe(Z)[Ohm/m]\tIm(Z)[Ohm/m]");
    
    return Z;
 

def resonator_wake(R,fr,Q,tau,save=None):

    '''resonator model wake function (transverse)
    R=shunt impedance (Ohm/m), fr=resonant frequency, Q=quality factor
    computes wake at the times tau (in ns) (must be different from 0) (tau>0 behind the source)
    NOTE: Q must be different from 0.5 !
    if save different from None and contains a string, write wake to this file name (2 columns)'''
    
    from scipy import sqrt,sin
    
    W=np.zeros((len(tau),1));tau=tau.reshape((-1,1))*1e-9;
    omegar=2*np.pi*fr;omegarbar=omegar*sqrt(1.-1./((2.*Q)**2));alpha=omegar/(2.*Q);
    ind=np.where(tau>0);
    W[ind]=(R*omegar**2/(Q*omegarbar))*np.exp(-alpha*tau[ind])*sin(omegarbar*tau[ind]);
    
    if (save!=None):
    	data=np.hstack((tau*1e9,1.e-15*W));
    	write_ncol_file(save,data,header="Tau[ns]\tWake[V/pC/mm]");
    
    return W;
    
    
class impedance_wake(object):
    '''class for impedance or wake fucntion model (impedances components vs. frequency,
    or wake function vs. distance z)'''
     
    def __init__(self,a=1,b=0,c=0,d=0,plane='x',var=10**np.arange(1,13.1,0.1),func=np.zeros(121)):

   	'''all in SI units here.
	
	a,b,c,d are the powers the source and test coordinates in the force (or kick) due to impedance/wake
	ex:
	 - a=b=c=d=0, plane ='z' : longitudinal constant impedance/wake
	 - a=b=c=d=0, plane ='x' : horizontal constant impedance/wake
	 - a=b=c=d=0, plane ='y' : vertical constant impedance/wake
	 - a=1, b=c=d=0, plane ='x' : horizontal dipolar impedance/wake
	 - b=1, a=c=d=0, plane ='y' : vertical dipolar impedance/wake
	 - c=1, a=b=d=0, plane ='x' : horizontal quadrupolar impedance/wake
	 - d=1, a=b=c=0, plane ='y' : vertical quadrupolar impedance/wake
	 - b=1, a=c=d=0, plane ='x' : horizontal dipolar coupled-term impedance/wake
	 - d=1, a=b=c=0, plane ='x' : horizontal quadrupolar coupled-term impedance/wake
	 - a=1, b=c=d=0, plane ='y' : vertical dipolar coupled-term impedance/wake
	 - c=1, a=b=d=0, plane ='y' : vertical quadrupolar coupled-term impedance/wake
	etc...'''

	self.a=a; # power of x_s (source horizontal position)
	self.b=b; # power of y_s (source vertical position)
	self.c=c; # power of x_t (test horizontal position)
	self.d=d; # power of y_t (test vertical position)
	self.plane=plane; # plane of the impedance/wake ('x' for horizontal impedance/wake, 'y' for vertical, 'z' for longitudinal)
	self.var=deepcopy(var); # table with varying variable (frequencies [Hz] for impedance, z [m] for wake function)
	self.func=deepcopy(func); # 2-columns table with the corresponding impedance or wake function - real and 
	# imag. parts (Ohm/m^(a+b+c+d) for impedances, V/C/m^(a+b+c+d) for wakes).
	# Usually imaginary part of wake is zero, but could be non zero (for e.g. feedback wakes).


def freqscan_from_fpar(fpar):
    '''gives frequency scan associated with a freq_param object'''
    
    if (np.mod(fpar.ftypescan,2)==0):
        delta=1./float(fpar.nflog);
    	freq=10**(np.arange(np.log10(fpar.fmin),np.log10(fpar.fmax)+delta,delta));
	if fpar.ftypescan==2:
	    delta=(fpar.fmaxrefine-fpar.fminrefine)/fpar.nrefine;
	    freq2=np.arange(fpar.fminrefine,fpar.fmaxrefine+delta,delta);
	    freq=np.concatenate((freq,freq2));
    else: # fpar.ftypescan==1
    	freq=np.arange(fpar.fmin,fpar.fmax+fpar.fsamplin,fpar.fsamplin);
    
    freq=np.concatenate((freq,fpar.fadded));
    
    return sort_and_delete_duplicates(freq);
    
    
def zscan_from_zpar(zpar):
    '''gives z (distance) scan associated with a z_param object'''
    
    if (np.mod(zpar.ztypescan,2)==0):
        delta=1./float(zpar.nzlog);
    	z=10**(np.arange(np.log10(zpar.zmin),np.log10(zpar.zmax)+delta,delta));
	if zpar.ztypescan==2:
	    z2=np.arange(zpar.zminrefine,zpar.zmaxrefine+zpar.zsamplin,zpar.zsamplin);
	    z=np.concatenate((z,z2));
    else: # zpar.ztypescan==1
	z=np.arange(zpar.zminrefine,zpar.zmaxrefine+zpar.zsamplin,zpar.zsamplin);

    z=np.concatenate((z,zpar.zadded));
    
    return sort_and_delete_duplicates(z);

    
def test_impedance_wake_comp(iw,a,b,c,d,plane):
    '''test if the impedance or wake in iw has these a,b,c,d and plane'''
    flag=((iw.a==a)and(iw.b==b))and((iw.c==c)and(iw.d==d));
    
    return (flag)and(iw.plane==plane);


def add_impedance_wake(iw_model,iw_added,weightx,weighty):

    '''add the model "iw_added" to an impedance or wake model in "iw_model".
    both iw_added and iw_model are impedance or wake models, i.e. a list of objects from the class
    impedance_wake. The two lists are not necessarily of same length, and the frequency/z scans
    are not necessarily all the same.
    weightx and weighty are resp. betax/avbetax and betay/avbetay (i.e. without the square root)
    at the place of the added impedance "iw_added".
    
    It does this "in-place" (on iw_model itself)'''
    
    # loop over all impedance objects in iw_added
    for iw in iw_added:   
     	
	# compute the impedance weight
	powx=iw.a/2.+iw.c/2.;
	powy=iw.b/2.+iw.d/2.;
    	if iw.plane.startswith('x'): powx +=1./2.;
    	elif iw.plane.startswith('y'): powy +=1./2.;
	
	weight=(weightx**powx) * (weighty**powy);
	
	# determine if in iw_model_new there is already the same kind of component as in iw.
	flag=False;
	for kiw_mod,iw_mod in enumerate(iw_model):
	    if test_impedance_wake_comp(iw_mod,iw.a,iw.b,iw.c,iw.d,iw.plane):
	    	flag=True;ksamecomp=kiw_mod;

	if flag:
	    # same kind of component as in imp already present in iw_model_new
	    iw_mod=iw_model[ksamecomp];
	    # concatenate the 2 var scans
	    vartot=sort_and_delete_duplicates(np.concatenate((iw.var,iw_mod.var)));
	    if (vartot[-1]>iw.var[-1])and(np.abs((iw.func[-1,0]+1j*iw.func[-1,1])/np.average(iw.func[:,0]+1j*iw.func[:,1]))>1e-2):
	        print " Warning in add_impedance_wake: added model extrapolated above",iw.var[-1],"by zeros; func[-1]=",iw.func[-1,0],iw.func[-1,1];
	    if (vartot[-1]>iw_mod.var[-1])and(np.abs((iw_mod.func[-1,0]+1j*iw_mod.func[-1,1])/np.average(iw_mod.func[:,0]+1j*iw_mod.func[:,1]))>1e-2):
	        print " Warning in add_impedance_wake: initial model extrapolated above",iw_mod.var[-1],"by zeros; func[-1]=",iw_mod.func[-1,0],iw_mod.func[-1,1];
	    # sum impedances
	    functot=np.zeros((len(vartot),2),dtype=float);
	    for icol in range(2):
	    	functot[:,icol]=np.interp(vartot,iw_mod.var,iw_mod.func[:,icol],right=0.)+weight*np.interp(vartot,iw.var,iw.func[:,icol],right=0.);
	    # include it in the new model
	    iw_mod.var=vartot;
	    iw_mod.func=functot;
	    
	else:
	    # same kind of component as in iw is not already present: 
	    # append to iw_model_new a new impedance/wake with this component (with the correct weight)
	    iw2=deepcopy(iw);iw2.func=weight*iw.func;
	    iw_model.append(iw2);
	    
    return iw_model;   
    

def multiply_impedance_wake(iw_model,factor):

    '''multiply all components of the model "iw_model" by the factor "factor"
    It does this "in-place" (on iw_model itself)'''
    
    # loop over all impedance objects in iw_model
    for iw in iw_model:   
     	
	iw.func *= factor;
	
    return iw_model;   
    

def imp_model_transverse_resonator(R,fr,Q,beta=1,wake_calc=False,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),plane=''):

    '''define an imp/wake model for a transverse resonator (SI units)
    R, fr, Q: shunt impedance, resonance frequency and quality factor
    beta: relativistic velocity factor (for z to tau conversion
    if wake_calc==True, calculate wake also
    assume axisymmetry (purely dipolar)
    plane can be 'x' (only to dipolar x component), 'y' (only to dipolar y component), 
    or '' (both x and y dipolar components)'''

    imp_mod=[];wake_mod=[];
    
    freq=freqscan_from_fpar(fpar);
    Zt=resonator_impedance(R,fr,Q,freq);
    if (plane=='')or(plane=='x'): imp_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=freq,func=Zt));
    if (plane=='')or(plane=='y'): imp_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=freq,func=Zt));
    
    if wake_calc:
        z=zscan_from_zpar(zpar);
	tau=z*1e9/(beta*299792458); # conversion in ns
	Wt=np.zeros((len(z),2));
    	tmp=resonator_wake(R,fr,Q,tau);
	Wt[:,0]=tmp[:,0];
	if (plane=='')or(plane=='x'): wake_mod.append(impedance_wake(a=1,b=0,c=0,d=0,plane='x',var=z,func=Wt));
	if (plane=='')or(plane=='y'): wake_mod.append(impedance_wake(a=0,b=1,c=0,d=0,plane='y',var=z,func=Wt));

    return imp_mod,wake_mod;
    

def imp_model_from_HOMfile(filename,beta=1,fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    '''Create an impedance or wake model from a file containing HOMs (high
    order modes) data, i.e. resonator models data. Impedance / wake is in 
    the end of sum of all those modes.
    In the file each line is one mode, and each column should have a header indicating:
    - Rl: longitudinal shunt impedance
    - Rxd(ip): horizontal dipolar impedance shunt impedance
    - same with Ryd(ip), Rxq(uad), Ryq(uad)
    - same with Ql, Qxd(ip), etc. for each quality factor
    - same with fl, fxd(ip), etc. for each resonance frequency.
    beta indicates the relativistic velocity factor.
    '''
    
    colname=['Rl','Rxd','Ryd','Rxq','Ryq'];
    Rlist=[];frlist=[];Qlist=[];complist=[];
    
    # read file
    for icol,col in enumerate(colname):
    
	R=read_ncol_file_identify_header(filename,col,dispflag=False);
	fr=read_ncol_file_identify_header(filename,col.replace('R','f'),dispflag=False);
	Q=read_ncol_file_identify_header(filename,col.replace('R','Q'),dispflag=False);
	
	if (len(R)*len(fr)*len(Q)>0):
	    complist.append(col.replace('R','Z'));
	    Rlist.append(R);frlist.append(fr);Qlist.append(Q);
    
    if (len(complist)==0): print "Pb in imp_model_from_HOMfile: no valid data in HOM file!";sys.exit();
    
    # construct model
    imp_mod=[];wake_mod=[];
    
    # sum all modes
    for i in range(len(Rlist[0])):
    	# compute resonator model for each mode
	Rlistcomp=[Rlist[icomp][i] for icomp in range(len(complist))];
	frlistcomp=[frlist[icomp][i] for icomp in range(len(complist))];
	Qlistcomp=[Qlist[icomp][i] for icomp in range(len(complist))];
	
	imp,wake=imp_model_resonator(Rlistcomp,frlistcomp,Qlistcomp,beta=beta,
		wake_calc=True,fpar=fpar,zpar=zpar,listcomp=complist);
	
	add_impedance_wake(imp_mod,imp,1,1);
    	add_impedance_wake(wake_mod,wake,1,1);
	
	
    return imp_mod,wake_mod;
    

def imp_model_from_file(filename,compname,ignored_rows=0,sign=1):
    '''Create an impedance or wake model from a file containing an impedance
    or wake function (typically from a simulation).
    File has 'ignored_rows' header lines.
    There should be one column with frequency [Hz] or distance [m] (first column)
    then one or two columns with real and imaginary part of impedance or wake,
    for the component given in 'compname' (component identification should respect 
    the standard given in function "identify_component").
    Type (wake or impedance) is detected from first letter (W or Z) of 'compname'
    (but actually it plays no role here).
    sign is the sign convention: 1 if keep the same sign as in the file,
    -1 otherwise. (use -1 for GdFidl wakes for instance).'''
    
    a,b,c,d,plane,wakeflag=identify_component(compname);
    imp_mod=[];
    
    s=read_ncol_file(filename,ignored_rows=ignored_rows);
    
    func=np.zeros((len(s[:,0]),2));
    if (len(s[0,:])>=3):
    	func=sign*s[:,1:3];
    elif (len(s[0,:])==2):
    	func[:,0]=sign*s[:,1];
    else:
    	print "imp_model_from_file: not enough columns in file "+filename;
	sys.exit();
    
    imp_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=s[:,0],func=func));
    
    return imp_mod;
    

def imp_model_from_files(filenamelist,scan,value,compname,ignored_rows=0,sign=1):
    '''Create an impedance or wake model from several files, each containing
    an impedance or wake function (typically from simulations).
    Each file has 'ignored_rows' header lines, and should have one column 
    with frequency [Hz] or distance [m] (first column)
    then one or two columns with real and imaginary part of impedance or wake,
    for the component given in 'compname' (component identification should respect 
    the standard given in function "identify_component").
    The different files represent a scan given by 'scan' (e.g. a half-gap scan)
    and we interpolate between them at the value 'value' (extrapolation is 
    possible: we use the closest value in scan).
    Type (wake or impedance) is detected from first letter (W or Z) of 'compname'
    (but actually it plays no role here).
    sign is the sign convention: 1 if keep the same sign as in the files,
    -1 otherwise. (use -1 for GdFidl wakes for instance).'''
    
    from scipy import interpolate as itp
    
    a,b,c,d,plane,wakeflag=identify_component(compname);
    
    # build the total scan interpolation table
    for ifile,filename in enumerate(filenamelist):
    
	# read file
	s=read_ncol_file(filename,ignored_rows=ignored_rows);
	
	# choose first file to define the distance or frequency scan
	if ifile==0: var=s[:,0];inttable=np.zeros((len(filenamelist),len(var),2));

	func=np.zeros((len(s[:,0]),2));
	if (len(s[0,:])>=3):
	    func=sign*s[:,1:3];
	elif (len(s[0,:])==2):
	    func[:,0]=sign*s[:,1];
	else:
	    print "imp_model_from_file: not enough columns in file "+filename;
	    sys.exit();
	    
	# construct final interpolation table
	inttable[ifile,:,0]=np.interp(var,s[:,0],func[:,0]);
	inttable[ifile,:,1]=np.interp(var,s[:,0],func[:,1]);
    
    # interpolated function
    fint=itp.interp1d(scan,inttable,axis=0);
    
    # first treat extrapolation cases
    if (value<scan[0]): funcfinal=inttable[0,:,:].squeeze();
    elif (value>scan[-1]): funcfinal=inttable[-1,:,:].squeeze();
    else: # general case
    	funcfinal=fint(value);
    
    imp_mod=[];
    imp_mod.append(impedance_wake(a=a,b=b,c=c,d=d,plane=plane,var=var,func=funcfinal));
    
    return imp_mod;
    

def write_wake_HEADTAIL(wake_mod,filename,beta=1,ncomp=6):

    '''write a wake file, in the HEADTAIL format, from the wake model 'wake_mod'
    reminder: units are ns and kV/pC(/m)
    beta: relativistic beta factor (for unit conversion - m to ns)
    ncomp= number of components to take (correspond to flag i_waketb in HEADTAIL input file)
    NOTE: HEADTAIL assumes xy=yx for coupling terms; we also make this assumption here'''
    
    c=299792458; # speed of light
    
    listcomp=['Wlong','Wxdip','Wydip','Wxquad','WZyquad','Wxydip','Wxyquad','Wxcst','Wycst']
    lista=[0,1,0,0,0,0,0,0,0];listb=[0,0,1,0,0,1,0,0,0];
    listc=[0,0,0,1,0,0,0,0,0];listd=[0,0,0,0,1,0,1,0,0];
    listplane=['z','x','y','x','y','x','x','x','y'];
    
    # components to select (depends on ncomp=i_waketb; if its odd Wlong is in the list)
    ind=np.arange(1-np.mod(ncomp,2),1-np.mod(ncomp,2)+ncomp);
    
    k_comp=[];
    for icomp in ind:
        flagfound=False;
	for kw,w in enumerate(wake_mod):
	    if test_impedance_wake_comp(w,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): k_comp.append(kw);flagfound=True;
	if (not flagfound): print "Pb in write_wake_HEADTAIL: component "+listcomp[icomp]+" not found !";sys.exit();
	
    # collect all distances z into one array
    z=[];
    for k in k_comp: z=sort_and_delete_duplicates(np.concatenate((z,wake_mod[k].var)));
    
    # unit conversion and build the final table
    data=z.reshape((-1,1))*1e9/(beta*c);
    for k in k_comp:
    	s=1e-15*np.interp(z,wake_mod[k].var,wake_mod[k].func[:,0]);
    	data=np.hstack((data,s.reshape((-1,1))));
    
    # write in file
    write_ncol_file(filename,data);
    
    return;


def write_imp_wake_mod(imp_mod,name,listcomp=['Zxdip','Zydip'],dire=''):

    '''write files with components of an impedance or wake model given in imp_mod.
    filenames will be: [component name][name].dat 
    listcomp= components to be written (each in a separate file, with 3 columns:
    frequency - or distance, real and imaginary parts)
    dire is the directory where to put the files (do not forget the "/")'''
    
    units=["","/m","/m^2","/m^3","/m^4"];
    
    for icomp,comp in enumerate(listcomp):
    
    	a,b,c,d,plane,wakeflag=identify_component(comp);
	unit=units[a+b+c+d];
    
	if wakeflag: header="Distance[m] Re_"+comp+"[V/C"+unit+"] Im_"+comp+"[V/C"+unit+"]";
	else: header="Frequency[Hz] Re_"+comp+"[Ohm"+unit+"] Im_"+comp+"[Ohm"+unit+"]";
	
	for iw in imp_mod:
	    flag=False;
	    if test_impedance_wake_comp(iw,a,b,c,d,plane): Z=deepcopy(iw.func);freq=deepcopy(iw.var);flag=True;
	    if flag:
	    	data=np.hstack((freq.reshape((-1,1)),Z.reshape((-1,2))));
	    	write_ncol_file(dire+comp+name+'.dat',data,header=header);
    
    return;


def identify_component(compname):
    '''identify a component from its name, and convert it to the impedance/wake class
    format in terms  of a, b, c, d, plane and wakeflag (wakeflag=True if this is a 
    wake component).
    first letter of compname should be 'Z' (for impedance) or 'W' (for wake potential)
    for e.g. impedance, compname can be: Zl(ong), Zxd(ip), Zyd(ip), Zxq(uad), Zyq(uad), Zxyd(ip), Zyxd(ip),
    Zxyq(uad), Zyxq(uad), Zxc(st), Zyc(st), or a name directly indicating a, b, c, d 
    and the plane, in the form e.g. "Z1000x" (Zxdip in this case)
    (for wakes, it is the same replacing 'Z' by 'W')'''
    
    if (not compname.startswith('Z'))and(not compname.startswith('W')):
    	print "Pb in identify_component: component name does not begin by Z or W";
	sys.exit();
    
    else:
    	first=compname[0]; # Z or W
	if first=='W': wakeflag=True;
	else: wakeflag=False;
    	
	if compname.startswith(first+'l'): a=0;b=0;c=0;d=0;plane='z';
	elif compname.startswith(first+'xd'): a=1;b=0;c=0;d=0;plane='x';
	elif compname.startswith(first+'yd'): a=0;b=1;c=0;d=0;plane='y';
	elif compname.startswith(first+'xq'): a=0;b=0;c=1;d=0;plane='x';
	elif compname.startswith(first+'yq'): a=0;b=0;c=0;d=1;plane='y';
	elif compname.startswith(first+'xyd'): a=0;b=1;c=0;d=0;plane='x';
	elif compname.startswith(first+'yxd'): a=1;b=0;c=0;d=0;plane='y';
	elif compname.startswith(first+'xyq'): a=0;b=0;c=0;d=1;plane='x';
	elif compname.startswith(first+'yxq'): a=0;b=0;c=1;d=0;plane='y';
	elif compname.startswith(first+'xc'): a=0;b=0;c=0;d=0;plane='x';
	elif compname.startswith(first+'yc'): a=0;b=0;c=0;d=0;plane='y';
	
	else: 
    	    try:
		a=int(compname[1]);b=int(compname[2]);
		c=int(compname[3]);d=int(compname[4]);
		plane=compname[5];
	    except IndexError:
		print "Pb in identify_component: bad length for component name";
		sys.exit();
	    except ValueError:
		print "Pb in identify_component: bad format for component name";
		sys.exit();
	    
    return a,b,c,d,plane,wakeflag;
    
    
def plot_compare_imp_model(imp_mod_list,leglist,listcomp=['Zlong','Zxdip','Zydip'],
	figimp=None,aximp=None,figratio=None,axratio=None,saveimp='',saveratio='',
	xlim=[1e3,1e11],ylim=[1e5,1e9],yliml=[1,1e6],bounds=[20e6,2e9],beta=1.,legpos=0,
	plotpercent=False,legpercentpos=(0.8,0.6),maxpercent=110.):
    
    '''plot on the same plot various impedance models provided in imp_mod_list
    (list of impedance models, each of them being itself a list of components).
    comparison is done for a set of components given in listcomp (see format in 
    function 'identify_component' just above).
    The impedance ratio is also plotted (w.r.t. to the first model in the list).
    - if figimp /aximp (resp. figratio / axratio) are provided and not None,
    imp. plots (resp. imp. ratio plots) are made on the corresponding axes 
    (if they are lists of the length of listcomp, then each component is plotted 
    on a separate plot).
    otherwise figures and axes are created here and all components are plotted
    separately.
    - leglist is the list of legends (same length as imp_mod_list)
    - if saveimp (resp. saveratio) is a string of non-zero length, save plot 
    on this filename.
    - xlim, ylim and yliml are the axes limits (only for impedance plot for ylim
    and only for longitudinal impedance for yliml)
    - bounds are the frequencies between which we compute the maximum ratio;
    they are converted to distances in the case of a wake (beta is used in 
    that case, for the conversion)
    - If plotpercent is True, we plot (on a single plot per component) the percentages
    w.r.t to the first model given, in a "filled" manner (see N. Mounet PhD),
    instead of curves
    - legpercentpos: position (bbox_to_anchor) of legend for percentage plot.
    - maxpercent: maximum of y-axis for the percentage plot (in percent).
    
    It also works for wakes, with component names beginning by "W" instead of "Z".
    Then xlim and ylim have to be changed accordingly (e.g. resp. [1e-5,1e6] and [1e8,1e19] for LHC total wakes)'''
    
    from matplotlib.patches import Rectangle
    from matplotlib.font_manager import FontProperties
    
    clight=299792458;
    
    # if figure and axes are not None, check if they are lists, if not replace them 
    # by a list with the same length as listcomp (repeating the same element)
    for name in ['imp','ratio']:
	if (eval('fig'+name+'!=None'))or(eval('ax'+name+'!=None')):
	    exec('fig'+name+'=create_list_for_figure(fig'+name+',n=len(listcomp))');
	    exec('ax'+name+'=create_list_for_figure(ax'+name+',n=len(listcomp))');
	else:
	    for b in ['fig','ax']: exec(b+name+'=[]');
	    for i in range(len(listcomp)):
	    	for ir in range(plotpercent+1): # do twice more figures (one for real, one for imag.) if "filled percentage plot"
		    if (len(imp_mod_list)<=6): fig,ax=init_figure();
		    else: fig,ax=init_figure(axes=[0.08,0.12,0.55,0.8],figsize=(20,12));legpos=(1.05,0.3);legpercentpos=(1.65,0.8);
	    	    for b in ['fig','ax']:eval(b+name+'.append('+b+')');
    
    pat=['-','--'];units=["","/m","/m$^2$","/m$^3$","/m$^4$"];
    #col=['b','r','g','m','k','c','y']; # colors
    col=build_colors(len(imp_mod_list),randomize=True);
    maxratio=np.zeros((len(listcomp),len(imp_mod_list)-1,2));
    partcomp=['real','imag'];
    
    for icomp,comp in enumerate(listcomp):
    	# identify the kind of component in terms of a,b,c,d and plane
	a,b,c,d,plane,wakeflag=identify_component(comp);
	unit=units[a+b+c+d];
	
	if wakeflag:
	    parts=[''];str1="wakes";str2="Wake";
	    xlab="Distance behind the source [m]";ylab="Wake [V/C";
	    boundmin=beta*clight/bounds[1];boundmax=beta*clight/bounds[0];
	else:
	    parts=['Real part, ','Imag part, '];str1="impedances";str2="Imp.";
	    xlab="Frequency [Hz]";ylab="Z [$\Omega$";
	    boundmin=bounds[0];boundmax=bounds[1];
	    	
	if plotpercent:
	    p=[];leglistless=[];
	    for ir,r in enumerate(parts): p.append([]);leglistless.append([]);
	
	kiw0_comp=-1;
	for imod,imp_mod in enumerate(imp_mod_list):
    	    
	    kiw_comp=-1;
	    for kiw,iw in enumerate(imp_mod):
		if test_impedance_wake_comp(iw,a,b,c,d,plane): kiw_comp=kiw;
	    
	    if (kiw_comp>=0):
	    	
		if imod==0:
	    	    kiw0_comp=kiw_comp;
		    if plotpercent: # initialize sums for filled percent plot
			xsum=imp_mod[kiw0_comp].var;ysum=np.zeros(len(xsum),dtype=complex);

		# impedance plot
		for ir,r in enumerate(parts):
       		    plot(imp_mod[kiw_comp].var,np.abs(imp_mod[kiw_comp].func[:,ir]),r+leglist[imod],pat[ir],ylab+unit+"]",aximp[icomp],3,xlab=xlab,colr=col[imod]);

		# impedance ratio plot
		if (imod>0)and(kiw0_comp>=0):
		    for ir,r in enumerate(parts):
			imp0=np.interp(imp_mod[kiw_comp].var,imp_mod_list[0][kiw0_comp].var,imp_mod_list[0][kiw0_comp].func[:,ir]);
			ratio=imp_mod[kiw_comp].func[:,ir]/imp0;

        		if plotpercent:
		    	    sumratio=np.interp(imp_mod[kiw_comp].var,xsum,getattr(ysum,partcomp[ir]))/imp0;
			    fillplot_percent(imp_mod[kiw_comp].var,sumratio+ratio,sumratio,xlab,leglist,col[imod],axratio[2*icomp+ir])
		    	    ysum += (1j)**ir *np.interp(xsum,imp_mod[kiw_comp].var,imp_mod[kiw_comp].func[:,ir])
			    # with fill_between we have to create the legend box ourselves
			    p[ir].append(Rectangle((0, 0), 1, 1,axes=axratio[2*icomp+ir],color=col[imod]));
			    leglistless[ir].append(leglist[imod]);

			else: plot(imp_mod[kiw_comp].var,np.abs(ratio),r+leglist[imod],pat[ir],str2+" ratio w.r.t "+leglist[0],axratio[icomp],1,xlab=xlab,colr=col[imod]);

			# find max ratio
			ind=pylab.mlab.find((imp_mod[kiw_comp].var<boundmax)*(imp_mod[kiw_comp].var>boundmin));
			maxratio[icomp,imod-1,ir]=np.max(ratio[ind]);
			print r+comp+": max. ratio between",leglist[imod],"and",leglist[0],str1,":",maxratio[icomp,imod-1,ir];

	# finish plots
	aximp[icomp].set_xlim(xlim);
	if (comp.find('l')!=-1): aximp[icomp].set_ylim(yliml);
	else: aximp[icomp].set_ylim(ylim);
	end_figure(figimp[icomp],aximp[icomp],save=(len(saveimp)>0)*(saveimp+'_'+comp),legpos=legpos);
	
	if imod>0:
	
	    if not(plotpercent):
		axratio[icomp].set_xlim(xlim);
		axratio[icomp].set_ylim([0,np.ceil(np.max(maxratio[icomp,:,:])*5)/5.]);
		end_figure(figratio[icomp],axratio[icomp],save=(len(saveratio)>0)*(saveratio+'_'+comp),legpos=legpos);

	    else:
		fontP = FontProperties();
		if len(imp_mod_list)<=6: fontP.set_size('medium');
		else: fontP.set_size('small');
		for ir,r in enumerate(parts):
		    axratio[2*icomp+ir].legend(p[ir][::-1],leglistless[ir][::-1],bbox_to_anchor=legpercentpos,prop=fontP);
		    axratio[2*icomp+ir].set_xscale('log');	    
		    axratio[2*icomp+ir].set_xlim(xlim);
		    axratio[2*icomp+ir].set_ylim([0,maxpercent]);
	    	    if len(r)>0: part='_'+r;part=part.replace(', ','').replace(' ','_').replace('_part','')
		    else: part='';
		    end_figure(figratio[2*icomp+ir],axratio[2*icomp+ir],save=(len(saveratio)>0)*(saveratio+'_percent_'+comp+part),legpos=legpercentpos);
	
    return maxratio;		

