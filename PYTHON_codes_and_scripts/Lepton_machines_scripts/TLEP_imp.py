#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

if len(sys.argv)>1: lxplusbatchImp=str(sys.argv[1]);
else: lxplusbatchImp=None;
print lxplusbatchImp   

import numpy as np
import pylab,os,re
path_here=os.getcwd()+"/";
from copy import deepcopy
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from particle_param import *
from Impedance import *
from DELPHI import *

def TLEP_param(E=45e9,option='Z',Qxfrac=0.2,Qyfrac=0.15):

    ''' generate typical TLEP parameters, from the beam energy E in eV, the TLEP option 
    ('Z', 'H', 'W', 't', 'tB' or 'Z4C') and the fractional parts of the tunes.
    Outputs:
    - machine: string with machine name(here 'TLEP'+option),
    - E: same as input (beam energy in eV),
    - gamma: relativistic mass factor,
    - sigmaz: RMS bunch length in m,
    - taub: total bunch length in s (4*RMS),
    - R: machine pysical radius (circumference/(2 pi)),
    - Qx: total horizontal tune (integer + fractional parts),
    - Qxfrac: fractional horizontal tune,
    - Qy: total vertical tune (integer + fractional parts),
    - Qyfrac: fractional vertical tune,
    - Qs: synchrotron tune,
    - eta: slippage factor (alpha_p-1/gamma^2),
    - f0: revolution frequency,
    - omega0: revolution angular frequency=2pi*f0,
    - omegas: synchrotron angular frequency=Qs*omega0,
    - dphase: phase of damper w.r.t. "normal" purely resistive damper (0),
    - Estr: string with energy (e.g. '45GeV').
    - syncdamp: synchrotron transverse damping time (transverse=2*longitudinal) in seconds.
    '''

    e,m0,c,E0=electron_param();
    # E is the energy in eV
    Estr=str(int(E/1e9))+'GeV';print Estr
    machine='TLEP'+option;
    gamma=E*e/E0
    beta=np.sqrt(1.-1./(gamma**2))
    circ=100e3; # total circumference in m
    R=circ/(2.*np.pi) # machine radius
    f0=c*beta/circ # rev. frequency
    omega0=2.*np.pi*f0;
    
    # various TLEP options
    if (option=='Z'):
    	if abs(E-45e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
    	Qs=0.77e3/f0;Qs=0.25;alphap=3.6e-5;sigmaz=1.16e-3;
	#sigmaz=2.93e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP rescaled with nb FODO cells (w.r.t. TLEPH)
	Qxint=4*78/3;Qyint=4*98/3;
    
    elif (option=='Z4C'): # modified Z option with 4 cells
    	if abs(E-45e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
    	Qs=0.52e3/f0;alphap=1.6e-5;sigmaz=0.73e-3;
	#sigmaz=2.93e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP rescaled with nb FODO cells (w.r.t. TLEPH)
	Qxint=4*78/2;Qyint=4*98/2;

    
    elif (option=='W'):
    	if abs(E-80e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
    	Qs=0.19e3/f0;alphap=0.4e-5;sigmaz=0.91e-3;
	#sigmaz=1.98e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP
	Qxint=4*78;Qyint=4*98;
    
    elif (option=='H'):
    	if abs(E-120e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
	Qs=0.27e3/f0;alphap=0.4e-5;sigmaz=0.98e-3;
	#sigmaz=2.11e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP
	Qxint=4*78;Qyint=4*98;

    elif (option=='t'):
    	if abs(E-175e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
    	Qs=0.14e3/f0;alphap=0.1e-5;sigmaz=0.68e-3;
	#sigmaz=0.77e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP rescaled with nb FODO cells (w.r.t. TLEPH)
	Qxint=4*78*2;Qyint=4*98*2;
    
    elif (option=='tB'):
    	if abs(E-175e9)>1e9: print "Pb in TLEP_param: energy not correct !";sys.exit()
    	Qs=0.29e3/f0;alphap=0.4e-5;sigmaz=1.35e-3;
	#sigmaz=1.95e-3; # with beam-beam (beamstralhung)
	# integer part of tune=4*LEP rescaled with nb FODO cells (w.r.t. TLEPH)
	Qxint=4*78;Qyint=4*98;
    
    else:
    	print "Pb in TLEP_param: option not recognized !";sys.exit()

    Qx=Qxfrac+Qxint;
    Qy=Qyfrac+Qyint;
    taub=4.*sigmaz/(beta*c); # full length in s
    omegas=Qs*omega0;
    eta=alphap-1./(gamma*gamma); # slip factor
    dphase=0.; # additional damper phase
    # synchrotron damping time (transverse=2*longitudinal) in seconds 
    # (if instability slower than this -> stable)
    syncdamp=2*1319/f0;
    
    return machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,syncdamp;


def TLEP_imp(gamma,circ,avbetax,avbetay,length_dip=11.,length_absorber=0.5,
	bending_radius=11000,b=0.015,dire='TLEP',wake_calc=False,lxplusbatch=None,option='all'):

    ''' Generates TLEP impedance and wake model
    all units are SI.
    gamma is the relativistic mass factor and circ the total circumference
    avbetax and avbetay are the average beta functions (where kick will be applied)
    length_dip is the dipole length, length_absorber the photon absorber length, 
    bending_radius the dipoles bending radius, b the vertical semi-axis
    dire is the directory in ImpedanceWake2D where to put the results
    wake_calc should be True to compute wake as well
    lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus
    		   if 'retrieve' -> retrieve outputs
    option:
     -'all': compute total impedance
     -'RW': resistive-wall only
     -'RWnoabs' : resistive-wall without absorber part
     -'BB' : total broad-band part
     -'RF': only broad-band from RF
     -'abs': only absorber part (RW+BB from tapers)'''
    
    Ax=0.045; # horizontal aperture of vacuum pipe
    # compute absorber aperture from dipole & absorber lengths
    L=bending_radius*(np.arccos(bending_radius*np.cos(1./gamma)/(bending_radius+Ax))-1./gamma)-length_dip;
    absorber_ax=bending_radius*(np.cos(-1./gamma)/np.cos((L+length_absorber)/bending_radius-1./gamma)-1.);
    print absorber_ax
    
    fpar=freq_param(fmin=1,fmax=1e14,ftypescan=2,nflog=40,fminrefine=1e11,fmaxrefine=5e12,nrefine=2000)
    zpar=z_param(ztypescan=2,zmin=0.1,zmax=1e7,nzlog=20,zminrefine=2e-6,zmaxrefine=0.02,zsamplin=2e-6)
    waketol=1.e11;
    
    beta=np.sqrt(1.-1./(gamma**2))
    lendip=2*np.pi*bending_radius; # total length of dipoles
    ndip=lendip/length_dip; # number of dipoles
    
    imp_mod=[];wake_mod=[];
    
    if wake_calc: queue='1nd';
    else: queue='1nd'
    
    if (option=='all')or(option.startswith('RW')):
	# vacuum pipe
	wvac=Ax;length=circ-ndip*length_absorber;
	layers_iw=[Al_layer(thickness=np.inf)];
	betax=2*avbetax;betay=2*avbetay;

	iw_input=impedance_wake_input(gamma=gamma,length=1,b=[b],layers=layers_iw,
    	    fpar=fpar,zpar=zpar,waketol=waketol,freqlinbisect=1e11,geometry='round',comment='_vacuum_pipe');

	imp_mod_vac,wake_mod_vac=imp_model_elliptic(iw_input,wvac,orientation='V',
    	    wake_calc=wake_calc,flagrm=True,lxplusbatch=lxplusbatch,queue=queue,dire=dire);

	multiply_impedance_wake(imp_mod_vac,length);
	multiply_impedance_wake(wake_mod_vac,length);

	# add to the model
	add_impedance_wake(imp_mod,imp_mod_vac,betax/avbetax,betay/avbetay);
	add_impedance_wake(wake_mod,wake_mod_vac,betax/avbetax,betay/avbetay);


    if (option=='all')or((option=='RW')or(option=='abs')):
	# photon absorber
	wabs=absorber_ax;length=ndip*length_absorber;
	
	if (b<wabs): small_axis=b;large_axis=wabs;orientation='V';
	else: small_axis=wabs;large_axis=b;orientation='H';print "Photon-absorber: H smaller",wabs,b
	
	layers_iw=[W_layer(thickness=np.inf)]
	betax=2*avbetax;betay=2*avbetay;

	iw_input=impedance_wake_input(gamma=gamma,length=1,b=[small_axis],layers=layers_iw,
    	    fpar=fpar,zpar=zpar,waketol=waketol,freqlinbisect=1e11,geometry='round',comment='_photon_absorbers');

	imp_mod_abs,wake_mod_abs=imp_model_from_IW2D(iw_input,wake_calc=wake_calc,
    	    flagrm=True,lxplusbatch=lxplusbatch,queue=queue,dire=dire);

	imp_mod_abs=apply_Yokoya_elliptic(imp_mod_abs,small_axis,large_axis,orientation=orientation)
	wake_mod_abs=apply_Yokoya_elliptic(wake_mod_abs,small_axis,large_axis,orientation=orientation)

	multiply_impedance_wake(imp_mod_abs,length);
	multiply_impedance_wake(wake_mod_abs,length);

	# add to the model
	add_impedance_wake(imp_mod,imp_mod_abs,betax/avbetax,betay/avbetay);
	add_impedance_wake(wake_mod,wake_mod_abs,betax/avbetax,betay/avbetay);
    
    
    if (option=='all')or((option=='BB')or((option=='RF')or(option=='abs'))):
	# broad-band model
	fcutoff=6e9;Q=1;betax=2*avbetax;betay=2*avbetay;
	# tapers for photon absorbers
	ntaper=2*ndip;
	Rtaper=ntaper*transverse_imp_taper_round_Yokoya(absorber_ax,Ax,(Ax-absorber_ax)/length_absorber);
	if (option=='RF'): Rtaper=0;
	# RF cavities: from BB approximation of Rama Calaga PhD thesis results (600m)
	Rcav=600*3e3;
	if (option=='abs'): Rcav=0;
	print Rtaper,Rcav

	# resonator BB model
	imp_mod_BB,wake_mod_BB=imp_model_resonator(Rtaper+Rcav,fcutoff,Q,beta=beta,wake_calc=wake_calc,
	    fpar=fpar,zpar=zpar,listcomp=['Zxdip','Zydip']);

	# add to the model
	add_impedance_wake(imp_mod,imp_mod_BB,betax/avbetax,betay/avbetay);
	add_impedance_wake(wake_mod,wake_mod_BB,betax/avbetax,betay/avbetay);
    

    return imp_mod,wake_mod
    
    
if __name__ == "__main__":

    e,m0,c,E0=electron_param();

    # fixed parameters
    E=45e9;option='Z';
    machine,E,gamma,sigmaz,taub,R,Qx,Qxfrac,Qy,Qyfrac,Qs,eta,f0,omega0,omegas,dphase,Estr,syncdamp=TLEP_param(E=E,option=option);
    beta=np.sqrt(1.-1./(gamma**2))
    bscan=[0.02,0.025,0.03,0.035,0.04]; # vertical aperture (semi-axis) [m]
    #bscan=[0.03]; # vertical aperture (semi-axis) [m]
    #ldipscan=[5.5]; # dipole length [m]
    ldipscan=[1.,2.,3.5,5.5,8.,11.]; # dipole length [m]

    print machine,R/Qx,R/Qy,Qs
    
    for ib,b in enumerate(bscan):
    
	for ildip,ldip in enumerate(ldipscan):

	    imp_mod,wake_mod=TLEP_imp(gamma,2*np.pi*R,R/Qx,R/Qy,length_dip=ldip,
		length_absorber=0.5,bending_radius=11000.,b=b,
		dire='TLEP/',wake_calc=True,lxplusbatch=lxplusbatchImp,option='all')

	    if (lxplusbatchImp==None)or(lxplusbatchImp.startswith('ret')):
		# output in a HEADTAIL kind of wake file
	        os.system("mkdir -p "+path_here+"../../../DELPHI_results/"+machine);
		wakefilename=path_here+'../../../DELPHI_results/'+machine+'/wake_for_hdtl_'+machine+Estr+'_'+float_to_str(1000*b)+'mm_ldip'+float_to_str(ldip)+'m_all_wake.dat';
		write_wake_HEADTAIL(wake_mod,wakefilename,beta=beta,ncomp=4);

