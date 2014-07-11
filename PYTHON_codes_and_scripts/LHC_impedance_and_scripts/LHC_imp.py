#!/usr/bin/python

import sys
from string import *
import numpy as np
import pickle as pick
import pylab,re
from plot_lib import plot,init_figure,end_figure
from io_lib import *
from string_lib import *
from particle_param import *
from Impedance import *
from LHC_coll_imp import *


def LHC_element_iw_model(name,layers,b,w,E,length,orientation,thickness=[],Bfield=None,
	wake_calc=False,fpar=freq_param(),zpar=z_param(),lxplusbatch=None,comment='',dire=''):
    
    # construct impedance/wake model for an LHC elliptic element
    # name is the name of the element, layers a list of layers or material names,
    # b and w the semi-axes (resp. small and large) in m.
    # E is the energy in eV.
    # orientation:'H' or 'V' for hor. or vert. orientation of b. 
    # if material names are given (not layer class object), then layers
    # thickness should be given as well in the list 'thickness'.
    # Bfield: magnetic field (if None, compute it from the energy with LHC bending radius)
    # wake_calc: flag for wake calculation
    #
    # fpar and zpar select the frequencies and distances (for the wake) scans 
    # (objects of the classes freq_param and z_param).
    #
    # lxplusbatch: if None, no use of lxplus batch system
    # 		   if 'launch' -> launch calculation on lxplus on queue 'queue'
    #		   if 'retrieve' -> retrieve outputs
    # comment is added to the IW2D output file names
    # dire contains the directory name where to put the outputs (default='./'=directory of IW2D)
    
    e,m0,c,E0=proton_param();
    gamma=E*e/E0;

    # find on which queue to launch calculation (case when lxplusbatch=='launch')
    queues=['8nm','1nh','1nh','8nh','1nd','2nd','1nw'];
    if wake_calc: iq=3;
    else:      
	nfreq=(np.mod(fpar.ftypescan,2)==0)*((np.log10(fpar.fmax)-np.log10(fpar.fmin))*fpar.nflog+1) + (fpar.ftypescan==1)*round((fpar.fmax-fpar.fmin)/fpar.fsamplin)+(fpar.ftypescan==2)*fpar.nrefine;
	iq=max(int(np.ceil(np.log10(nfreq/500.))),0); # =0 up to 500 frequencies; 1 up to 5000; 2 up to 50000; etc.
	iq+=1;
    queue=queues[iq];
    
    print layers[0], queue;
    
    # some hard coded parameters
    freqlin=1.e11;
    
    # construct layers list, then input file
    if Bfield==None: Bfield=E/(2803.9*299792458);
    layers_iw=construct_layers(layers,thickness=thickness,Bfield=Bfield);
    
    iw_input=impedance_wake_input(gamma=gamma,length=length,b=[b],layers=layers_iw,
    	fpar=fpar,zpar=zpar,freqlinbisect=freqlin,geometry='round',comment=name+comment);
    
    imp_mod,wake_mod=imp_model_elliptic(iw_input,w,orientation=orientation,
    	wake_calc=wake_calc,flagrm=True,lxplusbatch=lxplusbatch,queue=queue,dire=dire);
    
    return imp_mod,wake_mod;
    

def LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename,beta_filename,Bfield=None,
	wake_calc=False,ftypescan=2,nflog=100,namesref=None,power_loss_param=None,
	freq_dep_factor_file=None,dist_dep_factor_file=None,
	lxplusbatch=None,comment='',dire=''):

    # creates an impedance or wake model for many elements (cold beam screens
    # and warm pipes).
    # E is the energy in eV, avbetax and avbetay the average beta functions
    # used for the weighting, param_filename is the file with all parameters
    # except betas, beta_filename the file with beta functions (in m).
    # nmat is the number of material layers -> look for columns 'material1' to 'material[nmat]'
    # in the file param_filename, and for columns 'thickness1' to thickness[nmat-1]'
    #
    # Bfield: magnetic field (if None, compute it from the energy with LHC bending radius)
    #
    # wake_calc selects the wake computation if True, nflog is the number of frequencies
    # per decade, and ftypescan the type of frequency scan (0,1 or 2: logarithmic only, 
    # linear only or logarithmic with refinement around high-frequency resonance(s) ).
    #
    # namesref are the coll. names (from param_filename) to select (if None take all),
    #
    # if power_loss_param is not None, compute (and print) the power loss (in W and W/m) for each element,
    # using the parameters given in power_loss_param, which should be like [sigz(RMS in m),Nb(nb p+),M(nb of bunches)]
    # (assumes Gaussian equidistant bunches and protons)
    #
    # if freq_dep_factor_file is not None, each component is multiplied by a 
    # frequency dependent factor, from the file indicated by this. File should
    # contain (colums): frequency in Hz, then each component (factors for real
    # part, then imaginary part), WITH HEADERS indicating the component (format:
    # see function 'identify_component' in module Impedance + at the end re or 
    # Re for real part and im or Im for imag. part).
    # It is exactly the same for the wake with dist_dep_factor_file (if wake_calc==True).
    #
    # lxplusbatch: if None, no use of lxplus batch system
    # 		   if 'launch' -> launch calculation on lxplus on queue 'queue'
    #		   if 'retrieve' -> retrieve outputs
    # dire contains the directory name where to put the outputs (default='./'=directory of IW2D)

    # read files
    # file with materials, orientation and thickness
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);

    orientation=read_ncol_file_identify_header(param_filename,'[oO]rientation');
    b=read_ncol_file_identify_header(param_filename,'[sS]mall_axis');
    w=read_ncol_file_identify_header(param_filename,'[lL]arge_axis');
    orientation=[orientation[i] for i in ind];
    b=b[ind];w=w[ind];
    
    # read materials an thicknesses
    material,thick,nmat=read_materials(param_filename,ind);
    print material,thick,nmat
    
    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];length=length[ind];
    
    # main loop to construct model
    imp_mod=[];wake_mod=[];
    for iname,name in enumerate(namesref):
    
	fminrefine=1.e11;fmaxrefine=5.e12;nrefine=5000;
	
	materials=[material[n][iname] for n in range(nmat)];
	thicks=[thick[n][iname] for n in range(0,nmat-1)];thicks.append(np.inf);
	
	imp,wake=LHC_element_iw_model(name,materials,
		b[iname],w[iname],E,length[iname],orientation[iname],
		thickness=thicks,Bfield=Bfield,wake_calc=wake_calc,
		fpar=freq_param(ftypescan=ftypescan,nflog=nflog,
		fminrefine=fminrefine,fmaxrefine=fmaxrefine,nrefine=nrefine),
		lxplusbatch=lxplusbatch,comment=comment,dire=dire);
		
	if (power_loss_param!=None):
	    # compute power loss for this element
	    e,m0,c,E0=proton_param();
	    gamma=E*e/E0;
	    sigz=power_loss_param[0];Nb=power_loss_param[1];M=power_loss_param[2];
	    ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
	    print name,comment,", without factor: power loss=",ploss,"W, power loss per unit length=",ploss/length[iname],"W/m";

	if (freq_dep_factor_file!=None):
	    # apply freq. dependent factor to imp (typically from weld)
	    imp=multiply_impedance_wake_file(imp,freq_dep_factor_file);
	    if (power_loss_param!=None):
		# compute power loss for this element (with the factor this time)
		ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
		print name,comment,", with factor from "+freq_dep_factor_file+":";
		print "    power loss=",ploss,"W, power loss per unit length=",ploss/length[iname],"W/m";
	
	if (dist_dep_factor_file!=None):
	    # apply freq. dependent factor to wake (typically from weld)
	    wake=multiply_impedance_wake_file(wake,dist_dep_factor_file);
	
	add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);
	    
    return imp_mod,wake_mod;
	

def LHC_design_Broadband(squeeze=True,wake_calc=False,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param()):

    # LHC broad-band model (longitudinal and transverse)
    # based on LHC design report (chap. 5)
    # two different cases: squeezed optics or injection optics
    # fpar and zpar select the frequencies and
    # distances (for the wake) scans (objects of the classes freq_param and z_param).
       
    Q=1;fr=5e9; # quality factor and cutoff
    # shunt impedances (Ohm in long., Ohm/m in transverse)
    if squeeze: Rl=33.8e3;Rt=2.67e6;
    else: Rl=31.1e3;Rt=1.34e6;
    
    imp_trans,wake_trans=imp_model_transverse_resonator(Rt,fr,Q,beta=1,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar);
    
    imp_long,wake_long=imp_model_longitudinal_resonator(Rl,fr,Q,beta=1,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar);
    
    # total impedance and wake (weights=1 -> assume broad-band imp. at average beta functions)
    imp_mod=[];wake_mod=[];
    
    add_impedance_wake(imp_mod,imp_trans,1,1);
    add_impedance_wake(imp_mod,imp_long,1,1);

    add_impedance_wake(wake_mod,wake_trans,1,1);
    add_impedance_wake(wake_mod,wake_long,1,1);

    return imp_mod,wake_mod;
    

def LHC_design_Broadband_less(wake_calc=False,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),fcutoffBB=5e9):

    ''' LHC broad-band model (longitudinal and transverse)
    based on LHC design report (chap. 5), taking out collimators, pumping slots and RF.
    two different cases: squeezed optics or injection optics
    fpar and zpar select the frequencies and
    distances (for the wake) scans (objects of the classes freq_param and z_param).
    fcutoffBB: cutoff frequency of the broad-band model
    '''
       
    Q=1;fr=fcutoffBB; # quality factor and cutoff
    # shunt impedances (Ohm in long., Ohm/m in transverse)
    Rl=10.5e3*fcutoffBB/5e9; # longitudinal shunt impedance is scaled by cutoff frequency
    Rt=0.67e6; #(no space-charge here)
    
    imp_trans,wake_trans=imp_model_transverse_resonator(Rt,fr,Q,beta=1,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar);
    
    imp_long,wake_long=imp_model_longitudinal_resonator(Rl,fr,Q,beta=1,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar);
    
    # total impedance and wake (weights=1 -> assume broad-band imp. at average beta functions)
    imp_mod=[];wake_mod=[];
    
    add_impedance_wake(imp_mod,imp_trans,1,1);
    add_impedance_wake(imp_mod,imp_long,1,1);

    add_impedance_wake(wake_mod,wake_trans,1,1);
    add_impedance_wake(wake_mod,wake_long,1,1);

    return imp_mod,wake_mod;
    

def LHC_manyBB_resonator(avbetax,avbetay,param_filename,beta_filename,
	fcutoff=5e9,Q=1,beta=1,wake_calc=False,namesref=None,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=None):
	
    # add up many axisymmetric broad-band resonators into a single model.
    # only longitudinal and dipolar terms are considered.
    # avbetax and avbetay are the average beta functions
    # used for the weighting, param_filename is the file with all parameters
    # except betas, beta_filename the file with beta functions (in m).
    # fcutoff is the cutoff frequency of the broad-band model (from the beam pipe radius)
    # and Q the quality factor chosen (usually 1).
    # beta is the relativistic velocity factor
    # wake_calc selects the wake computation if True, fpar and zpar select the frequencies and
    # distances (for the wake) scans (objects of the classes freq_param and z_param).
    # namesref are the coll. names (from param_filename) to select (if None take all).
    # if power_loss_param is not None, compute (and print) the power loss (in W and W/m) for each element,
    # using the parameters given in power_loss_param, which should be like
    # [sigz(RMS in m),gamma,Nb(nb p+),M(nb of bunches)]
    # (assumes Gaussian equidistant bunches and protons)
    
    # read files
    # file with materials, orientation and thickness
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);

    Rt=read_ncol_file_identify_header(param_filename,'[rR][tT]');
    Rl=read_ncol_file_identify_header(param_filename,'[rR][lL]');
    length_ind=read_ncol_file_identify_header(param_filename,'[lL]ength'); # length of individual elements
    Rt=Rt[ind];Rl=Rl[ind];length_ind=length_ind[ind];
    
    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];length=length[ind]
    
    if (power_loss_param!=None):
    	for iname,name in enumerate(namesref):
	    # compute power loss for each element
	    imp,wake=imp_model_longitudinal_resonator(Rl[iname],fcutoff,Q,beta=beta,wake_calc=wake_calc,fpar=fpar,zpar=zpar)
	    sigz=power_loss_param[0];gamma=power_loss_param[1];
	    Nb=power_loss_param[2];M=power_loss_param[3];
	    ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
	    print name,": power loss=",ploss,"W, power loss per unit length=",ploss/length_ind[iname],"W/m";

    # sum all broad-band contributions with appropriate weighting
    Rtranstotx=np.sum(Rt*length*betax/(avbetax*length_ind));
    Rtranstoty=np.sum(Rt*length*betay/(avbetay*length_ind));
    Rlongtot=np.sum(Rl*length/length_ind);
    
    # create impedance/wake model
    imp_modx,wake_modx=imp_model_transverse_resonator(Rtranstotx,fcutoff,Q,beta=beta,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar,plane='x')
    imp_mody,wake_mody=imp_model_transverse_resonator(Rtranstoty,fcutoff,Q,beta=beta,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar,plane='y')
    imp_modl,wake_modl=imp_model_longitudinal_resonator(Rlongtot,fcutoff,Q,beta=beta,wake_calc=wake_calc,
	fpar=fpar,zpar=zpar)
	
    # add up to get total model
    imp_mod=[];wake_mod=[];
    
    add_impedance_wake(imp_mod,imp_modx,1,1);
    add_impedance_wake(imp_mod,imp_mody,1,1);
    add_impedance_wake(imp_mod,imp_modl,1,1);

    add_impedance_wake(wake_mod,wake_modx,1,1);
    add_impedance_wake(wake_mod,wake_mody,1,1);
    add_impedance_wake(wake_mod,wake_modl,1,1);

    return imp_mod,wake_mod;
    

def LHC_many_resonators(avbetax,avbetay,param_filename,beta_filename,
	beta=1,wake_calc=False,namesref=None,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=None):
	
    # add up many resonators into a single model.
    # each device is defined by an HOM file. 
    # avbetax and avbetay are the average beta functions
    # used for the weighting, param_filename is the file with all parameters
    # except betas, beta_filename the file with beta functions (in m).
    # beta is the relativistic velocity factor
    # wake_calc selects the wake computation if True, fpar and zpar select the frequencies and
    # distances (for the wake) scans (objects of the classes freq_param and z_param).
    # namesref are the coll. names (from param_filename) to select (if None take all).
    # if power_loss_param is not None, compute (and print) the power loss (in W and W/m) for each element,
    # using the parameters given in power_loss_param, which should be like
    # [sigz(RMS in m),gamma,Nb(nb p+),M(nb of bunches)]
    # (assumes Gaussian equidistant bunches and protons)
    
    # read files
    # file with HOMs file names and individual length of each device
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);

    filenames=read_ncol_file_identify_header(param_filename,'[fF]ile');
    filenames=np.array(filenames);
    length_ind=read_ncol_file_identify_header(param_filename,'[lL]ength'); # length of individual elements
    filenames=filenames[ind];length_ind=length_ind[ind];
    
    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];length=length[ind]
    
    # add up to get total model
    imp_mod=[];wake_mod=[];
    
    for iname,name in enumerate(namesref):
    
    	imp,wake=imp_model_from_HOMfile(filenames[iname],beta=beta,fpar=fpar,zpar=zpar);
    	
	if (power_loss_param!=None):
	    # compute power loss for each individual element
	    sigz=power_loss_param[0];gamma=power_loss_param[1];
	    Nb=power_loss_param[2];M=power_loss_param[3];
	    ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
	    print name,": power loss=",ploss,"W, power loss per unit length=",ploss/length_ind[iname],"W/m";

	# multiply imp & wake by number of such elements and add (with beta function weight)
	# to the model
	multiply_impedance_wake(imp,length[iname]/length_ind[iname]);
	multiply_impedance_wake(wake,length[iname]/length_ind[iname]);
	add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);
    

    return imp_mod,wake_mod;
    

def LHC_many_striplineBPMs(avbetax,avbetay,param_filename,beta_filename,
	beta=1,wake_calc=False,namesref=None,
	fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=None):
	
    # add up many stripline BPMS into a single model (using Ng simple fomulas).
    # each device is defined by an HOM file. 
    # avbetax and avbetay are the average beta functions
    # used for the weighting, param_filename is the file with all parameters
    # except betas, beta_filename the file with beta functions (in m).
    # beta is the relativistic velocity factor
    # wake_calc selects the wake computation if True, fpar and zpar select the frequencies and
    # distances (for the wake) scans (objects of the classes freq_param and z_param).
    # namesref are the coll. names (from param_filename) to select (if None take all).
    # if power_loss_param is not None, compute (and print) the power loss (in W and W/m) for each element,
    # using the parameters given in power_loss_param, which should be like
    # [sigz(RMS in m),gamma,Nb(nb p+),M(nb of bunches)]
    # (assumes Gaussian equidistant bunches and protons)
    
    # read files
    # file with stripline BPMs parameters and individual length of each device
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);

    strips_length=read_ncol_file_identify_header(param_filename,'[sS]trip'); # strips length
    angle=read_ncol_file_identify_header(param_filename,'[aA]ngle'); # azimuthal angle under which it is seen
    b=read_ncol_file_identify_header(param_filename,'[bB]'); # chamber radius
    length_ind=read_ncol_file_identify_header(param_filename,'[lL]ength'); # length of individual elements
    strips_length=strips_length[ind];length_ind=length_ind[ind];
    angle=angle[ind];b=b[ind];
    
    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];length=length[ind]
    
    # add up to get total model
    imp_mod=[];wake_mod=[];
    
    for iname,name in enumerate(namesref):
    
    	imp,wake=imp_model_striplineBPM_Ng(strips_length[iname],angle[iname],b[iname],beta=beta,
		wake_calc=wake_calc,fpar=fpar,zpar=zpar);
    	
	if (power_loss_param!=None):
	    # compute power loss for each individual element
	    sigz=power_loss_param[0];gamma=power_loss_param[1];
	    Nb=power_loss_param[2];M=power_loss_param[3];
	    ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
	    print name,": power loss=",ploss,"W, power loss per unit length=",ploss/length_ind[iname],"W/m";

	# multiply imp & wake by number of such elements and add (with beta function weight)
	# to the model
	multiply_impedance_wake(imp,length[iname]/length_ind[iname]);
	multiply_impedance_wake(wake,length[iname]/length_ind[iname]);
	add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);
    

    return imp_mod,wake_mod;
    

def LHC_many_holes(avbetax,avbetay,param_filename,beta_filename,
	fcutoff=5e9,namesref=None,fpar=freq_param(ftypescan=0,nflog=100),zpar=z_param(),power_loss_param=None):
	
    # add up many axisymmetric elements with pumping holes into a single model.
    # only longitudinal and dipolar terms are considered.
    # avbetax and avbetay are the average beta functions
    # used for the weighting, param_filename is the file with all parameters
    # except betas, beta_filename the file with beta functions (in m).
    # fcutoff is the cutoff frequency (from the beam pipe radius, default value=LHC).
    # fpar selects the frequencies scan (object of the classes freq_param).
    # namesref are the coll. names (from param_filename) to select (if None take all).
    # if power_loss_param is not None, compute (and print) the power loss (in W and W/m) for each element,
    # using the parameters given in power_loss_param, which should be like
    # [sigz(RMS in m),gamma,Nb(nb p+),M(nb of bunches)]
    # (assumes Gaussian equidistant bunches and protons)
    
    # read files
    # file with all parameters
    if (namesref==None): namesref=read_ncol_file_identify_header(param_filename,'[nN]ame');
    names=read_ncol_file_identify_header(param_filename,'[nN]ame');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);

    Lh=read_ncol_file_identify_header(param_filename,'[lL][hH]');
    Wh=read_ncol_file_identify_header(param_filename,'[wW][hH]');
    T=read_ncol_file_identify_header(param_filename,'[tT]');
    b=read_ncol_file_identify_header(param_filename,'[bB]');
    d=read_ncol_file_identify_header(param_filename,'[dD]');
    eta=read_ncol_file_identify_header(param_filename,'[eE]ta');
    rhob=read_ncol_file_identify_header(param_filename,'[rR]hob');
    rhod=read_ncol_file_identify_header(param_filename,'[rR]hod');
    length_ind=read_ncol_file_identify_header(param_filename,'[lL]ength'); # length of individual elements
    nb_holes_cs=read_ncol_file_identify_header(param_filename,'[nN]b_holes');
    
    Lh=Lh[ind];Wh=Wh[ind];T=T[ind];b=b[ind];d=d[ind];eta=eta[ind];
    rhob=rhob[ind];rhod=rhod[ind];length_ind=length_ind[ind];nb_holes_cs=nb_holes_cs[ind];
    
    # file with beta functions
    names=read_ncol_file_identify_header(beta_filename,'[nN]ame');
    betax=read_ncol_file_identify_header(beta_filename,'[bB]etax');
    betay=read_ncol_file_identify_header(beta_filename,'[bB]etay');
    length=read_ncol_file_identify_header(beta_filename,'[lL]ength');
    # reorder such that the names match with namesref
    ind=find_ind_names(namesref,names);
    betax=betax[ind];betay=betay[ind];length=length[ind]
    
    # main loop to construct model
    imp_mod=[];wake_mod=[];
    for iname,name in enumerate(namesref):
	
	imp,wake=imp_model_holes_Mostacci(Lh[iname],Wh[iname],T[iname],b[iname],d[iname],
		eta[iname],rhob[iname],rhod[iname],length_ind[iname],
		fpar=freq_param(ftypescan=0,nflog=100),zpar=zpar,fcutoff=fcutoff,
		nb_holes_per_crosssection=nb_holes_cs[iname],Cm=1.,Ce=1.)
	
	if (power_loss_param!=None):
	    # compute power loss for each element
	    sigz=power_loss_param[0];gamma=power_loss_param[1];
	    Nb=power_loss_param[2];M=power_loss_param[3];
	    ploss=power_loss(imp,sigz,gamma,Nb,M,26658.883);
	    print name,": power loss=",ploss,"W, power loss per unit length=",ploss/length_ind[iname],"W/m";

	# multiply imp & wake by number of such elements and add (with beta function weight)
	# to the model
	multiply_impedance_wake(imp,length[iname]/length_ind[iname]);
	multiply_impedance_wake(wake,length[iname]/length_ind[iname]);
	add_impedance_wake(imp_mod,imp,betax[iname]/avbetax,betay[iname]/avbetay);
	add_impedance_wake(wake_mod,wake,betax[iname]/avbetax,betay[iname]/avbetay);

    return imp_mod,wake_mod;
    

def select_LHC_names(names,pattern=''):
    # select from a list of LHC device names those that begin with "pattern"
    # By default, pattern is such that all names are selected.
    
    namesnew=[];
    for name in names:
	if (name.startswith(pattern)): namesnew.append(name);
	
    return namesnew;
    

def compare_imp_vs_zbase(imp_mod,root_zbase="../Impedances/PostLS1/total_impedances/",
	suffix_zbase="_Allthemachine_7TeV_B1_postLS1_baseline.dat",save=None,
	xlim=[1e2,1e12],ylim=[1e-3,1e10],wake_flag=False,ncomp=None):

    ''' read zbase impedances and compare with impedance model in imp_mod
    root_zbase and suffix_zbase are resp. the directory and suffix of
    zbase impedance files.
    to save plots: save = root of filename where to save
    xlim and ylim are the axes limits (only for impedance plot for ylim)
    
    if wake_flag==True, compare wake components instead of impedances
    Note: xlim and ylim have to be changed accordingly (e.g. resp. [1e-5,1e6] and [1e8,1e19] for LHC total wakes)
    
    if ncomp is not None, it is the number of compoments to compare (see listcomp)

    exemples: 
    root_zbase="/afs/cern.ch/user/z/zdata/public/zbase/data2/LHC/2012/Coll_BS_Warmpipe_MBW_MQW_BB_newwakesbeta_modelApril2012/";
    suffix_zbase="_AllColl_4TeV_B1.dat"; # strange: this case is not OK
    suffix_zbase="_AllColl_4TeV_B1_20120624_TCSGclosed.dat"; # perfect benchmark
    suffix_zbase="_Allthemachine_4TeV_B1_physics_fill_3265.dat"; # perfect benchmark (except long. -> normal: BB not implemented in long.)
    suffix_zbase="_TCP_4TeV_B1.dat"; # perfect benchmark
    root_zbase="../Impedances/PostLS1/total_impedances/";
    suffix_zbase="_Allthemachine_7TeV_B1_postLS1_baseline.dat"; # 
    suffix_zbase="_Allthemachine_7TeV_B1_postLS1_relaxed.dat"; # 
    suffix_zbase="_Allthemachine_7TeV_B1_postLS1_veryrelaxed.dat"; #
    '''

    if wake_flag:
        completter='W';parts=[''];strtype="Wake";
	xlab="Distance behind the source [m]";ylab="Wake [V/C";
    else:
    	completter='Z';parts=['Re','Im'];strtype="Imp.";
	xlab="Frequency [Hz]";ylab="Z [$\Omega$";
    
    listcomp=['long','xxdip','yydip','xxquad','yyquad','xydip','xyquad','xydip','xyquad'];
    listcomp=[completter+compname for compname in listcomp];
    # corresponding a,b,c,d and planes in imp_mod
    lista=[0,1,0,0,0,0,0,1,0];listb=[0,0,1,0,0,1,0,0,0];listc=[0,0,0,1,0,0,0,0,1];listd=[0,0,0,0,1,0,1,0,0];
    listplane=['z','x','y','x','y','x','x','y','y'];
    if (ncomp==None): ncomp=len(listcomp);

    # plot to compare with zbase
    pat=['-','--'];units=["","/m","/m$^2$","/m$^3$","/m$^4$"];
    for icomp,comp in enumerate(listcomp[:ncomp]):
        #fzbase,Zzbase=readZ(root_zbase+comp+suffix_zbase);
        s=read_ncol_file(root_zbase+comp+suffix_zbase);
	fzbase=s[:,0];Zzbase=s[:,1:];
	unit=units[int(lista[icomp]+listb[icomp]+listc[icomp]+listd[icomp])];
	# find corresponding term in imp_mod
	for kiw,iw in enumerate(imp_mod):
	    if test_impedance_wake_comp(iw,lista[icomp],listb[icomp],listc[icomp],listd[icomp],listplane[icomp]): kiw_comp=kiw;

        # plot impedance
	fig,ax=init_figure();
	for ir,r in enumerate(parts):
            plot(fzbase,Zzbase[:,ir],r+'('+comp+'), Zbase',pat[ir]+'b',ylab+unit+"]",ax,3,xlab=xlab);
            plot(imp_mod[kiw_comp].var,imp_mod[kiw_comp].func[:,ir],r+' '+comp+', PyZbase','x'+pat[ir]+'r',ylab+unit+"]",ax,3,xlab=xlab);

	ax.set_xlim(xlim);ax.set_ylim(ylim);
    	if (save==None): end_figure(fig,ax);
	else: end_figure(fig,ax,save=save+'_'+comp);

        # plot impedance ratio
	fig,ax=init_figure();
	for ir,r in enumerate(parts):
	    impzbase=np.interp(imp_mod[kiw_comp].var,fzbase,Zzbase[:,ir]);
            plot(imp_mod[kiw_comp].var,imp_mod[kiw_comp].func[:,ir]/impzbase,r+' '+comp+' ',pat[ir]+'r',strtype+" ratio w.r.t. Zbase",ax,1,xlab=xlab);

	ax.set_xlim(xlim);ax.set_ylim([-1,3]);
    	if (save==None): end_figure(fig,ax);
	else: end_figure(fig,ax,save=save+'_ratio_'+comp);


def LHC_imp_model_v1(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll=None,dire="../LHC_elements/",commentcoll='',direcoll='Coll/',lxplusbatch=None,
	beam='1',squeeze='0p6m_3m_0p6m_3m',wake_calc=False,ftypescan=0,nflog=100,zpar=z_param(),
	flagplot=False,root_result='../DELPHI_results/LHC',commentsave=''):

    ''' Total LHC impedance model as of 2012 -> oct. 2013
    - E: particle energy in eV
    - avbetax, avbetay: average beta function from smooth approximation (R/Q) in m
    - param_filename_coll: file with collimator parameters (including materials, 
    angle, length, beta functions)
    - settings_filename_coll: collimator settings file
    - beta_filename_coll: collimator beta functions file (if present - otherwise take param_filename_coll)
    - dire: directory where to find all files with parameters for the rest of the machine
    - commentcoll: comment for the collimators IW2D computation
    - direcoll: subdirectory of ImpedanceWake2D where to put collimators impedance
    - lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus
    		   if 'retrieve' -> retrieve outputs
    - beam: beam number ('1' or '2')
    - squeeze: suffix of filename with beta functions, for the rest of the 
    machine (everything except collimators). The first number gives the squeeze in IP1, 
    from which the broad-band impedance (from design report) is evaluated (squeezed if <2m,
    otherwise injection BB model is taken)
    - wake_calc: True to compute wake function as well
    - ftypescan, nflog and zpar: parameters for frequency and distance scans
    - flagplot: if True, plot impedances and percent of each part of the model
    - root_result: used only with flagplot: directory where to put impedance plots.
    - commentsave: used only with flagplot: additional comment for filename of impedance plots.
    '''

    imp_mod=[];wake_mod=[];
    
    # compute model for collimators
    if beta_filename_coll==None: beta_filename_coll=param_filename_coll;
    
    imp_mod_coll,wake_mod_coll=LHC_manycoll_iw_model(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=None,lxplusbatch=lxplusbatch,
	comment=commentcoll,dire=direcoll);
    
    
    # compute model for the rest
    if E>=1e12: Estr=float_to_str(E/1e12)+'TeV';
    else: Estr=float_to_str(E/1e9)+'GeV';
    param_filename_rest=dire+"beam_screens_warm_pipe_LHC_param.dat"
    beta_filename_rest=dire+"beam_screens_warm_pipe_LHC_beta_length_B"+str(beam)+"_sq"+squeeze+".dat"
    
    imp_mod_rest,wake_mod_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_rest,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=None,
	lxplusbatch=lxplusbatch,comment='_'+Estr,dire='Rest_v2_'+Estr+'/');
    
    
    # compute broad-band model
    flagsqueeze=(float(squeeze[:squeeze.find('m')].replace('p','.'))<2); # test if squeeze below 2m (arbitrary value)
    
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband(squeeze=flagsqueeze,wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    

    # add up
    imp_mod=[];wake_mod=[];
    add_impedance_wake(imp_mod,imp_mod_coll,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll,1,1);
    add_impedance_wake(imp_mod,imp_mod_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_rest,1,1);
    add_impedance_wake(imp_mod,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod,wake_mod_BB,1,1);
    
    if flagplot:
    	imp_mod_list=[imp_mod,imp_mod_coll,imp_mod_rest,imp_mod_BB];

	leg=['Total','RW from coll','RW from beam-screen & warm pipe',
		'Broad-band contributions (design report)'];	
	
	plot_compare_imp_model(imp_mod_list,leg,listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad'],
	    saveimp=root_result+'/plot_imp_LHC_v1_'+commentsave+'details',
	    saveratio=root_result+'/plot_imp_ratio_LHC_v1_'+commentsave+'details',
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3,plotpercent=True);
	
	# plot of the total
	plot_compare_imp_model([imp_mod],[''],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad','Zxcst','Zycst'],
	    saveimp=root_result+'/plot_imp_LHC_v1_'+commentsave,
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3);

	if wake_calc:

    	    wake_mod_list=[wake_mod,wake_mod_coll,wake_mod_rest,wake_mod_BB];

	    plot_compare_imp_model(wake_mod_list,leg,listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad'],
		saveimp=root_result+'/plot_wake_LHC_v1_'+commentsave+'details',
		saveratio=root_result+'/plot_wake_ratio_LHC_v1_'+commentsave+'details',
		xlim=[1e-1,1e6],ylim=[1e12,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3,plotpercent=True);

 	    # plot of the total
	    plot_compare_imp_model([wake_mod],[''],listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad','Wxcst','Wycst'],
		saveimp=root_result+'/plot_wake_LHC_v1_'+commentsave,
		xlim=[1e-2,1e6],ylim=[1e10,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3);
	
    return imp_mod,wake_mod;
    

def LHC_imp_model_v2(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll=None,TDIcoating='preLS1',dire="../LHC_elements/",commentcoll='',direcoll='Coll_v2/',lxplusbatch=None,
	BPM=False,beam='1',squeeze='0p6m_3m_0p6m_3m',wake_calc=False,ftypescan=0,nflog=100,zpar=z_param(),
	fcutoffBB=50e9,flagplot=False,root_result='../DELPHI_results/LHC',commentsave='',
	assymetry_factor_TCL6=1.):

    ''' Total LHC impedance model as of Nov. 2013
    - E: particle energy in eV
    - avbetax, avbetay: average beta function from smooth approximation (R/Q) in m
    - param_filename_coll: file with collimator parameters (including materials, 
      angle, length, beta functions)
    - settings_filename_coll: collimator settings file
    - beta_filename_coll: collimator beta functions file (if present - otherwise take param_filename_coll)
    - TDIcoating: kind of coating for first block of TDI: can be 'preLS1'
      (5mum Ti) or 'postLS1' (1mum NEG + 2mum CU + 0.3mum NEG + 5mum Ti) or 
      directly a list of layer objects
    - dire: directory where to find all files with parameters for the rest of the machine
    - commentcoll: comment for the collimators IW2D computation
    - direcoll: subdirectory of ImpedanceWake2D where to put collimators impedance
    - lxplusbatch: if None, no use of lxplus batch system
     		   if 'launch' -> launch calculation on lxplus
    		   if 'retrieve' -> retrieve outputs
    - BPM: coll. geometry contains a BPM cavity if BPM is True, otherwise old LHC coll. geometry
    - beam: beam number ('1' or '2')
    - squeeze: suffix of filename with beta functions, for the rest of the 
    machine (everything except collimators). The first number gives the squeeze in IP1, 
    from which the broad-band impedance (from design report) is evaluated (squeezed if <2m,
    otherwise injection BB model is taken)
    - wake_calc: True to compute wake function as well
    - ftypescan, nflog and zpar: parameters for frequency and distance scans,
    - fcutoffBB: cutoff frequency for broad-band models,
    - flagplot: if True, plot impedances and percent of each part of the model,
    - root_result: used only with flagplot: directory where to put impedance plots,
    - commentsave: used only with flagplot: additional comment for filename of impedance plots.
    - assymetry_factor_TCL6: assymetry factor between the distance beam-jaw
    in the case of the TCL6 (1 means that jaws are symmetric - default).
    '''

    imp_mod=[];wake_mod=[];
    
    # compute model for collimators
    if beta_filename_coll==None: beta_filename_coll=param_filename_coll;
    
    imp_mod_coll_RW,wake_mod_coll_RW,imp_mod_coll_geom,wake_mod_coll_geom=LHC_manycoll_iw_model_with_geom(E,avbetax,avbetay,param_filename_coll,settings_filename_coll,
	beta_filename_coll,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,zpar=zpar,namesref=None,
	TDIcoating=TDIcoating,BPM=BPM,fcutoffBB=fcutoffBB,lxplusbatch=lxplusbatch,comment=commentcoll,dire=direcoll,
	assymetry_factor_TCL6=assymetry_factor_TCL6);
    
    # beta functions for all the rest
    beta_filename_rest=dire+"LHC_beta_length_B"+str(beam)+"_sq"+squeeze+".dat"

    # compute model for the rest of the RW, with weld factor for beam screens
    if E>=1e12: Estr=float_to_str(E/1e12)+'TeV';
    else: Estr=float_to_str(E/1e9)+'GeV';
    param_filename_RW=dire+"LHC_RW_param.dat"
    weld_filename=dire+"from_Carlo_weld_factor/presentscreen/weld_factor_current_arcBS_from_CZannini.dat";
    namesRW=read_ncol_file_identify_header(param_filename_RW,'name');
    namesBS=select_LHC_names(namesRW,pattern='BS'); # only BS
    namesrest=invert_selection(namesRW,namesBS) # the rest
    
    imp_mod_RW_BS,wake_mod_RW_BS=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesBS,
	freq_dep_factor_file=weld_filename,
	lxplusbatch=lxplusbatch,comment='_'+Estr,dire='BS_v2_'+Estr+'/');

    imp_mod_RW_rest,wake_mod_RW_rest=LHC_manyelem_iw_model(E,avbetax,avbetay,param_filename_RW,
    	beta_filename_rest,wake_calc=wake_calc,ftypescan=ftypescan,nflog=nflog,namesref=namesrest,
	lxplusbatch=lxplusbatch,comment='_'+Estr,dire='Rest_v2_'+Estr+'/');
    
    # individual broad-band contributions
    param_filename_BB=dire+'triplets_LHC_BB_param.dat';
    namesBB=read_ncol_file_identify_header(param_filename_BB,'name');
    namestaper=select_LHC_names(namesBB,pattern='taper');
    namesBPMs=select_LHC_names(namesBB,pattern='BPM');
    
    imp_mod_triplets_BB_taper,wake_mod_triplets_BB_taper=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB,
    	beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=namestaper,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    
    imp_mod_triplets_BB_BPMs,wake_mod_triplets_BB_BPMs=LHC_manyBB_resonator(avbetax,avbetay,param_filename_BB,
    	beta_filename_rest,fcutoff=fcutoffBB,Q=1,beta=1,wake_calc=wake_calc,namesref=namesBPMs,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    
    # HOMs contributions
    # WARNING: check RF, and check convention for shunt impedance (omega/c...)
    param_filename_HOMs=dire+'LHC_HOMs_param.dat';
    imp_mod_HOMs,wake_mod_HOMs=LHC_many_resonators(avbetax,avbetay,param_filename_HOMs,
    	beta_filename_rest,beta=1,wake_calc=wake_calc,namesref=None,
	fpar=freq_param(ftypescan=1,fmin=3e8,fmax=3e9,fsamplin=2e4,fadded=np.concatenate((10**np.arange(2,7,0.2),np.arange(1e7,5e10+1e7,1e7),[1e12]))),
	zpar=zpar);
    

    # holes contributions (note: no wake yet from long. imp contribution due to wave behind the screen (probably small))
    param_filename_holes=dire+'LHC_pumpingholes_param.dat';
    imp_mod_holes,wake_mod_holes=LHC_many_holes(avbetax,avbetay,param_filename_holes,
    	beta_filename_rest,fcutoff=fcutoffBB,namesref=None,
	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar);
    
    # compute broad-band model from design (the contributions not already taken into account elsewhere)
    imp_mod_BB,wake_mod_BB=LHC_design_Broadband_less(wake_calc=wake_calc,
    	fpar=freq_param(ftypescan=ftypescan,nflog=nflog),zpar=zpar,fcutoffBB=fcutoffBB);

    # add up
    imp_mod=[];wake_mod=[];
    add_impedance_wake(imp_mod,imp_mod_coll_RW,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll_RW,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_coll_geom,1,1);
    add_impedance_wake(wake_mod,wake_mod_coll_geom,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_RW_BS,1,1);
    add_impedance_wake(wake_mod,wake_mod_RW_BS,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_RW_rest,1,1);
    add_impedance_wake(wake_mod,wake_mod_RW_rest,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_triplets_BB_taper,1,1);
    add_impedance_wake(wake_mod,wake_mod_triplets_BB_taper,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_triplets_BB_BPMs,1,1);
    add_impedance_wake(wake_mod,wake_mod_triplets_BB_BPMs,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_HOMs,1,1);
    add_impedance_wake(wake_mod,wake_mod_HOMs,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_holes,1,1);
    add_impedance_wake(wake_mod,wake_mod_holes,1,1);
    
    add_impedance_wake(imp_mod,imp_mod_BB,1,1);
    add_impedance_wake(wake_mod,wake_mod_BB,1,1);
        
    if flagplot:
    	imp_mod_list=[imp_mod,imp_mod_coll_RW,imp_mod_coll_geom,imp_mod_RW_BS,
		imp_mod_RW_rest,imp_mod_triplets_BB_taper,
		imp_mod_triplets_BB_BPMs,imp_mod_HOMs,imp_mod_holes,imp_mod_BB];

	leg=['Total','RW from coll','Geom. from coll','RW from beam-screen',
		'RW from warm pipe','Tapers in triplets','BPMs in triplets',
		'RF, CMS, ALICE & LHCb','Pumping holes','Other broad-band contributions'];	
	
	plot_compare_imp_model(imp_mod_list,leg,listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad'],
	    saveimp=root_result+'/plot_imp_LHC_v2_'+commentsave+'details',
	    saveratio=root_result+'/plot_imp_ratio_LHC_v2_'+commentsave+'details',
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3,plotpercent=True);
	
	# plot of the total
	plot_compare_imp_model([imp_mod],[''],listcomp=['Zlong','Zxdip','Zydip','Zxquad','Zyquad','Zxydip','Zxyquad','Zxcst','Zycst'],
	    saveimp=root_result+'/plot_imp_LHC_v2_'+commentsave,
	    xlim=[1e3,5e9],ylim=[1e5,1e10],bounds=[40e6,2e9],legpos=3);

	if wake_calc:

    	    wake_mod_list=[wake_mod,wake_mod_coll_RW,wake_mod_coll_geom,wake_mod_RW_BS,
		    wake_mod_RW_rest,wake_mod_triplets_BB_taper,
		    wake_mod_triplets_BB_BPMs,wake_mod_HOMs,wake_mod_holes,wake_mod_BB];

	    plot_compare_imp_model(wake_mod_list,leg,listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad'],
		saveimp=root_result+'/plot_wake_LHC_v2_'+commentsave+'details',
		saveratio=root_result+'/plot_wake_ratio_LHC_v2_'+commentsave+'details',
		xlim=[1e-1,1e6],ylim=[1e12,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3,plotpercent=True);
	
 	    # plot of the total
	    plot_compare_imp_model([wake_mod],[''],listcomp=['Wlong','Wxdip','Wydip','Wxquad','Wyquad','Wxydip','Wxyquad','Wxcst','Wycst'],
		saveimp=root_result+'/plot_wake_LHC_v2_'+commentsave,
		xlim=[1e-2,1e6],ylim=[1e10,1e19],yliml=[1e6,1e15],bounds=[40e6,2e9],legpos=3);
	
    return imp_mod,wake_mod;
    
