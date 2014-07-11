#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);
from string import *
from DELPHI import *
import pickle as pick;

# script encapsulating a DELPHI calculation with scans
# launched in command line with one argument, the name of the file containing
# all the parameters (pickle format)


if __name__ == "__main__":

    if len(sys.argv)>1: filename=str(sys.argv[1]);
    print filename;
    
    # read parameters
    fileall=open(filename,'r');
    listnames=['Qpscan','nxscan','dampscan','Nbscan','omegasscan','dphasescan','M',
	    'omega0','Q','gamma','eta','a','b','taub','g','Z','freq',
	    'particle','flagnorm','flag_trapz','flagdamperimp','d','freqd',
	    'kmax','kmaxplot','crit','abseps','flagm0']
    for name in listnames[:-1]: exec(name+"=pick.load(fileall)");
    # case of flagm0 is particular (was not there in first versions of pickle files)
    for name in listnames[-1:]:
	try: exec(name+"=pick.load(fileall)");
	except EOFError: exec(name+"=False");
    fileall.close();
    
    # remove initial pickle file
    os.system("rm -vf "+filename);
    
    if flagm0:
	tuneshift_most,tuneshiftnx,tuneshiftm0=eigenmodesDELPHI_converged_scan(Qpscan,nxscan,dampscan,Nbscan,
    	    omegasscan,dphasescan,M,omega0,Q,gamma,eta,a,b,taub,g,Z,freq,
	    particle=particle,flagnorm=flagnorm,flag_trapz=flag_trapz,
	    flagdamperimp=flagdamperimp,d=d,freqd=d,kmax=kmax,kmaxplot=kmaxplot,
	    crit=crit,abseps=abseps,flagm0=flagm0)
    else:
	tuneshift_most,tuneshiftnx=eigenmodesDELPHI_converged_scan(Qpscan,nxscan,dampscan,Nbscan,
    	    omegasscan,dphasescan,M,omega0,Q,gamma,eta,a,b,taub,g,Z,freq,
	    particle=particle,flagnorm=flagnorm,flag_trapz=flag_trapz,
	    flagdamperimp=flagdamperimp,d=d,freqd=d,kmax=kmax,kmaxplot=kmaxplot,
	    crit=crit,abseps=abseps,flagm0=flagm0)

    print tuneshift_most.shape,tuneshiftnx.shape
    print tuneshift_most#,tuneshiftnx
    if flagm0: print tuneshiftm0.shape,tuneshiftm0
    
    tuneshiftnx=[]; # 7/12/2013: I erase it because it is never used... (to save space in multibunch)
    # save results
    fileall=open('out'+filename,'w');
    pick.dump(tuneshift_most,fileall);
    pick.dump(tuneshiftnx,fileall);
    if flagm0: pick.dump(tuneshiftm0,fileall);
    fileall.close();
