#!/usr/bin/python

# From Xavier BUFFAT

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);
    py_scipy=commands.getoutput("echo $PY_SCIPY");sys.path.insert(1,py_scipy);

sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
from Timber import parseout,extractnew,toseconds
import numpy as np
import scipy as sp
import math
import matplotlib.pyplot as plt
import csv
from optparse import OptionParser
from io_lib import list_files


def parsse():
    parser = OptionParser()
    parser.add_option("-s", "--slotminmax",type=int,
                      help="Specify the minimum and maximum slots to plot (default=0 & 3564)",
                      metavar="SLOT", default=None,dest="SLOT")
    parser.add_option("-f", "--file",action="append",
                      help="Specify the _fits.dat names of the files to plot. Several files are possible - either with regular expression (using *, [], etc. and BETWEEN QUOTES \"\") or with several -f options.",
                      metavar="FILE", default=None,dest="FILE")
    (opt, args) = parser.parse_args()
    return opt, args


def gateDelayToSlotNumber(gateDelay,beamNumber):
    if beamNumber==1:
        return int(math.floor(0.5*gateDelay-356.5));
    else:
        return int(math.floor(0.5*gateDelay-1249));

#use data[emittance variable name][slot number][0=time,1=varValue][valueNumber]
def getData(fileName):
    emitVars = ["LHC.BSRTS.5R4.B1:BEAM_NORM_EMIT_H","LHC.BSRTS.5R4.B1:BEAM_NORM_EMIT_V","LHC.BSRTS.5L4.B2:BEAM_NORM_EMIT_H","LHC.BSRTS.5L4.B2:BEAM_NORM_EMIT_V"];
    gateVars = ["LHC.BSRCTL.B1:GATE_DELAY","LHC.BSRCTL.B2:GATE_DELAY"]
    data = {};
    timberData = parseout(fileName);
    for i in np.arange(len(emitVars)):
        print emitVars[i];
        if "B1" in emitVars[i]:
            beamNumber = 1;
        else :
            beamNumber=2;
        if(emitVars[i] in timberData and gateVars[beamNumber-1] in timberData):
            if not emitVars[i] in data:
                data[emitVars[i]] = {};
            tGate,dGate = timberData[gateVars[beamNumber-1]];
            tEmit,dEmit = timberData[emitVars[i]]
            for j in np.arange(len(tGate)-1):
                tGated,dGated = extractnew(tEmit,dEmit,tGate[j],tGate[j+1]);
                slot = gateDelayToSlotNumber(dGate[j][0],beamNumber);
                if not slot in data[emitVars[i]]:
                    t,v = [],[];
                    data[emitVars[i]][slot]=(t,v);
                if len(tGated)>0:
                    avg = 0;
                    for k in np.arange(len(tGated)):
                        avg+=dGated[k][0];
                    avg=avg/len(tGated);
                    data[emitVars[i]][slot][0].append(tGated[0]);
                    data[emitVars[i]][slot][1].append(avg);
#                for k in np.arange(len(tGated)):
#                    data[emitVars[i]][slot][0].append(tGated[k]);
#                    data[emitVars[i]][slot][1].append(dGated[k][0]);
    return data;

def toCustomizeFile(data,outputFolder):
    for var in data:
        for slot in data[var]:
                csvWriter = csv.writer(open(outputFolder+"/"+var+"_slot="+str(slot),'wb'));
                t,v = data[var][slot];
                csvWriter.writerows(zip(t,v));

def fromCustomizeFile(inputFolder):
    #TODO
    pass;


def makePlots(data,slotmin=0,slotmax=3564,ifig=0):
    for var in data.keys():
        for slot in data[var].keys():
            if "B1" in var:
                if "_H" in var:
                    plt.figure(ifig+1);
                elif "_V" in var : 
                    plt.figure(ifig+2);
                else:
                    plt.figure(ifig+0);
            elif "B2" in var:
                if "_H" in var:
                    plt.figure(ifig+3);
                elif "_V" in var:
                    plt.figure(ifig+4);
                else:
                    plt.figure(ifig+0);
            else:
                plt.figure(ifig+0);
            if (slot<=slotmax)and(slot>=slotmin):
	    	plt.plot([toseconds(value)-toseconds(data[var][slot][0][0]) for value in data[var][slot][0]],data[var][slot][1]);
		plt.title(var);
		plt.xlabel('Time (seconds) since '+data[var][slot][0][0]);
		plt.ylabel('Normalized RMS emittance');


def makeGrowthPlot(data):
    growthRates = {};
    t,d = [],[];
    growthRates["HB1"] = (t,d);
    t,d = [],[];
    growthRates["VB1"] = (t,d);
    t,d = [],[];
    growthRates["HB2"] = (t,d);
    t,d = [],[];
    growthRates["VB2"] = (t,d);
    for var in data.keys():
        for slot in data[var].keys():
            if len(data[var][slot][0])==0:
                print var+", "+str(slot);
            else:
                (a,b)=sp.polyfit([toseconds(value) for value in data[var][slot][0]], data[var][slot][1], 1);
                if "B1" in var:
                    if "_H" in var:
                        key = "HB1";
                    elif "_V" in var : 
                        key = "VB1";
                elif "B2" in var:
                    if "_H" in var:
                        key = "HB2";
                    elif "_V" in var:
                        key = "VB2";
                growthRates[key][0].append(slot);
                growthRates[key][1].append(a);
    plt.figure(0);
    plt.bar([value-0.4 for value in growthRates["HB1"][0]],growthRates["HB1"][1],0.2,color='b',edgecolor='b');
    plt.bar([value-0.2 for value in growthRates["VB1"][0]],growthRates["VB1"][1],0.2,color='g',edgecolor='g');
    plt.bar([value for value in growthRates["HB2"][0]],growthRates["HB2"][1],0.2,color='r',edgecolor='r');
    plt.bar([value+0.2 for value in growthRates["VB2"][0]],growthRates["VB2"][1],0.2,color='y',edgecolor='y');
    plt.figure(1);
    plt.subplot(221);
    plt.bar([value-0.4 for value in growthRates["HB1"][0]],growthRates["HB1"][1],0.2,color='b',edgecolor='b');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance growth rate HB1 [Hz]")
    plt.subplot(222);
    plt.bar([value-0.2 for value in growthRates["VB1"][0]],growthRates["VB1"][1],0.2,color='g',edgecolor='g');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance growth rate VB1 [Hz]")
    plt.subplot(223);
    plt.bar([value for value in growthRates["HB2"][0]],growthRates["HB2"][1],0.2,color='r',edgecolor='r');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance growth rate HB2 [Hz]")
    plt.subplot(224);
    plt.bar([value+0.2 for value in growthRates["VB2"][0]],growthRates["VB2"][1],0.2,color='y',edgecolor='y');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance growth rate VB2 [Hz]")
    
def makeAvgPlot(data):
    growthRates = {};
    t,d = [],[];
    growthRates["HB1"] = (t,d);
    t,d = [],[];
    growthRates["VB1"] = (t,d);
    t,d = [],[];
    growthRates["HB2"] = (t,d);
    t,d = [],[];
    growthRates["VB2"] = (t,d);
    for var in data.keys():
        for slot in data[var].keys():
            if len(data[var][slot][0])==0:
                print var+", "+str(slot);
            else:
                a=np.average(data[var][slot][1]);
                if "B1" in var:
                    if "_H" in var:
                        key = "HB1";
                    elif "_V" in var : 
                        key = "VB1";
                elif "B2" in var:
                    if "_H" in var:
                        key = "HB2";
                    elif "_V" in var:
                        key = "VB2";
                growthRates[key][0].append(slot);
                growthRates[key][1].append(a);
    plt.figure(0);
    plt.bar([value-0.4 for value in growthRates["HB1"][0]],growthRates["HB1"][1],0.2,color='b',edgecolor='b');
    plt.bar([value-0.2 for value in growthRates["VB1"][0]],growthRates["VB1"][1],0.2,color='g',edgecolor='g');
    plt.bar([value for value in growthRates["HB2"][0]],growthRates["HB2"][1],0.2,color='r',edgecolor='r');
    plt.bar([value+0.2 for value in growthRates["VB2"][0]],growthRates["VB2"][1],0.2,color='y',edgecolor='y');
    plt.figure(1);
    plt.subplot(221);
    plt.bar([value-0.4 for value in growthRates["HB1"][0]],growthRates["HB1"][1],0.2,color='b',edgecolor='b');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance HB1 [\mu m]")
    plt.subplot(222);
    plt.bar([value-0.2 for value in growthRates["VB1"][0]],growthRates["VB1"][1],0.2,color='g',edgecolor='g');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance VB1 [\mu m]")
    plt.subplot(223);
    plt.bar([value for value in growthRates["HB2"][0]],growthRates["HB2"][1],0.2,color='r',edgecolor='r');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance HB2 \mu m]")
    plt.subplot(224);
    plt.bar([value+0.2 for value in growthRates["VB2"][0]],growthRates["VB2"][1],0.2,color='y',edgecolor='y');
    plt.xlabel("Slot number");
    plt.ylabel("Emittance VB2 [\mu m]")


if __name__ == "__main__":

    opt,args=parsse();
    
    if (opt.SLOT!=None):
        bmin=opt.SLOT[0];
        bmax=opt.SLOT[1];
    else:
        bmin=0;bmax=3564;

    # create list of filenames to analyse
    listname=list_files(opt.FILE);

    for filename in listname:
    
    	data = getData(filename);
	#toCustomizeFile(data,"/home/xbuffat/Data/BBMD_LR_2/BSRT/bbb");
	makePlots(data,slotmin=bmin,slotmax=bmax,ifig=2);
	#makeGrowthPlot(data);
    	makeAvgPlot(data);
    
    plt.show();
                
    
    
