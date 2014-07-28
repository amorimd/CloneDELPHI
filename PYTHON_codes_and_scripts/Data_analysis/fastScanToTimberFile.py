#!/usr/bin/python

# From Xavier BUFFAT

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import os,string,csv
from Timber import toseconds
import matplotlib.pyplot as plt
import matplotlib as mpl

def parseName(fileName):
    tokens = fileName.split('_');
    return tokens[0],string.replace(tokens[1],'.','-')+' '+tokens[2]+'.000';

def getVarName(beam,plane,comment=''):
    return 'BSRT_'+comment+'_'+plane+'B'+str(beam);

def loadDatafromDirectory(sourceDir):
    fileNames = os.listdir(sourceDir);
    fileNames.sort();
    data = {};
    for fileName in fileNames:
        if 'fits.dat' in fileName:
            var,time = parseName(fileName);
            beam = 1;
            if 'B2' in var:
                beam = 2;
            myFile = open(sourceDir+fileName,'r');
            var0 = None;
            for line in myFile.readlines():
                if var0 == None:
                    if 'RECORD_START' in line:
                        if "Horizontal Size" in line:
                            var0 = getVarName(beam,'H',comment='SIZE');
                        elif "Vertical Size" in line:
                            var0 = getVarName(beam,'V',comment='SIZE');
                        elif 'Horizontal Emittance' in line:
                            var0 = getVarName(beam,'H',comment='EMIT');
                        elif 'Vertical Emittance' in line:
                            var0 = getVarName(beam,'V',comment='EMIT');
                        else:
                            print 'ERROR parsing file '+fileName+' at line : '+line;
                            return None;
                    else:
                        print 'ERROR parsing file '+file+', RECORD_START expected';
                        return None;
                else:
                    if 'RECORD_END' in line:
                        var0 = None;
                        continue;
                    tokens = line.split(',');
                    try:
                        slot = int(float(tokens[0]));
                        value = float(tokens[1]);#TODO vectornumeric
                    except ValueError:
                        print 'ERROR parsing file '+file+': '+line+', is not parsable';
                    if var0 not in data.keys():
                        data[var0] = {};
                    if slot not in data[var0].keys():
                        data[var0][slot] = ([],[]);
                    data[var0][slot][0].append(time);
                    data[var0][slot][1].append(value);
    return data;
                        
def plotData(data):
    t0 = None;
    for var in data.keys():
        if 'EMIT' in var:
            if 'H' in var:
                fig = plt.figure(0);
            else:
                fig = plt.figure(1);
        elif 'SIZE' in var:
            if 'H' in var:
                fig = plt.figure(2);
            else:
                fig = plt.figure(3);
        else:
            fig = plt.figure(4);
        if 'B1' in var:
            fig.add_subplot(211);
        else:
            fig.add_subplot(212);
        print var
        for slot in data[var].keys():
            if t0 == None:
                t0 = data[var][slot][0][0];
                t0s = toseconds(t0);
            plt.plot([toseconds(data[var][slot][0][k])-t0s for k in range(len(data[var][slot][0]))],data[var][slot][1],label = str(slot));
    fig=plt.figure(0);
    fig.add_subplot(211);
    plt.ylabel('Norm. emittance HB1');
    fig.add_subplot(212);
    plt.xlabel('Time [s since '+t0+']');plt.ylabel('Norm. emittance HB2');
    fig=plt.figure(1);
    fig.add_subplot(211);
    plt.ylabel('Norm. emittance VB1');
    fig.add_subplot(212);
    plt.xlabel('Time [s since '+t0+']');plt.ylabel('Norm. emittance VB2');
    fig=plt.figure(2);
    fig.add_subplot(211);
    plt.ylabel('Beam size HB1');
    fig.add_subplot(212);
    plt.xlabel('Time [s since '+t0+']');plt.ylabel('Norm. Beam size HB2');
    fig=plt.figure(3);
    fig.add_subplot(211);
    plt.ylabel('Beam size VB1');
    fig.add_subplot(212);
    plt.xlabel('Time [s since '+t0+']');plt.ylabel('Beam size VB2');
    


def dataToTimberFile(targetDir):
    for var in data.keys():
        for slot in data[var].keys():
                myFile = open(targetDir+var+"_"+str(slot),'wb');
                myFile.write('VARIABLE: '+var+'_'+str(slot)+'\n\nTimestamp (LOCAL_TIME),Value\n');
                csvWriter = csv.writer(myFile);
                t,v = data[var][slot];
                csvWriter.writerows(zip(t,v));

if __name__ == '__main__':
    mpl.rcParams['font.size'] = 30.0;
    sourceDir = '/home/nmounet/Documents/Mesures_LHC/MD_block2_2012/MD_octupole_instability_threshold_19062012/BSRT_fastscan/';
    targetDir = '/home/nmounet/Documents/Mesures_LHC/MD_block2_2012/MD_octupole_instability_threshold_19062012/BSRT_fastscan/Extract/';
    data = loadDatafromDirectory(sourceDir);
    plotData(data);
    dataToTimberFile(targetDir);
    plt.show();
