#!/usr/bin/python2.6

import sys
sys.path.append("Sussix")
from optparse import OptionParser
from string import *
from numpy import *
from numpy.fft import *
from sys import *
import os
import re
from SussixNM import *


def parsse():
    parser = OptionParser()
    parser.add_option("-q", "--tunes",nargs=3,
                      help="Specify the tunes: e.g. -q 0.13 0.18 0.07",
                      metavar="TUNES", default=[0.31,0.32,0.07],dest="TUNES")
    parser.add_option("-k", "--bunch",type=int,
                      help="Specify bunch number",action="append",
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-n", "--nbunch",type=int,
                      help="Specify total number of bunches",
                      metavar="NB", default=36,dest="NB")
    parser.add_option("-t", "--turns",type=int,
                      help="Specify number of turns to take into account (from the end of data)",
                      metavar="TURNS", default=8192,dest="TURNS")
    parser.add_option("-i", "--istun",
                      help="Specify istun option (default 0.02)",
                      metavar="ISTUN", default=0.02,dest="ISTUN")
    parser.add_option("-r", "--ir",type=int,
                      help="Specify ir option (default=0:consider both x and x'; 1:consider only x)",
                      metavar="IR", default=0,dest="IR")
    parser.add_option("-b", "--beta",nargs=2,
    		      help="Specify the average beta functions: e.g. -b 42 42",
                      metavar="BETA", default=[66,72],dest="BETA")
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default=None,dest="FILE")
    (opt, args) = parser.parse_args()
    return opt, args


###########################
def gettune(filename,bunches,nbunch,betax,betay,Qx,Qy,Qz,istun=0.02,ir=1,nturns=8192):
###########################
    x=[]
    y=[]
    px=[]
    py=[]
    for line in open(filename):
        x.append(float(split(line)[1]))
        px.append(float(split(line)[2]))
	y.append(float(split(line)[3]))
	py.append(float(split(line)[4]))
    x=array(x)
    px=betax*array(px)
    y=array(y)
    py=betay*array(py)
    xfft=[];yfft=[];a=[];

    for j,bnum in enumerate(bunches):
    
	x1=x[-nbunch*nturns+bnum-1::nbunch]
    	px1=px[-nbunch*nturns+bnum-1::nbunch]
    	y1=y[-nbunch*nturns+bnum-1::nbunch]
    	py1=py[-nbunch*nturns+bnum-1::nbunch]
	tunex=argmax(abs(fft.fft(x1))[1:])*1.0/len(x1)
	tuney=argmax(abs(fft.fft(y1))[1:])*1.0/len(y1)
	xfft.append(abs(fft.fft(x1))[1:]*1.0/len(x1))
	yfft.append(abs(fft.fft(y1))[1:]*1.0/len(y1))
	sussix_inp(ir=ir, turns=nturns, tunex=Qx, tuney=Qy, tunez=Qz, istun=istun, idam=2, narm=300)
	a.append(sussixBS(x1,y1,px1,py1))
	print 'fft tunes bunch no', bnum, ' : ', tunex, tuney
	print 'sussix tunes bunch no', bnum, ' : ', a[j].tunexy[0],a[j].tunexy[1]
	print 'sussix tune shifts bunch no', bnum, ' : ', a[j].tunexy[0]-Qx,a[j].tunexy[1]-Qy
        #print a[j].amplitude
	#for i in range(len(a[j].ox)):
	#    print a[j].ox[i], a[j].ax[i],a[j].oy[i], a[j].ay[i]

    
    return a,xfft,yfft



#print "your file is: ",argv[1]

opt,args=parsse();
if opt.BNUM==None: opt.BNUM=range(1,opt.NB+1);

susstab,xfft,yfft=gettune(opt.FILE,opt.BNUM,opt.NB,betax=float(opt.BETA[0]),betay=float(opt.BETA[1]),
	Qx=float(opt.TUNES[0]),Qy=float(opt.TUNES[1]),Qz=float(opt.TUNES[2]),istun=float(opt.ISTUN),
	ir=opt.IR,nturns=opt.TURNS)


for j,suss in enumerate(susstab):

    #Frequency and amplitude of lines:

    outf=open(opt.FILE+'.sussix_x_bunch'+str(opt.BNUM[j]),'w')

    dicx=dict(transpose([suss.ox,suss.ax]))

    for i in sort(suss.ox):
	outf.write( str(i)+' '+str(dicx[i])+'\n')

    outf.close()

    dicy=dict(transpose([suss.oy,suss.ay]))

    outf=open(opt.FILE+'.sussix_y_bunch'+str(opt.BNUM[j]),'w')

    for i in sort(suss.oy):
	outf.write( str(i)+' '+str(dicy[i])+'\n')

    outf.close()

    outf=open(opt.FILE+'.fft_bunch'+str(opt.BNUM[j]),'w')

    for i in range(len(xfft[j])):
	outf.write( str(i*1.0/len(xfft[j]))+'  '+str(xfft[j][i])+' '+str(yfft[j][i])+'\n')

    outf.close()
