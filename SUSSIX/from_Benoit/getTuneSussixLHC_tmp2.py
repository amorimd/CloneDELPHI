#!/afs/cern.ch/eng/sl/lintrack/Python-2.5/bin//python2.5

from string import *
from numpy import *
from numpy.fft import *
from sys import *
import os
import re
from SussixBS import *




###########################
def gettune(filename):
###########################
    # Definition of the betatron tune (in m)
    betax=42.0941 
    betay=42.0137 
    x=[]
    y=[]
    px=[]
    py=[]
    for line in open(filename ):
        x.append(float(split(line)[1]))
        px.append(float(split(line)[2]))
	y.append(float(split(line)[3]))
	py.append(float(split(line)[4]))
    x=array(x)
    px=betax*array(px)
    y=array(y)
    py=betay*array(py)
    tunex=argmax(abs(fft.fft(x))[1:])*1.0/len(x)
    tuney=argmax(abs(fft.fft(y))[1:])*1.0/len(y)
    xfft=abs(fft.fft(x))[1:]*1.0/len(x)
    yfft=abs(fft.fft(y))[1:]*1.0/len(y)
    print 'fft tunes', tunex, tuney, len(x), len(y), len(xfft), len(yfft),'\n'
    sussix_inp(ir=0, turns=8192, tunex=0.13, tuney=0.18, istun=0.02, idam=2, narm=300)
    a=sussixBS(x,y,px,py)
    return a,xfft,yfft



#print "your file is: ",argv[1]

suss=gettune(argv[1])

#Frequency and amplitude of lines:

outf=open(argv[1]+'.sussix_x','w')

dicx=dict(transpose([suss[0].ox,suss[0].ax]))

for i in sort(suss[0].ox):
    outf.write( str(i)+' '+str(dicx[i])+'\n')

outf.close()

dicy=dict(transpose([suss[0].oy,suss[0].ay]))

outf=open(argv[1]+'.sussix_y','w')

for i in sort(suss[0].oy):
    outf.write( str(i)+' '+str(dicy[i])+'\n')

outf.close()

outf=open(argv[1]+'.fft','w')

for i in range(len(suss[1])):
    outf.write( str(i*1.0/len(suss[1]))+'  '+str(suss[1][i])+' '+str(suss[2][i])+'\n')

outf.close()
#Tunes are here
#print a.tunexy
