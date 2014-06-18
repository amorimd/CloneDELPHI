import os
import gzip
import time

import numpy as np
from numpy import *
import matplotlib.pyplot as pl

# Author: R. de Maria

#from workpkg import *

namemaps={
  'ADTBposHorQ7LB1':'q7hb1',
  'ADTBposHorQ9LB1':'q9hb1',
  'ADTBposVerQ7RB1':'q7vb1',
  'ADTBposVerQ7RB1':'q9vb1',
  'ADTBposHorQ7LB2':'q7vb2',
  'ADTBposHorQ9LB2':'q9vb2',
  'ADTBposVerQ7RB2':'q7hb2',
  'ADTBposVerQ7RB2':'q9hb2',
}

def c2db(c):
  """complex signal to decibel"""
  return 20.*np.log10(abs(c))

def t2f(t):
  n=len(t)/2+1
  fs=1/(t[1]-t[0])
  return np.arange(0.,n)*fs/len(t)

def movavg2(x,n=100):
  """moving average with rect"""
  return np.convolve(np.ones(float(n))/n,x,mode='same')





class BPM(object):
  def __init__(self,pos,bunchNumbers=None,acqStamp=None,deviceName=None,sigma=None,calib=None,**nargs):
    self.pos=pos
    self.turns,self.bunches=pos.shape
    self.bunchNumbers=bunchNumbers
    self.acqStamp=acqStamp
    self.calib=calib
    self.deviceName=deviceName
    self.sigma=None
    self.__dict__.update(nargs)
  def plotTune(self,bunch,smooth=None):
    v=self.pos[bunch]
    vf=c2db(np.fft.rfft(v))
    f=t2f(np.arange(len(v)))
    if smooth is not None:
      wdw=smooth
      vf=movavg2(vf,wdw)
    pl.plot(f,vf,label="b%d"%self.bunchNumbers[bunch],alpha=.4)
  def plotTuneAll(self,smooth=None):
    for bunch in range(len(self.bunchNumbers)):
      self.plotTune(bunch,smooth=smooth)

  def plotPos(self,bunch):
    pl.plot(self.pos[bunch])
  def plotTitle(self,bunch=None,timestamp=None,bunchNumbers=None):
    title=self.deviceName
    if bunch:
      title+=' bunch %s'%(self.bunchNumbers[bunch])
    if bunchNumbers:
      title+='\n%s'%(','.join(map(str,self.bunchNumbers)))
    if timestamp and self.acqStamp is not None:
      timestamp=time.asctime(time.localtime(self.acqStamp))
      title+='\n%s'%timestamp
    pl.title(title)


class MultiTurnData(object): # for ADT "Verena kind of data" (i.e. from the multiturn application)
  namemaps={
    'ADTBposHorQ7LB1':'q7hb1',
    'ADTBposHorQ9LB1':'q9hb1',
    'ADTBposVerQ7RB1':'q7vb1',
    'ADTBposVerQ7RB1':'q9vb1',
    'ADTBposHorQ7LB2':'q7vb1',
    'ADTBposHorQ9LB2':'q9vb1',
    'ADTBposVerQ7RB2':'q7hb1',
    'ADTBposVerQ7RB2':'q9hb1',
  }

  def __init__(self,fn):
    self.filename=fn
    if hasattr(fn,'readline'):
      fh=fn
    else:
      if fn.endswith('.gz'):
        fh=gzip.open(fn)
      else:
        fh=file(fn)
    self.beam=int(fh.readline().split()[-1])
    self.bunchNumbers=map(int,fh.readline().split()[1:])
    self.monitors=fh.readline().split()[1:]
    self.turns=int(fh.readline().split()[-1])
    self.bunches=262144/self.turns
    self.timestamp=float(fh.readline().split()[-1][4:])/1e9
    fh.readline()
    data={}
    self.deviceNames=[]
    for i in range(len(self.monitors)*2):
      line=fh.readline()
      if line[7]==line[16]:
        deviceName=line[:15]
        arr=np.fromstring(line[17:],dtype=float,sep=" ")
        arr=arr.reshape(self.bunches,self.turns)
        np.array(arr,dtype=int)
        bpm=BPM(arr,acqStamp=self.timestamp,calib=1.0,
                deviceName=deviceName,bunchNumbers=self.bunchNumbers)
        setattr(self,deviceName,bpm)
        self.deviceNames.append(deviceName)
  def plotTuneAll(self,xlim=[0.25,0.33],ylim=[60,120]):
    iplot=len(self.deviceNames)/2*100+20
    for deviceName in self.deviceNames:
      iplot+=1
      pl.subplot(iplot)
      bpm=getattr(self,deviceName)
      bpm.plotTuneAll()
      bpm.plotTitle()
      pl.xlim(*xlim)
      pl.ylim(*ylim)
      #pl.legend()
    tlt='%s\n%s'%(self.getBaseName(),','.join(map(str,self.bunchNumbers)))
    pl.suptitle(tlt)
  def getBaseName(self):
    return os.path.basename(self.filename).split('.')[0]
  def getTimeStamp(self,fmt="'%Y-%m-%d_%H-%M-%S"):
    pass



class LHCDamperBPMData(object): # for ADT "Riccardo kind of data" (i.e. from his scripts)
  def log(self,msg):
    print "%s:"% self.filename.split('/')[-1],msg
  def __init__(self,fn):
    self.filename=fn
    if hasattr(fn,'readline'):
      fh=fn
    else:
      if fn.endswith('.gz'):
        fh=gzip.open(fn)
      else:
        fh=file(fn)
    devices=int(fh.readline().split()[-1])
    self.deviceNames=[]
    if devices==0:
      log("Warning: wrong device name")
      devices=16
    for  i in range(devices):
      deviceName=fh.readline().split()[-1]
      channelName=fh.readline().split()[-1]
      timestamp=float(fh.readline().split()[-1])/1e3
      turns=int(fh.readline().split()[-1])
      bunchNumbers=map(int,fh.readline().split()[1:])
      bunchNumbers=bunchNumbers
      fh.readline()
      arr=np.fromstring(fh.readline(),dtype=int,sep=" ")
      if bunchNumbers==[0,0,0,0,0,0,0,0]:
        self.mode=0
        self.bunches=3564
        self.turns=len(arr)/3564
        narr=arr[:self.turns*3564].reshape(self.turns,3564)
        if channelName=='OBS_RPOS_BUNCH':
          bpm=BPM(narr,acqStamp=timestamp,calib=1.0,
                   deviceName=deviceName,posflat=arr)
          setattr(self,deviceName,bpm)
          self.deviceNames.append(deviceName)
        elif channelName =='OBS_SIGMA_MAG':
          bpm=getattr(self,deviceName)
          bpm.sigma=narr
          bpm.sigmaflat=arr
      else:
        self.turns=turns
        self.mode=[1,3,5,7][int(log2(262144/turns))]
        self.bunches=len(bunchNumbers)
        if channelName=='OBS_RPOS_BUNCH':
          narr=arr.reshape(self.bunches,turns).T
          bpm=BPM(narr,acqStamp=timestamp,calib=1.0,
                   deviceName=deviceName,bunchNumbers=bunchNumbers)
          setattr(self,deviceName,bpm)
          self.deviceNames.append(deviceName)
  def checkData(self):
    good=True
    for deviceName in self.deviceNames:
      bpm=getattr(self,deviceName)
      if bpm.bunchNumbers:
        empty_slot=[]
        for i in range(bpm.bunches):
          if bpm.pos[:,i].min()==0 and bpm.pos[:,i].max()==0:
            empty_slot.append(bpm.bunchNumbers[i])
            good=False
        if empty_slot:
          bl=','.join(map(str,empty_slot))
          self.log('Error: Bunch %s empty'%(bl))
          break
      else:
        if bpm.pos.min()==0 and bpm.pos.max()==0:
          self.log('Error: Data empty')
          good=False
    if good:
      self.log('Checked')
    return good
  def plotTuneAll(self,xlim=[0.25,0.33],ylim=[60,120]):
    iplot=len(self.deviceNames)/2*100+20
    for deviceName in self.deviceNames:
      iplot+=1
      pl.subplot(iplot)
      bpm=getattr(self,deviceName)
      bpm.plotTuneAll()
      bpm.plotTitle(bunchNumbers=True)
      pl.xlim(*xlim)
      pl.ylim(*ylim)
      #pl.legend()
    tlt='%s'%(self.getBaseName())
    pl.suptitle(tlt)
  def plotPosNmodeAll(self,step=100):
    iplot=len(self.deviceNames)/2*100+20
    for deviceName in self.deviceNames:
      iplot+=1
      pl.subplot(iplot)
      bpm=getattr(self,deviceName)
      idx=find_12_bunches(bpm.posflat)
      delta=bpm.posflat[idx:]
      turns=len(delta)/3564
      delta=delta[:turns*3564].reshape(turns,3564)
      sigma=bpm.sigma
      delta=diff(delta,axis=0)
      for i in range(turns-1):
        plot(delta[i]+step*i)
      bpm.plotTitle()
    tlt='%s'%(self.getBaseName())
    pl.suptitle('Nmode position '+tlt)
  def getBaseName(self):
    return os.path.basename(self.filename).split('.')[0]

def find_12_bunches(delta):
  b12=zeros(24,dtype=bool)
  zz=zeros(48,dtype=bool)
  b12[::2]=True
  b12=r_[zz,b12,zz]
  for i in range(10000):
    res=all((r_[delta][i:i+len(b12)]!=0)==b12)
    if res:
      break
  return i






def get_slots_number_fromdb(t1,t2,beam=1,last=10):
  data=cernlogdb.dbget('LHC.BQM.B%d:FILLED_BUCKETS'%beam,t1,t2,conf='mdb.conf')
  data=cernlogdb.combine_data(data,vtype=int)
  t,v=data[0]
  v=(v-1)/10
  for tt,vv in zip(t,v):
    print rdmdate.dumpdate(tt/1000)
    nbunch=where(vv==-1)[0][0]
    if last is not None:
      vv=vv[nbunch-last:nbunch]
    else:
      vv=vv[nbunch-last:nbunch]
    print ','.join(array(vv,dtype='str'))
    print








