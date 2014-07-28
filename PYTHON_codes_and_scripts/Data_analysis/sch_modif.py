#!/usr/bin/python

import sys
import commands
# import local libraries if needed
pymod=commands.getoutput("echo $PYMOD");
if pymod.startswith('local'):
    py_numpy=commands.getoutput("echo $PY_NUMPY");sys.path.insert(1,py_numpy);
    py_matpl=commands.getoutput("echo $PY_MATPL");sys.path.insert(1,py_matpl);

sys.path.append("/afs/cern.ch/eng/sl/lintrack/Python_Classes4MAD/")
import pylab,re, time,dateutil,datetime,random
import numpy as np
from string import split, replace
from optparse import OptionParser
from Timber import parseout
from plot_lib import cmap
from datetime_lib import tt,tb

# interpret data from LHC Schottky monitors
# From R. Calaga

def permute(inputData, outputSoFar):
    for elem in inputData:
        if elem not in outputSoFar:
            outputSoFar.append(elem)
            if len(outputSoFar) == len(inputData):print outputSoFar
            else: permute(inputData, outputSoFar)
            outputSoFar.pop()

def parsse():
    parser = OptionParser()
    parser.add_option("-a", action="store_true",help="compute average",
                      default=False,dest="AVG")
    parser.add_option("-o", "--output",help="Specify output filename",
                      default=None,dest="OUT")
    parser.add_option("-b", "--beam",
                      help="Specify the beam 1 or 2",
                      metavar="BEAM", default=None,dest="BEAM")
    parser.add_option("-p", "--plot",
                      help="Specify 1 for Qx Vs. Qy, 2 for Qx/Qy vs. bunch number",type=int,
                      metavar="PLOT", default=0,dest="PLOT")
    parser.add_option("-t", "--time",
                      help="Specify time gate -t 0 90",nargs=2,
                      type=int,metavar="TIME", default=[15,50],dest="TIME")
    parser.add_option("-f", "--file",
                      help="Specify the file name",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-k", "--bunch",type=int,
                      help="Specify bunch number",action="append", 
                      metavar="BNUM", default=None,dest="BNUM")
    parser.add_option("-c", "--batch",type=int,nargs=3,
                      help="Specify batch to look at and number of bunches in first batch and in all the others: -c3 12 36", 
                      metavar="BATCH", default=None,dest="BATCH")
    parser.add_option("-i", "--first",type=int,
                      help="Specify first bunch to look at", 
                      metavar="FIRST", default=None,dest="FIRST")
    parser.add_option("-l", "--last",type=int,
                      help="Specify last bunch to look at", 
                      metavar="LAST", default=None,dest="LAST")
    parser.add_option("-s", "--sample",
                      help="Sample every nth bunch",type=int,
                      metavar="SAMPLE", default=1,dest="SAMPLE")
    parser.add_option("-x", "--tunex",type=float,nargs=2,
                      help="Specify Qx tune and tolerance at which we accept meas.: -x0.28 0.1", 
                      metavar="QX", default=None,dest="QX")
    parser.add_option("-y", "--tuney",type=float,nargs=2,
                      help="Specify Qy tune and tolerance at which we accept meas.: -y0.31 0.1", 
                      metavar="QY", default=None,dest="QY")
    (opt, args) = parser.parse_args()
    print "Selected File:", opt.FILE
    return opt, args

def plot(a,b,c,be,ce):
    #dt=[dateutil.parser.parse(j) for j in x0];pylab.subplot(2,1,1)
    #pylab.plot_date(pylab.date2num(dt),x1,yerr=x1s,fmt='o',label=str(bnum))
    #pylab.plot_date(pylab.date2num(dt),y1,fmt='o',label=str(bnum))
    #pylab.legend();
    pylab.subplot(2,1,1)
    pylab.errorbar(np.array(a),b,yerr=be,fmt='o');#,label=str(bnum))
    pylab.xlabel("Time from beginning ( "+ tbeg + " ) in sec.")
    pylab.ylabel("Qx");
    #pylab.legend();
    pylab.subplot(2,1,2)
    pylab.errorbar(np.array(a),c,yerr=ce,fmt='o');#,label=str(bnum))
    pylab.xlabel("Time from beginning ( "+ tbeg + " ) in sec.")
    pylab.ylabel("Qy");
    #pylab.legend()

def bplot(b,c,be,ce):
    pylab.subplot(2,1,1)
    pylab.errorbar(bnum,b,yerr=be,fmt='o');
    pylab.xlabel("Bunch 25ns-bucket number");
    pylab.ylabel("Qx");
    pylab.subplot(2,1,2)
    pylab.errorbar(bnum,c,yerr=ce,fmt='o');
    pylab.xlabel("Bunch 25ns-bucket number");
    pylab.ylabel("Qy");

def wpplot(a,b,c,d):
    t=np.arange(2)/2.0;t1=t[::-1];t2=np.ones(2)/2.0
    for m in range(1,6):
        pylab.plot(t,m*t,c='y');pylab.plot(m*t,t,c='y');
        pylab.plot(t,m*t1,c='y'); pylab.plot(m*t,t1,c='y')
        pylab.plot(t2/m,c='y');pylab.plot(t2/m,t,c='y')
        pylab.plot(0.5-t2/m,c='y');pylab.plot(0.5-t2/m,t,c='y')
       # pylab.plot(0.5-t2/m,c='y');pylab.plot(0.5-t2/m,t,c='y')
     
    pylab.errorbar(a,b,xerr=c,yerr=d,\
                   label=str(bnum),marker='o',color=cl);
    pylab.ylim(0.29,0.33);pylab.xlim(0.28,0.32)
    pylab.xlabel("Qx");pylab.ylabel("Qy");

def dct(dt,bm):
    kx='LHC.BQS.'+bm+':BUNCH_SEL';idd={}
    for i,j in enumerate(dt[kx][1]):
        for k,l in enumerate(j):
            if l=='1':
                if k not in idd: idd[k]=[]
                idd[k].append(dt[kx][0][i])
                break
    print "Bunches Available:",np.sort(idd.keys())
    return idd

def tunesAVG(dt,bm,id):
    x0=[];x1=[];y1=[];start=0;x1s=[];y1s=[]
    kh='LHC.BQS.'+bm+':TUNE_H';kv='LHC.BQS.'+bm+':TUNE_V'
    kx='LHC.BQS.'+bm+':BUNCH_SEL';
    tx0=np.array([j for j in dt[kh][0]]);
    ty0=np.array([j for j in dt[kv][0]]);
    tx1=np.array([float(j[0]) for j in dt[kh][1]])
    ty1=np.array([float(j[0]) for j in dt[kv][1]])
    
    print 'To match',len(id),'time stamps...'
    for dy in id:
        s0x=[];s0y=[];s1=[];s2=[];
	j=dt[kx][0].index(dy);#sys.stdout.write('\n%d \n '% j);
	if (j<len(dt[kx][0])-1): timeend=min(opt.TIME[1],tt(dt[kx][0][j+1][:-4])-tt(dy[:-4]));
	else: timeend=min(opt.TIME[1],tt(dy[:-4])-tt(dt[kx][0][j-1][:-4]));
	#sys.stdout.write('%s \n '% timeend);
        sys.stdout.write('\r Matching timestamp --> %s, endtime= %d s'% (dy, timeend))
	sys.stdout.flush();
	for dx in range(start,len(tx0)):
            if opt.TIME[0]<tt(tx0[dx][:-4])-tt(dy[:-4])<timeend:
		if (opt.QX==None) or ((opt.QX[0]-opt.QX[1])<tx1[dx]<(opt.QX[0]+opt.QX[1])):
                	s0x.append(tx0[dx]);s1.append(tx1[dx]);
			#sys.stdout.write('%f \n '% tx1[dx]);
		if (opt.QY==None) or (opt.QY[0]-opt.QY[1]<=ty1[dx]<=opt.QY[0]+opt.QY[1]):
                	s0y.append(ty0[dx]);s2.append(ty1[dx]);
		start=dx;
        
        x0.append(tt(dy[:-4])-tt(tbeg[:-4]));
        x1.append(np.average(s1));x1s.append(np.std(s1));
        y1.append(np.average(s2));y1s.append(np.std(s2));
    if opt.PLOT==1: wpplot(x1,y1,x1s,y1s)
    elif opt.PLOT==2:
        for i in range(0,len(x1)): bplot(x1[i],y1[i],x1s[i],y1s[i]);
    else: plot(x0,x1,y1,x1s,y1s);
    return x0,x1,y1,x1s,y1s
    
def tunes(dt,bm,id):
    x0=[];x1=[];y1=[];start=0;
    kh='LHC.BQS.'+bm+':TUNE_H';kv='LHC.BQS.'+bm+':TUNE_V'
    kx='LHC.BQS.'+bm+':BUNCH_SEL';
    tx0=np.array([j for j in dt[kh][0]]);
    ty0=np.array([j for j in dt[kv][0]]);
    tx1=np.array([float(j[0]) for j in dt[kh][1]])
    ty1=np.array([float(j[0]) for j in dt[kv][1]])
    print 'To match',len(id),'time stamps...'
    for dy in id:
	j=dt[kx][0].index(dy);
	if (j<len(dt[kx][0])-1): timeend=min(opt.TIME[1],tt(dt[kx][0][j+1][:-4])-tt(dy[:-4]));
	else: timeend=min(opt.TIME[1],tt(dy[:-4])-tt(dt[kx][0][j-1][:-4]));
        sys.stdout.write('\r Matching timestamp --> %s, endtime= %d s'% (dy, timeend))
        sys.stdout.flush();
        for dx in range(start,len(tx0)):
            if opt.TIME[0]<tt(tx0[dx][:-4])-tt(dy[:-4])<timeend:
                x0.append(tt(tx0[dx][:-4]));start=dx;
		if (opt.QX==None) or ((opt.QX[0]-opt.QX[1])<tx1[dx]<(opt.QX[0]+opt.QX[1])): x1.append(tx1[dx]);
		else : x1.append(nan)
		if (opt.QY==None) or (opt.QY[0]-opt.QY[1]<=ty1[dx]<=opt.QY[0]+opt.QY[1]): y1.append(ty1[dx]);
		else : y1.append(nan)
    x1s=np.zeros(len(x1));y1s=np.zeros(len(y1))
    if opt.PLOT==1: wpplot(x1,y1,x1s,y1s)
    elif opt.PLOT==2:
        for i in range(0,len(x1)): bplot(x1[i],y1[i],x1s[i],y1s[i]);
    else: plot(x0,x1,y1,x1s,y1s)
    return x0,x1,y1,x1s,y1s

def wrt(x0,x1,y1,x1s,y1s):
    file=open(opt.OUT+'.bunch'+str(bnum),'w');
    print >> file, "@  FILE %s ",opt.OUT
    print >> file, "*  DATE   TIME    Q1    Q2     sQ1    sQ2"
    print >> file, "$   %s     %s     %le   %le    %le    %le"
    for j in range(len(x0)):
        print >> file, tb(x0[j]), x1[j],y1[j],x1s[j],y1s[j]
    file.close()

def cmap1():
    maps=np.array((np.linspace(0,1,len(opt.BNUM)/3),
                   np.linspace(0,1,len(opt.BNUM)/3),
                   np.linspace(1,0,len(opt.BNUM)/3)))
    return np.transpose(maps)
                  
if __name__ == "__main__":
    opt,args=parsse(); data=parseout(opt.FILE);

    if opt.BEAM=="1": beam="B1"
    elif opt.BEAM=="2": beam="B2"
    else: print "specify beam 1 or 2"; sys.exit()
    kx='LHC.BQS.'+beam+':BUNCH_SEL';
    tbeg=data[kx][0][0];
    id=dct(data,beam);iter=0;
    if opt.BNUM==None: 
	opt.BNUM=np.sort(id.keys())[::opt.SAMPLE]
	if opt.BATCH!=None:
		ba1=opt.BATCH[1];
		ba2=opt.BATCH[2];
		ban=opt.BATCH[0];
		if ban==1: opt.BNUM=opt.BNUM[1:ba1+1];
		else: opt.BNUM=opt.BNUM[1+ba1+ba2*(ban-2):ba1+ba2*(ban-1)+1];
	else:
		if opt.FIRST==None: first=0;
		else: 
			first=0;
			while (opt.BNUM[first] < opt.FIRST): first=first+1;
		if opt.LAST==None: last=len(opt.BNUM)-1;
		else: 
			last=len(opt.BNUM)-1;
			while (opt.BNUM[last] > opt.LAST): last=last-1;
		opt.BNUM=opt.BNUM[first:last+1];
    maps=cmap(len(opt.BNUM))
    
    print "Default time range:",opt.TIME[0],'s -->',opt.TIME[1],'s'
    print "Matching",len(opt.BNUM),"bunches out of",len(id.keys())    
    for bnum in opt.BNUM:
        cl=maps[iter];iter+=1
        if bnum in id.keys():
            print "\n Matching bunch",bnum
            if opt.AVG:
                a1,b1,c1,b1s,c1s=tunesAVG(data,beam,id[bnum])
            else:
                a1,b1,c1,b1s,c1s=tunes(data,beam,id[bnum])
            if opt.OUT!=None: wrt(a1,b1,c1,b1s,c1s)
        else: "Bunch", bnum, 'not in TIMBER data'
    pylab.show();sys.exit()

