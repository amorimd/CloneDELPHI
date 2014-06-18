#!/usr/bin/python2.6
#import sys
from string import split, replace
from optparse import OptionParser


def parsse():
    parser = OptionParser()
    parser.add_option("-f", "--file",
                      help="Specify the file name (with extenstion)",
                      metavar="FILE", default=None,dest="FILE")
    parser.add_option("-s", "--string",
                      help="Specify the string to look for",
                      metavar="STR", default=None,dest="STR")
    (opt, args) = parser.parse_args()
    #print "Selected Files:", opt.FILE
    return opt, args


def read_cfg(filename,string):
    
    fh=open(filename,"r")
    for l in fh.readlines():
    	if l.startswith(string):
	    ll=l.strip().split();
	    value=ll[1]

    fh.close()
    
    return value
                  

if __name__ == "__main__":

    opt,args=parsse();
    
    print read_cfg(opt.FILE,opt.STR)
