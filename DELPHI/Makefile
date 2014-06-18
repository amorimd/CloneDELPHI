# Makefile for DELPHI -- LINUX GCC

.SUFFIXES: .cc

CC = g++

OBJS = DELPHI.o

COPTS = -O2 -fPIC -I/usr/local/include/gsl -llapack -lblas -L/usr/lib
LOPTS = 

# -Wall -> warning all

all: DELPHI lib


DELPHI.o: DELPHI.cc
	$(CC) $(LOPTS) -fPIC -o DELPHI.o -c DELPHI.cc

DELPHI: $(OBJS)
	$(CC) $(LOPTS) -o DELPHI.x $(OBJS) -O3 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas -L/usr/lib -L/usr/local/lib/

clean:
	rm -f $(OBJS) libDELPHI.so DELPHI.x
	
lib: $(OBJS)
	$(CC) $(LOPTS) -shared -Wl -o libDELPHI.so $(OBJS) -O3 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas -L/usr/lib -L/usr/local/lib/
