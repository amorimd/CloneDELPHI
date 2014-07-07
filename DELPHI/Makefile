# Makefile for DELPHI -- LINUX GCC

.SUFFIXES: .cc

CC = g++

OBJS = DELPHI.o

LOPTS = 

# -Wall -> warning all

all: DELPHI lib


DELPHI.o: DELPHI.cc
	$(CC) $(LOPTS) -fPIC -o DELPHI.o -c DELPHI.cc

DELPHI: $(OBJS)
	$(CC) $(LOPTS) -o DELPHI.x $(OBJS) -O2 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas

clean:
	rm -f $(OBJS) libDELPHI.so DELPHI.x
	
lib: $(OBJS)
	$(CC) $(LOPTS) -shared -o libDELPHI.so $(OBJS) -O2 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas
