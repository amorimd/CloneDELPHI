# Makefile for HFMM beam-beam package -- LINUX G77

.SUFFIXES: .cc

CC = g++ -c

OBJS = laclare.o

COPTS = -O2 -I/usr/local/include/gsl -llapack -lblas -L/usr/lib
LOPTS = 

# -Wall -> warning all

.c.o: Makefile
	$(CC) $(COPTS) $<

all: $(OBJS) Makefile
	g++ $(LOPTS) -o laclare.x $(OBJS) -O3 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas -L/usr/lib -L/usr/local/lib/
