# Makefile for HFMM beam-beam package -- LINUX G77

.SUFFIXES: .cc

CC = g++ -c

OBJS = laclare.o

COPTS = -O2 -llapack -lblas
LOPTS = 

# -Wall -> warning all

.c.o: Makefile
	$(CC) $(COPTS) $<

all: $(OBJS) Makefile
	g++ $(LOPTS) -o laclare.x $(OBJS) -O2 -lm -lgfortran -lgsl -lgslcblas -llapack -lblas

clean:
	rm -f $(OBJS) laclare.x

