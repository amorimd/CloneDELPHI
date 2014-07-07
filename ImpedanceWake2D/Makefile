.SUFFIXES: .cc

CC = g++

PATH_ALGLIB = ../../../ALGLIB/mpfr/out

OPTS = -O2 -I$(PATH_ALGLIB) -L$(PATH_ALGLIB)

LIBOPTS = -lgmp -lmpfr -lmpalglib -lm -lstdc++ -lgsl -lgslcblas -L$(PATH_ALGLIB)

LOPTS = 



# -Wall -> warning all

all: round wakeround flat wakeflat lib Makefile

roundchamber.o: roundchamber.cc
	$(CC) $(LOPTS) -fPIC -o roundchamber.o -c roundchamber.cc -I$(PATH_ALGLIB)

round: roundchamber.o
	$(CC) $(LOPTS) -o roundchamber.x roundchamber.o -O2 $(LIBOPTS)

wake_roundchamber.o: wake_roundchamber.cc
	$(CC) $(LOPTS) -fPIC -o wake_roundchamber.o -c wake_roundchamber.cc -I$(PATH_ALGLIB)

wakeround: wake_roundchamber.o
	$(CC) $(LOPTS) -o wake_roundchamber.x wake_roundchamber.o -O2 $(LIBOPTS)

flatchamber.o: flatchamber.cc
	$(CC) $(LOPTS) -fPIC -o flatchamber.o -c flatchamber.cc -I$(PATH_ALGLIB)

flat: flatchamber.o
	$(CC) $(LOPTS) -o flatchamber.x flatchamber.o -O2 $(LIBOPTS)

wake_flatchamber.o: wake_flatchamber.cc
	$(CC) $(LOPTS) -fPIC -o wake_flatchamber.o -c wake_flatchamber.cc -I$(PATH_ALGLIB)

wakeflat: wake_flatchamber.o
	$(CC) $(LOPTS) -o wake_flatchamber.x wake_flatchamber.o -O2 $(LIBOPTS)

clean:
	rm -f *.o libIW2D.so *.x

lib: wake_roundchamber.o
	$(CC) ${LOPTS} -shared -o libIW2D.so wake_roundchamber.o -O2 $(LIBOPTS) 
