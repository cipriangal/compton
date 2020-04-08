EXE = main

#ROOTCFLAGS    = $(shell root-config --cflags)
#ROOTLIBS      = $(shell root-config --libs)
#ROOTGLIBS     = $(shell root-config --glibs)

#INCFLAGS = -I$(ROOTSYS)/include 
#LDFLAGS = -L$(ROOTSYS)/lib -lMathMore
INCFLAGS = -I/afs/rhic.bnl.gov/eic/restructured/env/EIC2018a/PACKAGES/FairSoft-may16.build/include/root
LDFLAGS = -L/afs/rhic.bnl.gov/eic/restructured/env/EIC2018a/PACKAGES/FairSoft-may16.build/lib/root

#FC = gfortran -g -m64
#CXX = g++ -fno-inline -m64
FC = gfortran -g -m32
CXX = g++ -g -fno-inline -m32
FLAGS = -Wall $(INCFLAGS) $(LDFLAGS)

LIB = -lgfortran -lGui -lCore -lCint -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -pthread -lm -ldl -rdynamic  #-lmathlib -lkernlib -lpacklib_noshift -ldl -lm# $(ROOTLIBS) 
#LIB = -lgfortran -L/cern/pro/lib -lmathlib -lpthread -lkernlib -lpacklib_noshift -ldl -lm $(ROOTLIBS) 

#lpthread  added

COMPILE = $(CXX) $(FLAGS) -c #-O0

all: $(EXE)

comrad.o: comrad.f
	$(FC) -c comrad.f
ddilog.o: ddilog.f
	$(FC) -c ddilog.f
ranmar.o: ranmar.f
	$(FC) -c ranmar.f

#main_in_C.o: main_in_C.cpp
#	$(COMPILE) main_in_C.cpp

main: ddilog.o ranmar.o comrad.o 
	$(FC)   -o $(EXE) ddilog.o ranmar.o comrad.o 

clean:
	rm *.o main
