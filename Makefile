FC = ifort
OMP = -openmp 
FFLAGS = -fast -fpp -mkl $(OMP) #-fpe0 
LFLAGS = $(OMP) -mkl
OBJECTS = commonModules.o dataLoader.o pointNeighbor.o \
		  mLaplace.o biconjGrad.o smoothing.o calcRHS.o \
		  subrutinas.o meshMove.o ns2DComp.ALE.o 

.PHONY: clean

all: ns

xhost: FFLAGS += -xHost 
xhost: ns

debug: FFLAGS += -check bounds -traceback -warn
debug: ns

idb: FFLAGS = -debug -O0 -fpp $(OMP)
idb: LFLAGS = -debug $(OMP)
idb: ns

profile: FFLAGS += -p
profile: ns

ns: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o ns 

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) ns *.mod
