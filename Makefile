FC = ifort
OMP = -openmp 
FFLAGS = -fast -fpp $(OMP) 
LFLAGS = $(OMP) 
OBJECTS = commonModules.o dataLoader.o pointNeighbor.o laplace.o biconjugateGradient.o subrutinas.o ns2DComp.ALE.o 
#MODULES = mvelocidades.mod mvariabgen.mod mvariables.mod mestabilizacion.mod mnewmark.mod matrices.mod mat2.mod timers.mod melementsurrpoint.mod mpointsurrpoint.mod mlaplace.mod mbiconjgrad.mod inputdata.mod meshdata.mod

.PHONY: clean

all: ns

xhost: FFLAGS += -xHost 
xhost: ns

debug: FFLAGS += -check bounds -traceback -warn nodeclarations -warn unused
debug: ns

profile: FFLAGS += -p
profile: ns

mkl: FFLAGS += -mkl
mkl: LFLAGS += -mkl
mkl: ns

ns: $(MODULES) $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o ns 

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

#%.mod: %.f90
#	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) ns *.mod
