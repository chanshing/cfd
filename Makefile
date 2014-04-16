FC = ifort
OMP = 
FFLAGS = -fast $(OMP) #-check bounds -traceback
LFLAGS = $(OMP)
OBJECTS = commonModules.o pointNeighbor.o laplace.o biconjugateGradient.o ns2DComp.ALE.o 
.PHONY: clean

ns: $(OBJECTS)
	$(FC) $(LFLAGS) $(OBJECTS) -o ns 

%.o: %.f90
	$(FC) $(FFLAGS) -c $<

clean:
	rm -f $(OBJECTS) ns *.mod
