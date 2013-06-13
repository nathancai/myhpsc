
OBJECTS3 = functions.o quadrature.o test3.o
FFLAGS = 
NUM_PROCS ?= 4   # default if not specified on command line or env variable

.PHONY: test3 clean 

test3: test3.exe
	mpiexec -n $(NUM_PROCS) ./test3.exe

test3.exe: $(OBJECTS3)
	mpif90 $(FFLAGS) $(OBJECTS3) -o test3.exe

%.o : %.f90
	mpif90 $(FFLAGS) -c  $< 

clean:
	rm -f *.o *.exe *.mod

