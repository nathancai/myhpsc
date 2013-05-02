# $MYHPSC/homework3/Makefile
# Nathan Juswanto

OBJECTS  = functions.o newton.o test1.o
OBJECTS2 = functions.o newton.o intersections.o
MODULES  = functions.mod newton.mod

FFLAGS = -g

.PHONY: test1 intersections clean 

test1: test1.exe
	./test1.exe

test1.exe: $(MODULES) $(OBJECTS)
	gfortran $(FFLAGS) $(OBJECTS) -o test1.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

intersections: intersections.exe
	./intersections.exe

intersections.exe: $(MODULES) $(OBJECTS2)
	gfortran $(FFLAGS) $(OBJECTS2) -o intersections.exe

%.o : %.f90
	gfortran $(FFLAGS) -c  $< 

%.mod: %.f90
	gfortran $(FFLAGS) -c $<

clean:
	rm -f *.o *.exe *.mod

