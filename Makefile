SHELL = /bin/sh

FC = ifort
FCFLAGS = -O3 -module .mod -I.mod 
debug: FCFLAGS += -g -check all -fpe0 -warn -traceback -debug extended
#FC = gfortran 
#FCFLAGS = -O3 -J .mod -I.mod 
#debug: FCFLAGS += -g -Wall -Wextra -Warray-temporaries \
	-fimplicit-none -fbacktrace -ffree-line-length-0 -fcheck=all \
	-ffpe-trap=invalid,zero,overflow,underflow -finit-real=nan
MKLLIB = -qmkl

UTILS := utils/constants.o utils/params.o utils/utils.o
ROUTINES := #{}.o

all: clean main
debug: clean main

main: main.f90 $(UTILS) $(ROUTINES)
	$(FC) $(FCFLAGS) -o $@ $^

$(ROUTINES): %.o : %.f90 $(UTILS)
	$(FC) $(FCFLAGS) -c -o $@ $<

$(UTILS) : %.o : %.f90
	$(FC) $(FCFLAGS) -c -o $@ $<

test.x: test.f90 $(UTILS)
	$(FC) $(FCFLAGS) -o $@ $^

.PHONY: clean delint debug

clean: # clean_molecules
	rm -f main .mod/* *.o */*.o */*/*.o 
	rm -f *.?-1 *-

delint: 
	rm -f *.i90 */*.i90 */*/*.i90
	rm -f .genmod/*