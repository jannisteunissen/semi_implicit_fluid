FC := gfortran
FFLAGS := -Wall -O2 -fcheck=all

.PHONY: all clean

all: semi_implicit

clean:
	$(RM) semi_implicit

%: %.f90
	$(FC) -o $@ $< $(FFLAGS)
