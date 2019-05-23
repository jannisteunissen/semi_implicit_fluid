FC := gfortran
FFLAGS := -Wall -O2 -fcheck=all

.PHONY: all clean

all: semi_implicit

clean:
	$(RM) semi_implicit

semi_implicit: m_config.o

%: %.f90
	$(FC) -o $@ $^ $(FFLAGS)

%.o: %.f90
	$(FC) -c -o $@ $< $(FFLAGS)
