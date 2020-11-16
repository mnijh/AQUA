FC = gfortran
FCOPTIONS = -O3 -pedantic

default: aqua_eos.o

aqua_eos.o: aqua_eos.f95
	$(FC) $(FCOPTIONS) -c $<

sample_implementation: sample_implementation.f95 aqua_eos.o
	$(FC) $(FCOPTIONS) -o $@ $^ 

test: test.f95 aqua_eos.o
	$(FC) $(FCOPTIONS) -o $@ $^

PHONY.: clean

clean :
	-rm -rf *.mod *.o sample_implementation
