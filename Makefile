#
# Makefile to build example MPI programs
#

CC=mpicc
CXX=mpiCC
FC=mpif77
F90=mpif90

COMP=GNU

ifeq ($(COMP), GNU)
  CFLAGS=-Wall
  FFLAGS=-Wall
endif

EXES=stencil


hello_world_c: stencil.c
	$(CC) $(CFLAGS) -o $@ $^

all: $(EXES)

.PHONY: clean all

clean:
	\rm -f $(EXES)
	\rm -f *.o
