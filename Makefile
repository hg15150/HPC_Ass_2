stencil: stencil.c
	mpicc -Wall -g -DDEBUG -o stencil stencil.c
