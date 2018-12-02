stencil: stencil.c
	icc -std=c99 -Wall -Ofast -xHOST stencil.c -o stencil
	./stencil 1024 1024 100

big: stencil.c
	icc -std=c99 -Wall -Ofast stencil.c -o stencil
	./stencil 4096 4096 100

gpr: stencil.c
	icc -std=c99 -pg -fopenmp -g -o stream.gprof stencil.c
	./stream.gprof 1024 1024 100
	gprof -b -l stream.gprof gmon.out

gcc: stencil.c
	gcc -std=c99 -Wall -Ofast stencil.c -o stencil
	./stencil 1024 1024 100

vec: stencil.c
	icc -std=c99 -Wall $^ -o $@ -qopt-report=1 -qopt-report-phase=vec
	vim stencil.optrpt
check: 
	python check.py --ref-stencil-file stencil_8000_8000_100.pgm --stencil-file stencil.pgm --verbose
