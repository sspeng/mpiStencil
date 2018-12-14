stencil: stencil.c
	mpicc -std=c99 -Ofast $^ -o $@

