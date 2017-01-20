all:
	 gcc -std=c99 -o -O3 -navx -march=native cgSolver cgSolver.c -lm
clean: 
	rm *.o c
