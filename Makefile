all:
	 gcc -std=c99 -o cgSolver cgSolver.c -lm
clean: 
	rm *.o c
