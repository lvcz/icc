#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

#define TAMANHO_MAXIMO 32768
#define FALSE 0
#define TRUE 1
#define ERRO -1
#define erro 000.1


void metodo_gradiente(double *A, double *b, double *x, int n);
void calcula_residuo(double *A, double *b, double *x, double *r, int n);
double calcula_norma_residual(double *r, int n);
int generateRandomPositiveDefiniteLinearSystem( int N, double *A, double *b );
void ler_parametros(int argc, char *argv[], double erro, int n, int max_iter,FILE *arquivo_entrada,FILE arquivo_saida );
double timestamp(void);

#endif
