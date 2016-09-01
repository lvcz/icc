#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>

double *generateSquareRandomPositiveDefiniteMatrix( unsigned int n );
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int kMax, double *diag );
void imprime_sistema(double *A, double *b, int n) ;





/***********************
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * kMax: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/


int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int kMax, double *diag )
{
  if ( !diag || N < 3 || kMax > N/2 || k < 0 || k > kMax )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator = (k == 0) ? ((double)(kMax-1)) : (0.0);

  double invRandMax = 1.0 / (double)RAND_MAX;

  for (int i=0; i < N-k; ++i)
  {
    diag[i] = fator + (double)rand() * invRandMax;
  }

  return (0);
}




void main (int argc, int** argv){

	int N = 1;
	double *A = NULL, 
	*b = NULL , *x = NULL;
	double erro;


	printf("/nDigite o tamanho da matriz: ");
	scanf("%d",&N);
	srand(20142);
	//printf("/n Digite o erro:");
	//scanf("%la",&erro);


	A = generateSquareRandomPositiveDefiniteMatrix(N);
	x = (double*) malloc(N*sizeof(double));
	b = (double*) malloc(N*sizeof(double));

	double invRandMax = 1.0 / (double)RAND_MAX;
	for (int i = 0 ; i<N ; ++i){
		b[i] = (double)rand()*invRandMax;
	}
	
	imprime_sistema(A,b,N);
}


double *generateSquareRandomPositiveDefiniteMatrix( unsigned int n )
{
  double *mat = NULL;

  /* return NULL if memory allocation fails */
  if ( ! (mat = (double *) malloc(n*n*sizeof(double))) )
    return (NULL);

  /* generate a randomly initialized matrix in row-major order */
  double *ptr = mat;
  double *end = mat + n*n;

  double invRandMax = 1.0/(double)RAND_MAX;

  while( ptr != end ) {
    *ptr++ = (double)rand() * invRandMax;
  }

  /* Now we want to make this matrix positive definite. Since all
     values are positive and <1.0, we have to ensure it is symmetric
     and diagonally dominant.

     A = A + transpose(A)
     A = A + I*n                        */
  unsigned int i,j;
  for (i=0; i<n; ++i)
    for (j=i+1; j<n; ++j)
    {
      double aux = mat[i*n+j];
      mat[i*n+j] += mat[j*n+i];
      mat[j*n+i] += aux;
    }

  for (i=0; i<n; ++i)
    mat[i*n+i] += mat[i*n+i] + n;

  return (mat);
}


void imprime_sistema(double *A, double *b, int n) {
  
  
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      printf( "%lfx%d ", A[i*n+j], j);
    }
    printf( "= %lf\n", b[i]);
  }
}


