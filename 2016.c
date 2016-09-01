 #include <stdlib.h>



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




void main (int argc int** argv){

	int N=1;
	double *a=NULL, *b=NULL,*x=NULL;
	double erro;


	printf("/nDigite o tamanho da matriz: ");
	scanf("%d",&N);
	srand(20142);
	printf("/n Digite o erro:")
	scanf("%la",&erro);


	A=generateSquareRandoPositiveDifneMatrix(N);
	x= (double*) malloc(N*sizeof(double));
	b= (double*) malloc(N*sizeof(double));

	invRandMax= 1.0/(double)RAND_MAX;
	for (i=0;i<N,++i){
		b[i]=(double)rand()*invRandMax;
	}
}

int generateRandomPositiveDefiniteLinearSystem( unsigned int N, double *A, double *b )
{
  if ( !A || !b )
    return (-1);

  /* generate a randomly initialized matrix in row-major order */
  double *ptr = A;
  double *end = A + N*N;

  double invRandMax = 1.0/(double)RAND_MAX;

  while( ptr != end ) {
    *ptr++ = (double)rand() * invRandMax;
  }

  /* Now we want to make this matrix positive definite. Since all
     values are positive and <1.0, we have to ensure it is symmetric
     and diagonally dominant.

     A = A + transpose(A)
     A = A + I*N                        */
  unsigned int i,j;
  for (i=0; i<N; ++i)
    for (j=i+1; j<N; ++j)
    {
      double aux = A[i*N+j];
      A[i*N+j] += A[j*N+i];
      A[j*N+i] += aux;
    }

  for (i=0; i<N; ++i)
    A[i*N+i] += A[i*N+i] + N;

  /* create the vector of independent terms (b) */
  for (i=0; i<N; ++i)
    b[i] = (double)rand() * invRandMax

  return (0);
}'
