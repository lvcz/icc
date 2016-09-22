#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#define M_PI 3.1415926535



double *generateSquareRandomPositiveDefiniteMatrix( unsigned int n );
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int kMax, double *diag );
void imprime_sistema(double *A, double *b, int n) ;
double *calcula_func_b(int n);
void *GaussSeidel(double *A,double *b,double *x, int n);
double vetorT_x_vetor(double *a, double *b ,int n);
void copia_vetor(double *a, double *b,int n);
double *aloca_vetor(int n);
void matriz_x_vetor(double *a, double *v,double *z, int n);
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




void main (int argc, char** argv){

	int n;
	int kMax = 0;
	int k = 0;
    int maxIter = 0;
	//int retorno =0;
	double *A = NULL; 
	double *diag = NULL;
	double *b = NULL;
	double *x = NULL;
	double tol;
    FILE *arquivo_saida;
   srand(20162);	

	 if (!strcmp("n", argv[1]))
        {
            int aux;
            aux = atoi(argv[2]);
            if (aux > 0 ) {
                n = aux;
            } else {
                fprintf (stderr, "\nN deve ser n > 0 ");
            return;
          }
        }
		else {
            fprintf (stderr, "\nN deve ser n > 0 ");
            return;
        }		
    
    
    
    
	 for (int i = 3; i < argc; i += 2) {
    
       
    
        if (!strcmp("-i", argv[i]))
        {
            int aux;
            aux = atoi(argv[i+1]);
            if (aux > 0) {
                maxIter = aux;
            } else {
                fprintf (stderr, "\n Erro ao definir o numero maximo de iterações.");
            return ;
          }
        } 
        
        if (!strcmp("-t", argv[i])) {
            double aux;
            aux = atof(argv[i+1]);
            if (aux > 0 && aux < 1) {
                tol = aux;
            } else {
                fprintf (stderr,"Erro no valor da tolerancia");
                return ;
            }
        }
        
        if (strcmp("-o",argv[i]) == 0) {
            arquivo_saida = fopen(argv[i+1], "w");
            if (arquivo_saida == NULL) {
                fprintf (stderr,"\nErro no arquivo de saída!\n");
                return ;  
            }    
        }
    }

	
	
	A = generateSquareRandomPositiveDefiniteMatrix(n);
	//retorno = generateRandomDiagonal(n, n, kMax, diag);
	
	x = (double*) malloc(n*sizeof(double));
	for (int i = 0; i < n; ++i) {
		x[i] = 0.0;
	}

	

	//double invRandMax = 1.0 / (double)RAND_MAX;
	b = calcula_func_b(n);
	
	
	
	
	
	
	
	imprime_sistema(A,b,n);
	
	GaussSeidel(A,b,x,n);
	
	for (int i = 0; i < n; ++i) {
		printf( "x(%d)= %lf \n", i,x[i]);
	}
	
}

double *calcula_func_b(int n)
{
	
	double pi4 = 4 * ( M_PI * M_PI);
	
	 
	double *bx = NULL; 
	if ( ! ( bx = (double*) malloc(n*sizeof(double))) )
			return (NULL);
		
	
	for (int i = 0 ; i < n ; ++i)
		{ 
			
			double x = (i* M_PI/n);
			bx[i] = pi4 * (sin ( 2 * M_PI * x) + sin ( 2* M_PI *(M_PI - x)));
		
		}
		
	return (bx);
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


void imprime_sistema(double *A, double *b, int n) 
{
  
  
  for (int i = 0; i < n; ++i) {
	for (int j = 0; j < n; ++j) {
		printf( "%lfx%d ", A[i*n+j], j);
		}
    printf( "= %lf \n", b[i]);
  }
}


void *GaussSeidel(double *A,double *b,double *x, int n)
{
	int k =1;
	double norma=10;
	double *x_ant = (double*) malloc(n*sizeof(double));
	do
	{	
	for (int i = 0; i < n; ++i) {
		
		x[i] = b[i];
		
		for (int j = 0; j < i; ++j) 
		{
			if( j !=i)
			{				
			x[i] = x[i] -  (A[i*n+j] * x[j]);
			}
		}
		for (int j = i ; j<n; ++j)
		{
			if(j != i)
			{
			x[i] = x[i] - A[i*n+j] * x[j];
			}
		}
		x[i] =x[i] / A[i*n+i];
	}
	
	k++;
	
	
	 for (int i = 0; i < n; ++i) {
		if(k>1)		norma = x_ant[i] - x[i];
		if (norma <0 ) norma = norma *(-1);
		x_ant[i]=x[i];
		printf( "x(%d)= %lf \n", i,x[i]);
		
	}
	printf("\n norma = %lf,n = %d \n", norma, k);
	
}while (norma > 0.0d);
		return NULL;
}

double vetorT_x_vetor(double *a, double *b, int n)
{
    double aux =0.0d;
    for (int i = 0 ; i<n ;++i) 
        aux +=a[i]*b[i];
    
    return aux;
}

void copia_vetor(double *a, double *b,int n)
{
    for (int i = 0 ; i<n ;++i)
        a[i] = b[i];
}

double *aloca_vetor(int n)
{
    double *ret;
	if ( ! ( ret = (double*) malloc(n*sizeof(double))) )return NULL;
    for (int i = 0 ; i<n ;++i) ret[i] =0.0d;
    return ret;
}

void matriz_x_vetor(double *a, double *v,double *z, int n)
{
    for (int i =0; i< n ;++i)
	{
		for(int j = 0;j<n;j++)
		{
			z[i] += a[i*n+j] * v[j];
		}			
	}
}


/*
		Calcula metodo do gradiente seguindo o livro.
*/
void metodo_gradiente(double *a, double *b, double *x, int n, double tol, int maxIter)
{
	double *x_result = aloca_vetor(n);
	
	double *x_ant = aloca_vetor(n);
	
	double *r = aloca_vetor(n);
	
	double aux = vetorT_x_vetor(b,b,n);
    double aux1 = 0.0d;
    
	double *v= aloca_vetor(n);
    copia_vetor(v,b,n);
	
	double *z = aloca_vetor(n);	
	
	double s = 0.0d;
	
	double m = 0.0d;
	
	double norma = 0.0d;
	
	int k=0;
	
	
	
	while( k < maxIter )
	{	
		//calcula z	 = A * v
		matriz_x_vetor(a,v,z,n);
		
		//calcula  s= aux / v * z		
		for (int i =0; i< n ;++i)
		{
			//s[i] = aux[i] / v[i] * z[i];
		}
		//calcula x(k+1) = x(k) + s*v
		for (int i = 0 ; i < n;++i)
		{
			//x_ant[i] =  x_result[i];
			//x_result[i] +=  + s[i] *v[i];			
		}
		norma = 0.0d;
		//calcula  r = r- sz
		for (int i = 0 ; i < n;++i)
		{
			//r[i] -= s[i] * z[i];
			//aux1[i] = r[i] * r[i];
			//norma += aux1[i];
		}
		// caso convergiu retorna x-result
		//if (norma <tol){return;} 
		
		//caso não converja
		// m= aux1/aux; aux = aux1; v = r+m*r
		
		for( int i = 0;i<n ;++i){
			//m[i] = aux1[i] / aux[i];
			//aux[i] = aux1[i];
			
			//v[i] = r[i] + m[i] * r[i];
			
		}
		
		
		
		
		
		
	}
	
	
	
	
	
	
}
