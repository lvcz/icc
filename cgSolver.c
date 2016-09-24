#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#define M_PI 3.1415926535

	


int metodo_gradiente(double *a, double *b, int n, double tol, int maxIter,int nBandas,double *v_norma,double *erro, double *timeMin,double *timeMax,double *x_result);
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag );
double timestamp();
double *calcula_func_b(int n);
double vetorT_x_vetor(double *a, double *b ,int n);
void copia_vetor(double *a, double *b,int n);
double *aloca_vetor(int n);
void matriz_x_vetor(double *a, double *x,double *b, int n, int nBandas);
double vetor_x_num(double num, double *a,int n);
void vetor_mais_vetor(double *a, double *b, double *res, int n);
void imprime_saida(FILE *arquivo_saida,int k ,double *v_norma,double *erro,double timemin,double timemax,double timemedio,int n,double *x);


/***********************
 * N: tamanho do sistema linear
 * k: numero da diagonal, 0 = diagonal principal, 1 = acima/abaixo da diagonal, 2 = ...
 * nBandas: numero de bandas do sistema linear
 * diag: vetor para armazenar os valores da diagonal. Deve ser alocado por quem chama a função.
 ***********************/

int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag )
{
  if ( !diag || N < 3 || nBandas > N/2 || k < 0 || k > nBandas )
    return (-1);

  /* garante valor dominante para diagonal principal */
  double fator = (k == 0) ? ((double)(nBandas-1)) : (0.0);

  double invRandMax = 1.0 / (double)RAND_MAX;

  for (int i=0; i < N-k; ++i)
  {
    diag[k*N+i] = fator + (double)rand() * invRandMax;
  }

  return (0);
}
double timestamp(void){
    struct timeval tp;
    gettimeofday(&tp, NULL);
    return((double)(tp.tv_sec*1000.0 + tp.tv_usec/1000.0));
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
double vetorT_x_vetor(double *a, double *b, int n)
{
    double aux =0.0d;
    for (int i = 0 ; i<n ;++i) 
        aux +=a[i]*b[i];    
    return aux;
}
double vetor_x_num(double num, double *a,int n)
{
	for (int i = 0 ; i<n ;++i)
        a[i] = a[i]*num;
}
void vetor_mais_vetor(double *a, double *b, double *res, int n){
	
for (int i = 0 ; i<n ;++i) 
        res[i] = a[i] + b[i];
}

void vetor_menos_vetor(double *a, double *b, double *res, int n){
	
for (int i = 0 ; i<n ;++i) 
        res[i] = a[i] - b[i];
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

void matriz_x_vetor(double *a, double *x,double *b, int n, int nBandas)
{
	
    for (int i =0; i< n ;++i){
		b[i]= a[i] * x[i];
		
		for(int k=1; k <= nBandas; ++k){
			if (i+k < n)
				b[i] += a[i+k*n] *x[k+i]; 
			if(i-k >= 0)
				b[i] += a[(i-k)+(k*n)] +x[i-k];
		}
	}
}


void main (int argc, char** argv){

	
	int n;
	int nBandas = 0;
	int k = 0;
    int maxIter = 0;
	double *A = NULL; 
	double *diag = NULL;
	double *b = NULL;
	double *x = NULL;
	double tol;	
    FILE *arquivo_saida;
	srand(20162);	
	
	if(argc >2){
		
            int aux;
            aux = atoi(argv[1]);
            if (aux > 0 ) {
                n = aux;
            } else {
                fprintf (stderr, "\nN deve ser n > 0 ");
            return;
			}
				aux = atoi(argv[2]);
				if (aux >= 0) {
					nBandas = aux;
				} else {
					fprintf (stderr, "\n Erro ao definir o número de Bandas.");
				return ;
				if (nBandas > n/2)
				{
					fprintf (stderr, "\n Erro ao definir o número maximo de .");
				return ;
				}	
			}
        
				
		for (int i = 3; i < argc; i += 2) {
		
		   
		
			if (!strcmp("-i", argv[i]))
			{
				int aux;
				aux = atoi(argv[i+1]);
				if (aux >= 0) {
					maxIter = aux;
				} else {
					fprintf (stderr, "\n Erro ao definir o número maximo de iterações.");
				return ;
			  }
			} 
			
			if (!strcmp("-t", argv[i])) {
				double aux;
				aux = atof(argv[i+1]);
				if (aux > 0 && aux < 1) {
					tol = aux;
				} else {
					fprintf (stderr," \n Erro no valor da tolerancia");
					return ;
				}
			}
			
			if (strcmp("-o",argv[i]) == 0) {
				arquivo_saida = fopen(argv[i+1], "w");
				
				if (arquivo_saida == NULL) {
					fprintf (stderr,"\n Erro no arquivo de saída!\n");
					return ;  
				}    
			}
		}
	}
	else{
		 fprintf (stderr, "\n N deve ser n > 0 e nBandas < n/2 \n");
            return;
	}
	if(maxIter==0) maxIter = n;
	double *erro = aloca_vetor(n);
	double *v_norma = aloca_vetor(n);
	
		
	A= aloca_vetor(n+(nBandas*n));
	int i =0;
	do{
		generateRandomDiagonal(n,i,nBandas,A);
		++i;
		
	}while(i<nBandas+1);
	
	
	for (int i = 0 ; i<n+(nBandas*n);++i){
		printf("# A[%d]: %lf",i,A[i]);
		
	}
	
	
	b = calcula_func_b(n);
	double timeMin=1.0d;
	double timeMax=0.0d;
	
		double *x_result = aloca_vetor(n);
	double timeini = timestamp();
	
	k = metodo_gradiente(A,b,n,tol,maxIter,nBandas,v_norma,erro,&timeMin,&timeMax,x_result);
	double timefim = timestamp();
	
	double timemedio= (timefim-timeini)/k;
	
	imprime_saida(arquivo_saida,k ,v_norma,erro,timeMin,timeMax,timemedio,n,x_result);

	
}
void imprime_saida(FILE *arquivo_saida,int k ,double *v_norma,double *erro,double timeMin,double timeMax,double timemedio,int n,double *x){
	
	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"# Tempo Método CG: %lf %lf %lf \n" ,timeMin,timemedio,timeMax);
	fprintf(arquivo_saida,"#falta residuo \n");
	fprintf(arquivo_saida,"#\n");
	fprintf(arquivo_saida,"# Norma Euclidiana do Residuo e Erro aproximado\n");
	for( int i = 0; i < k ; ++i)
	{
		fprintf(arquivo_saida,"# i=%d: %lf %lf \n" ,i, v_norma[i],erro[i]);
	}
	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"%d\n",n);
	for( int i = 0; i < n ; ++i)
	{
		fprintf(arquivo_saida,"%lf " ,x[i]);
	}
	
}


/*
		Calcula metodo do gradiente seguindo o livro.
*/
int metodo_gradiente(double *a, double *b, int n, double tol, int maxIter,int nBandas,double *v_norma,double *erro, double *timeMin,double *timeMax,double *x_result)
{
	
	
	double *x_ant = aloca_vetor(n);	
	double *r = aloca_vetor(n);	
	double aux = vetorT_x_vetor(b,b,n);
    double aux1 = 0.0d;    
	double *v= aloca_vetor(n);
    copia_vetor(v,b,n);	
	double *z = aloca_vetor(n);	
	double *mr = aloca_vetor(n);	
	
	
	double s = 0.0d;	
	double m = 0.0d;
	double norma = 0.0d;
	double norma_ant = 0.0d;
	
	double v_x_z=0.0d;
	int k=0;
	while( k < maxIter )
	{	
		double timeIni = timestamp();
		
		//calcula z	 = A * v
		matriz_x_vetor(a,v,z,n,nBandas);

		//calcula  s= aux / v * z	
		//
		v_x_z = vetorT_x_vetor(v,z,n);	
		s = aux/v_x_z;
		//


		//calcula x(k+1) = x(k) + s*v		
		//
		copia_vetor(x_ant,x_result,n);		
		
		vetor_x_num(s,v,n);		
		vetor_mais_vetor(x_result,v,x_result,n);

		//

		norma = 0.0d;
		//calcula  r = r- sz
		vetor_x_num(s,z,n);
		vetor_menos_vetor(r,z,r,n);

		//aux1= r*r
		aux1 = vetorT_x_vetor(r,r,n);
		
		norma_ant= norma;
		
		norma = sqrt(aux1);
		v_norma[k]=norma;
	
		
		if (k>0){
			erro[k] = fabs(v_norma[k] - v_norma[k-1]);
		}
		else{
			erro[k] = norma;
		}	
		// caso convergiu retorna x-result
		
		if (erro[k] <tol){
			return k;
			} 
			
		m = aux1 / aux;
		aux = aux1;

		// v = r+m*r
		copia_vetor(mr,r,n);
		vetor_x_num(m,mr,n);
		vetor_mais_vetor(r,mr,v,n);
		k++;
		
		
		double timeFim =timestamp();
		
		double timediff= timeFim-timeIni;
		
		if( timediff < *timeMin){
			
			*timeMin = timediff;
			
			}
			
		if(timediff > *timeMax){
			
			*timeMax=timediff;
			
		}
	}
	
	return k;
}

