#include <stdlib.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <sys/time.h>
#include <float.h>

#define M_PI 3.1415926535

	


int metodo_gradiente(double *a, double *b, int n, double tol, int maxIter,int nBandas,double *v_norma,double *erro, double *timeMin,double *timeMax,double *x_result,double *timeResMin,double *timeResMax,double *timeResMedio);
void imprime_saida(FILE *arquivo_saida,int k ,double *v_norma,double *erro,double timeMin,double timeMax,double timemedio,int n,double *x,double timeResMin, double timeResMax, double timeResMedio);
int generateRandomDiagonal( unsigned int N, unsigned int k, unsigned int nBandas, double *diag );
double timestamp();
double *calcula_func_b(int n);
double vetorT_x_vetor(double *a, double *b ,int n);
void copia_vetor(double *a, double *b,int n);
double *aloca_vetor(int n);
void matriz_x_vetor(double *a, double *x,double *b, int n, int nBandas);
void vetor_x_num(double num, double *a,int n);
void vetor_mais_vetor(double *a, double *b, double *res, int n);
void vetor_menos_vetor(double *a, double *b, double *res, int n);



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
	double pi2= (2 * M_PI);
	 
	double *bx = NULL; 
	if ( ! ( bx = (double*) malloc(n*sizeof(double))) )
			return (NULL);
		
	
	for (int i = 0 ; i < n ; ++i)
		{ 
			
			double x = (i* M_PI/n);
			bx[i] = pi4 * (sin ( pi2 * x) + sin ( pi2 *(M_PI - x)));
		
		}
		
	return (bx);
}

double vetorT_x_vetor(double *a, double *b, int n)
{	//LIKWID_MARKER_INIT;
	//LIKWID_MARKER_START("MVV");
    double aux =0.0d;
    for (int i = 0 ; i<n ;++i){ 
        aux +=a[i]*b[i];}
    //LIKWID_MARKER_STOP("MVV");
    //LIKWID_MARKER_CLOSE;    
    return aux;
}

void vetor_x_num(double num, double *a,int n)
{
	for (int i = 0 ; i<n ;++i)
        a[i] = a[i]*num;
}

void vetor_mais_vetor(double *a, double *b, double *res, int n)
{	
	for (int i = 0 ; i<n ;++i) 
        res[i] = a[i] + b[i];
}

void vetor_menos_vetor(double *a, double *b, double *res, int n)
{	
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
    //for (int i = 0 ; i<n ;++i) ret[i] =0.0d;
    return ret;
}

void matriz_x_vetor(double *a, double *x,double *b, int n, int nBandas)
{
	//LIKWID_MARKER_START("MMV");
    for (int i =0; i< n ;++i){
		b[i]= a[i] * x[i];
		//printf(" (b[%d]) = a[%d] * x[%d] +", i,i,i);
        
		for(int k=1; k <= nBandas; ++k){
			if (i+k < n){
				b[i] += a[i+k*n] * x[k+i]; 
                //printf(" a[%d] * x[%d] +", i+k*n,k+i);
            }
			if(i-k >= 0){
				b[i] += a[(i-k) + (k*n)] * x[i-k];
				//printf(" a[%d] * x[%d] ", (i-k)+(k*n),i-k);
            }
		}
         //printf(")\n");
	}
	//LIKWID_MARKER_STOP("MMV");
}


void main (int argc, char** argv){
	
	
	int n;
	int nBandas = 0;
	int k = 0;
    int maxIter = 0;
	double *A = NULL; 
	double *b = NULL;
	double tol=0.0d;	
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
        if (aux > 0 && aux % 2==1 && aux< n/2 ) {
            nBandas =(int) aux;
        } else {
            fprintf (stderr, "\n Erro ao definir o número de Bandas( 0 < nBandas  <  n/2,  nBandas deve ser impar). \n");
        return ;
        }
        
		for (int i = 3; i < argc; i += 2) {
		
		   
		
			if (!strcmp("-i", argv[i]))
			{
				int aux;
				aux = atoi(argv[i+1]);
				if (aux >= 0 ) {
					maxIter = aux;
				} else {
					fprintf (stderr, "\n Erro ao definir o número maximo de iterações.\n");
				return ;
			  }
			} 
			
			if (!strcmp("-t", argv[i])) {
				double aux;
				aux = atof(argv[i+1]);
				if (aux > 0 && aux < 1) {
					tol = aux;
					
				} else {
					fprintf (stderr," \n Erro no valor da tolerancia\n");
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
		 fprintf (stderr, "\n Uso: cgSolver n nBandas -i <maxIter> -t <tolerancia> -o <arquivo_saida>   \n");
            return;
	}
	maxIter = maxIter==0? maxIter = n:maxIter;	
	double *erro = aloca_vetor(maxIter);
	double *v_norma = aloca_vetor(maxIter);
	nBandas=((nBandas - 1)/2);
  
    
	A= aloca_vetor(n+(nBandas*n));
	
	for( int i =0; i <= nBandas; ++i){
		generateRandomDiagonal(n,i,(nBandas+1)*2,A);
	}
	
	b = calcula_func_b(n);
	double timeMin=1000000.0d;
	double timeMax=0.0d;
	double *x_result = aloca_vetor(n);
	double timeini = timestamp();	
    double timeResMin; 
    double timeResMax;
    double timeResMedio;
	k = metodo_gradiente(A,b,n,tol,maxIter,nBandas,v_norma,erro,&timeMin,&timeMax,x_result,&timeResMin,&timeResMax,&timeResMedio);
	double timefim = timestamp();
	
	double timemedio= (timefim-timeini)/k;
	timeResMedio=timeResMedio/k;
	imprime_saida(arquivo_saida,k ,v_norma,erro,timeMin,timeMax,timemedio,n,x_result,timeResMin,timeResMax,timeResMedio);

    //free()
   //LIKWID_MARKER_CLOSE;
	
}

void imprime_saida(FILE *arquivo_saida,int k ,double *v_norma,double *erro,double timeMin,double timeMax,double timemedio,int n,double *x,double timeResMin, double timeResMax, double timeResMedio){
	
	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"# Tempo Método CG: %.14g %.14g %.14g \n" ,timeMin,timemedio,timeMax);
	fprintf(arquivo_saida,"# Tempo Resíduo: %.14g %.14g %.14g \n" ,timeResMin,timeResMedio,timeResMax);
	fprintf(arquivo_saida,"#\n");
	fprintf(arquivo_saida,"# Norma Euclidiana do Residuo e Erro aproximado\n");
	
	for( int i = 0; i < k ; ++i)
	{
		fprintf(arquivo_saida,"# i=%d: %.14g %.14g \n" ,i, v_norma[i],erro[i]);
	}
	fprintf(arquivo_saida,"###########\n");
	fprintf(arquivo_saida,"%d\n",n);
	for( int i = 0; i < n ; ++i)
	{
		fprintf(arquivo_saida,"%.14g " ,x[i]);
	}
	
}

int metodo_gradiente(double *a, double *b, int n, double tol, int maxIter,int nBandas,double *v_norma,double *erro, double *timeMin,double *timeMax,double *x_result,double *timeResMin,double *timeResMax,double *timeResMedio)
{
	
	
	double *x_ant = aloca_vetor(n);	
	//r=b
	double *v= aloca_vetor(n);
	copia_vetor(v,b,n);	
	
	//r=b
	double *r = aloca_vetor(n);	
	copia_vetor(r,b,n);
	
	double aux = vetorT_x_vetor(b,b,n);
    double aux1 = 0.0d;    
    
    
	double *z = aloca_vetor(n);	
	double *mr = aloca_vetor(n);	
	double *sv = aloca_vetor(n);	
	
	
	double s = 0.0d;	
	double m = 0.0d;
	double norma = 0.0d;
	//double norma_ant = 0.0d;
	
	double v_x_z=0.0d;
	int k=0;
	
	while( k < maxIter )
	{	
		double timeIni = timestamp();		
        //calcula z	 = A * v       
        matriz_x_vetor(a,v,z,n,nBandas);

		//calcula  s= aux / v * z	

		v_x_z = vetorT_x_vetor(v,z,n);	

		if (v_x_z < DBL_EPSILON || v_x_z > DBL_EPSILON )
			s = aux / v_x_z;
		else{
			fprintf(stderr,"Divisao por 0");
			return k;
		}
        //
    	//calcula x(k+1) = x(k) + s*v		
		//
        copia_vetor(x_ant,x_result,n);		
		copia_vetor(sv,v,n);
        vetor_x_num(s,sv,n);	
        vetor_mais_vetor(x_ant,sv,x_result,n);
        //
		//norma = 0.0d;
		double timeIniRes= timestamp();
        //calcula  r = r- sz
		
        vetor_x_num(s,z,n);
        vetor_menos_vetor(r,z,r,n);
		
		//aux1= r*r
        
        aux1 = vetorT_x_vetor(r,r,n);
		
		if (k == 0){
			// aux na primeira iteracao contem b² = r² originial e nao modificado 
			norma = sqrt(aux);
		}
		else 
		{		
			//norma_ant = norma;
			norma = sqrt(aux1);			
		}
		
        v_norma[k]=norma;	

        double timeFimRes =timestamp();
		
		if (k>0){
			erro[k] = fabs(v_norma[k] - v_norma[k-1]);
		}
		else{
			erro[k] = norma;
		}	
		// caso convergiu retorna x-result
		
		if(tol>0.0d){
			if (fabs(erro[k] - tol) < DBL_EPSILON){  
				
				return k;
			} 
		}
		// m=aux1/aux; aux =aux1;	
        
		if (aux < DBL_EPSILON || aux > DBL_EPSILON )
			m = aux1 / aux;
		else{
			fprintf(stderr,"Divisao por 0");
			return k;
		}
		
		aux = aux1;
		// v = r+m*r
		copia_vetor(mr,r,n);
		vetor_x_num(m,mr,n);	
		vetor_mais_vetor(r,mr,v,n);
		
		k++;
		
		
		double timeFim =timestamp();
		double timediff= timeFim-timeIni;
		double timeResdiff= timeFimRes-timeIniRes;
		// iteracao
        *timeResMedio += timeResdiff;
		if( timediff <= *timeMin){
			
			*timeMin = timediff;
			
			}
			
		if(timediff >= *timeMax){
			
			*timeMax=timediff;
			
		}
        
        
        //res
        if( timeResdiff <= *timeMin){
			
			*timeResMin = timeResdiff;
			
			}
			
		if(timeResdiff >= *timeResMax){
			
			*timeResMax=timeResdiff;
			
		}
        
        
		
	}
           
	return k;
}
