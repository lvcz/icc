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


double timestamp(void){
      struct timeval tp;
      gettimeofday(&tp, NULL);
      return((double)(tp.tv_sec + tp.tv_usec/1000000.0));
}

void calcula_residuo(double *A, double *b, double *x, double *r, int n) {

  double soma;

  for (int i = 0; i < n; ++i) {
    soma = 0.0;
    for (int j = 0; j < n; ++j) {
      soma += A[i*n+j]*x[j];
    }
    r[i] = b[i]-soma;
  }
}

double calcula_norma_residual(double *r, int n) {

  double soma = 0.0;

  for (int i = 0; i < n; ++i) {
    soma += r[i]*r[i];
  }

  return sqrt(soma);
}

int metodo_gradiente(double *A, double *b, double *x,  int n, FILE *arquivo_saida,double erros, double max_iter) {

  int k = 0;
fprintf(stdout,"23");
  // residuo
  double *r = (double *) malloc(n*sizeof(double));
  double norma_residual_ant, norma_residual_atual;
  // variaveis referentes a solucao do sistema
  double s = 0.0, soma = 0.0, erro_relativo = 0.0, erro_ant = 0.0;
  // variaveis referentes aos tempos deexecucao do metodo
  double start = 0.0, end = 0.0, start_erro = 0.0, total_erro = 0.0;

  if (!r) {
    fprintf(stderr, "\nProblemas ao alocar memória 47!\n");
    return ERRO;
  }

  // calculo do residuo e norma do residuo
  start = timestamp();
  start_erro = timestamp();
  calcula_residuo(A, b, x, r, n);
  norma_residual_atual = calcula_norma_residual(r, n);
  total_erro += timestamp() - start_erro;
  

  do {
    // calculo do lambida
    s = 0.0;
    soma = 0.0;

    for (int i = 0; i < n; ++i) {
      for (int j = 0; j < n; ++j) {
        soma += A[i*n+j]*r[j]*r[i];
      }
      s += r[i]*r[i];
    }
    
    // Evitar divisão por zero
    if (soma > -DBL_MAX) {
      s = s/soma;

    } else {
      fprintf(stderr, "\nDivisão por Zero!\n");
      return ERRO;
    }

    // calculo do novo x
    for (int i = 0; i < n; ++i) {
      x[i] = x[i] + s*r[i];
    }

    norma_residual_ant = norma_residual_atual;

    // calculo do residuo e norma do residuo
    start_erro = timestamp();
    calcula_residuo(A, b, x, r, n);
    norma_residual_atual = calcula_norma_residual(r, n);
    total_erro += timestamp() - start_erro;

    // erro relativo = |rn-1 - rn|/rn
    erro_relativo = sqrt((norma_residual_ant-norma_residual_atual)*(norma_residual_ant-norma_residual_atual))/norma_residual_atual;
    //fprintf(parametros.arquivo_saida, "erro_relativo %.17g\n", erro_relativo);
	erro_ant=erro_relativo;
  } while ((erro_relativo > erros) & (++k < max_iter));

  end = timestamp();
  free(r);
  fprintf(arquivo_saida, "\n");
  fprintf(arquivo_saida, " Erro: %.17g\n", norma_residual_ant);
  fprintf(arquivo_saida, " Tempo Grad: %.17g\n", end-start);
  fprintf(arquivo_saida, " Tempo Erro: %.17g\n", total_erro/k);
  fprintf(arquivo_saida, "\n");
  fprintf(arquivo_saida, "%d\n", n);
  for (int i = 0; i < n; ++i) {
    fprintf(arquivo_saida, "%.17g ", x[i]);
  }
  fprintf(arquivo_saida, "\n");
  return 0;
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
    b[i] = (double)rand() * invRandMax;

  return (0);
}

 void ler_sistema(double *A, double *b, FILE *arquivo_entrada, int n) {

  // ler matriz representada em uma dimensao
  for (int i = 0; i < n; ++i) {
    for (int j = 0; j < n; ++j) {
      fscanf(arquivo_entrada, "%lf", &A[i*n+j]);
    }
  }

  // leitura do vetor b
  for (int i = 0; i < n; ++i) {
    fscanf(arquivo_entrada, "%lf", &b[i]);
  }
}




int main(int argc, char *argv[]){
      double *A = NULL;
      double *b = NULL;
      double *x = NULL;
      FILE *arquivo_entrada, *arquivo_saida;
	int n,flag;
	int max_iter;
	double erros;
	
      srand( 20142 );
  double invRandMax = 1.0/(double)RAND_MAX;
      for (int i = 1; i < argc; i += 2) {
    
    if (strcmp("-r", argv[i]) == 0) {
			int aux;
      aux = atoi(argv[i+1]);
      if (aux > 0 && aux <= TAMANHO_MAXIMO) {
        n = aux;
      } else {
        fprintf (stderr, "\nN deve ser 0 < n <= %d\n", TAMANHO_MAXIMO);
        return ERRO;
      }
    } 
    
    
    if (strcmp("-e", argv[i]) == 0) {
	double aux;
      aux = atof(argv[i+1]);
      if (aux > 0 && aux < 1) {
        erros = aux;
      } else {
        fprintf (stderr,"\nErro deve ser 0 < erro < 1\n");
        return ERRO;
      }
    }
    
    
    if (strcmp("-k", argv[i]) == 0) {
		int aux;
      aux = atoi(argv[i+1]);
      if (aux >= 0 && aux < (2*TAMANHO_MAXIMO)) {
        max_iter = aux;
      } else {
        fprintf (stderr, "\nMáximo de Iterações deve ser 0 <= max_iter < %d\n", 2*TAMANHO_MAXIMO);
        return ERRO;
      }
    }
    
   
    if (strcmp("-i",argv[i]) == 0) {
      arquivo_entrada=fopen(argv[i+1], "r");
        
      if (arquivo_entrada != NULL){
		  fscanf(arquivo_entrada, "%d", &n);
     	}
     }
    
		
		
 
      
     
     
    
        if (strcmp("-o",argv[i]) == 0) {
        (arquivo_saida = fopen(argv[i+1], "w"));
      	if (arquivo_saida == NULL) {
        fprintf (stderr,"\nErro no arquivo de saída!\n");
        return ERRO;  
		}    
		}
}
      
      
	
	A = (double *) malloc(n*n*sizeof(double));
    b = (double *) malloc(n*sizeof(double));
    x = (double *) calloc(n,sizeof(double));
      if (!A || !b || !x) {
            fprintf(stderr, "\nProblemas ao alocar memóriaa!\n");
            return ERRO;
      }
      	
      if(arquivo_entrada){
		  printf("passei");
      ler_sistema(A,b,arquivo_entrada,n);}
      else
      printf("pass2/n");
       flag=generateRandomPositiveDefiniteLinearSystem(n,A, b); 
       

    flag=metodo_gradiente(A, b, x,  n,arquivo_saida,erros,max_iter);
    
return flag;
}
