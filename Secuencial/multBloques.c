#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>
#include <stdbool.h>
#include <unistd.h>

#define BASENAME_CSV "result-secuencial.csv"
#define HEADER_CSV "time,sizeMatrix,sizeBlock\n"
#define MIN(a, b) ((a)<(b) ? (a) : (b))
#define MAX(a, b) ((a)>(b) ? (a) : (b))

/* Time in seconds from some point in the past */
double dwalltime();
void producto(double *A,double *B,double *C, int N);
void productoBloques(double *A,double *B,double *C, int N, int r);
void inicializarMatrix(double *S, int sizeMatrix);
void crearMatriz(double *S, int sizeMatrix);
void imprimeMatriz(double *S,int N,int r);
void imprimeMatriz2(double *M, int N);
void imprimeVector(double *S, int sizeMatrix);
long get_integer_arg(int argc, char* argv[], int arg_index, long min_val, const char* description, const char* usage_msg, bool print_flag, void (*fun) (void) );
FILE* makeOutfile(char* basename);

int main (int argc, char *argv[]){

 double *A; // Matriz A
 double *B; // Matriz B
 double *C; // Matriz C
 double t1, t2;

 /* argumentos de linea de comandos */
 char *usage_msg = "Use %s <sizeMatrix> <sizeBlock> <print?>\n";
 int N = get_integer_arg(argc, argv, 1, 1, "sizeMatrix", usage_msg, true, NULL);
 int r = get_integer_arg(argc, argv, 2, 1, "sizeBlock", usage_msg, true, NULL);
 bool imprimir = get_integer_arg(argc, argv, 3, 0, "print", usage_msg, true, NULL);

 int sizeMatrix = N*N; //cantidad total de datos matriz
 int i;

 A= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para A
 B= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para B
 C= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para C

 crearMatriz(A, sizeMatrix);	//Inicializa A 
 crearMatriz(B, sizeMatrix);  //Inicializa B
 inicializarMatrix(C, sizeMatrix); 

 t1 = dwalltime();
 productoBloques(A,B, C, N, r);
 t2 = dwalltime();

 if (imprimir){
    printf("\n\n  A: \n" );
    imprimeMatriz2(A, N);

    printf("\n\n B: \n" );
    imprimeMatriz2(B, N);

    printf("\n\n  C: \n" );
    imprimeMatriz2(C, N);
 }

 printf("\nDuración total de la multiplicacion de matrices %4f segundos\n", t2-t1);

 FILE* outfile = makeOutfile(BASENAME_CSV);
 fprintf(outfile, "%4f,%d,%d\n", t2-t1, N, r);
 fclose(outfile);

 free(A);
 free(B);
 free(C);

 return 0;
} //FIN MAIN

void producto(double *A,double *B,double *C, int N){
  int i, j, k;
  for(i=0;i<N;i++){
    for(j=0;j<N;j++){
      C[i*N+j] = 0.0;
      for(k=0;k<N;k++){
        C[i*N+j] += A[i*N+k] * B[k*N+j];
      }
    }
  }   
}

void productoBloques(double *A,double *B,double *C, int N, int r){
  int kk, jj, i, j, k; 
  double temp;
  for(jj=0;jj<N;jj+= r){
        for(kk=0;kk<N;kk+= r){
                for(i=0;i<N;i++){
                        for(j = jj; j < MIN(N, jj+r); j++){
                                temp = 0.0;
                                for(k = kk; k <  MIN(N, kk+r); k++){
                                        temp += A[i*N+k]*B[k*N+j];
                                }
                                C[i*N+j] += temp;
                        }
                }
         }
   } 
}

void multBloques(double *A, double *B, double *C, int N, int p, int r) {
  int kk, jj, i, j, k;
  double temp;
  for (jj = 0; jj < N; jj += r) {
    for (kk = 0; kk < N; kk += r) {
      for (i = 0; i < N / p; i++) {
        for (j = jj; j < MIN(N, jj + r); j++) {
          temp = 0.0;
          for (k = kk; k <  MIN(N, kk + r); k++) {
            temp += A[i * N + k] * B[k * N + j];
          }
          C[i * N + j] += temp;
        }
      }
    }
  }
}


void inicializarMatrix(double *S, int sizeMatrix){
  int i;
  for (i=0; i<sizeMatrix ;i++){
    S[i]=0.0;
  };
}

void crearMatriz(double *S, int sizeMatrix){
  int i;
  for(i=0 ;i<sizeMatrix;i++){
  	S[i] = i;
  };//end i
}

void imprimeVector(double *S, int sizeMatrix){
  int i;
  printf("\n ");
  for(i=0 ;i<sizeMatrix;i++)
	   printf(" %.2f " ,S[i]);
 
  printf("\n\n ");
}

void imprimeMatriz2(double *M, int N){
  int i,j;
  printf("Contenido de la matriz:  \n ");
  for(i=0 ;i<N;i++){
    for(j=0; j<N;j++)
      printf(" %.2f " ,M[i*N+j]);
    printf("\n ");
  }
  printf("\n");
}

void imprimeMatriz(double *S,int N,int r){
// Imprime la matriz pasada por parametro
// N es la cantidad de bloques, r dimension de cada bloque
  int i,j,I,J,despB;

  printf("Contenido de la matriz: \n" );
  for (I= 0; I< N; I++){
    //para cada fila de bloques (I)
    for (i= 0; i< r; i++){
       for(J=0;J<N;J++){
          despB=(I*N+J)*r*r;
          for (j=0;j<r;j++){
              printf("%.2f ",S[despB+ i*r+j]);
          };//end for j
        };//end for J
        printf("\n ");
    };//end for i

  };//end for I
  printf(" \n\n");
}

long get_integer_arg(int argc, char* argv[], int arg_index, 
        long min_val, const char* description, const char* usage_msg,
        bool print_flag, void (*fun) (void) )
{
    if (arg_index >= argc) {
        if (print_flag)
            fprintf(stderr, usage_msg, argv[0]);
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    char* end;
    long rval = strtol(argv[arg_index], &end, 10);
    if (*end != '\0') {
        if (print_flag) {
            fprintf(stderr, "%s debe ser un número entero\n", description);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    if (rval < min_val) {
        if (print_flag) {
            fprintf(stderr, "%s debe ser al menos %ld\n", description, min_val);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    return rval;
}

FILE* makeOutfile(char* basename) {
    FILE *outfile;
    if( access( basename, F_OK ) != -1 ) {
       outfile = fopen(basename, "a");
    } else {
       outfile = fopen(basename, "w");
       if(outfile != NULL)
            fprintf(outfile, HEADER_CSV);
    }
    if (outfile == NULL)
        fprintf(stderr, "Cannot open outfile %s\n", basename);
    return outfile;
}

/*****************************************************************/

#include <sys/time.h>

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv,NULL);
	sec = tv.tv_sec + tv.tv_usec/1000000.0;
	return sec;
}