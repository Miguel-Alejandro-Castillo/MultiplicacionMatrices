#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include <math.h>

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

int main (int argc, char *argv[]){

 double *A; // Matriz A
 double *B; // Matriz B
 double *C; // Matriz C
 double timetick;

if (argc < 4){
  printf("\n Faltan parametros ");
  //printf("\n 1. Cantidad de bloques por dimension ");
  printf("\n 1. Tamaño de la matriz ");
  //printf("\n 2. Dimension de cada bloque ");
  printf("\n 2. Tamaño de cada bloque por dimensión ");
  printf("\n 3. 0 no imprime el resultado, 1 imprime el resultado ");
  printf("\n\n");
  return 0;
}

 int N = atoi(argv[1]);
 int r = atoi(argv[2]);
 int imprimir = atoi(argv[3]);

 int sizeMatrix = N*N; //cantidad total de datos matriz
 int i;

 A= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para A
 B= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para B
 C= (double *)malloc(sizeMatrix*sizeof(double)); //aloca memoria para C

 crearMatriz(A, sizeMatrix);	//Inicializa A 
 crearMatriz(B, sizeMatrix);  //Inicializa B
 inicializarMatrix(C, sizeMatrix); 

 timetick = dwalltime();
 productoBloques(A,B, C, N, r);
 printf("Tiempo en segundos %f \n", dwalltime() - timetick);

 if (imprimir == 1){
    printf("\n\n  A (como esta almacenada): \n" );
    imprimeVector(A, sizeMatrix);

    printf("\n\n  B (como esta almacenada): \n" );
    imprimeVector(B, sizeMatrix);

    printf("\n\n  C (como esta almacenada): \n" );
    imprimeVector(C, sizeMatrix);

    printf("\n\n  A: \n" );
    imprimeMatriz2(A, N);

    printf("\n\n B: \n" );
    imprimeMatriz2(B, N);

    printf("\n\n  C: \n" );
    imprimeMatriz2(C, N);
 }

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
  for(jj=0;jj<N;jj+= r){
        for(kk=0;kk<N;kk+= r){
                for(i=0;i<N;i++){
                        for(j = jj; j < ((jj+r)>N?N:(jj+r)); j++){
                                double temp = 0.0;
                                for(k = kk; k < ((kk+r)>N?N:(kk+r)); k++){
                                        temp += A[i*N+k]*B[k*N+j];
                                }
                                C[i*N+j] += temp;
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