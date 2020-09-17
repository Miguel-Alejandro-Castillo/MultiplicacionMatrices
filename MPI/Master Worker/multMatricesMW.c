#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>

#define USAGE_MSG "Use mpirun -np <processes> --hostfile <hostfile> %s <sizeMatrix> <sizeBlock> <threads> <print?> \n"
#define BASENAME_CSV "result_master-worker-parallel.csv"
#define HEADER_CSV "sizeMatrix,sizeBlock,threads,duration\n"
#define MIN(a, b) ((a)<(b) ? (a) : (b))

void cargarMatriz(double *M, int N);
void imprimirMatriz(double *M, int F, int C);
void imprimirMatrices(double *A, double *B, double * C, int N);
void guardarEjecucion(int N, int r, int t, double d);
void multBloques(double *A, double *B, double *C, int N, int p, int r);
void multBloquesOpenMP(double *A, double *B, double *C, int N, int p, int r, int t);
long get_integer_arg(int argc, char* argv[], int arg_index, long min_val, const char* description, const char* usage_msg, bool print_flag, int id, void (*fun) (void) );


int main(int argc, char **argv) {
	int nproc, myid, sizeMatrix, sizeBlock, threads;
	bool print;
	double *A, *B, *C, *BUFFER_A, *BUFFER_C, t1, t2;
	MPI_Comm comm = MPI_COMM_WORLD;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(comm, &nproc);
	MPI_Comm_rank(comm, &myid );
	
    /* argumentos de linea de comandos */
	sizeMatrix = get_integer_arg(argc, argv, 1, 1, "sizeMatrix", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
	sizeBlock = get_integer_arg(argc, argv, 2, 1, "sizeBlock", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
    threads = get_integer_arg(argc, argv, 3, 1, "threads", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
    print = get_integer_arg(argc, argv, 4, 0, "print", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
    int NxN = sizeMatrix * sizeMatrix;

    B = (double *)malloc(NxN*sizeof(double));
    BUFFER_A = (double *)malloc((NxN / nproc)*sizeof(double));
    BUFFER_C = (double *)malloc((NxN / nproc)*sizeof(double));

    if(myid == 0){
    	/*Se inicializan las matrices */
		printf("Inicializacion ...\n");
		A = (double *)malloc(NxN*sizeof(double));
    	C = (double *)malloc(NxN*sizeof(double));
    	cargarMatriz(A, sizeMatrix);
    	cargarMatriz(B, sizeMatrix);
    	printf("Ok!!\n\n");

		printf("Multiplicacion ...\n");
		t1 = MPI_Wtime();
    }

    /*  Se realiza la distribucion de las matrices  */
    MPI_Scatter(A, NxN /  nproc, MPI_DOUBLE, BUFFER_A, NxN / nproc, MPI_DOUBLE, 0, comm);
    MPI_Bcast(B, NxN, MPI_DOUBLE, 0, comm);	
    /* Se realiza la multiplicacion  */
    if( threads == 1 )
    	multBloques(BUFFER_A, B, BUFFER_C, sizeMatrix, nproc, sizeBlock);
    else
    	multBloquesOpenMP(BUFFER_A, B,BUFFER_C, sizeMatrix, nproc, sizeBlock, threads);
    /* Se envian y combinan los resultados */
    MPI_Gather( BUFFER_C, NxN / nproc, MPI_DOUBLE, C, NxN / nproc, MPI_DOUBLE, 0, comm);
    
    if(myid == 0) {
    	t2 = MPI_Wtime();
		printf("Ok!!\n\n");
    	if(print) {
    		imprimirMatrices(A, B, C, sizeMatrix);
		}
		printf("\nDuración total de la multiplicacion de matrices %4f segundos\n", t2-t1);
		guardarEjecucion(sizeMatrix, sizeBlock, threads, t2-t1);
		free(A);
		free(C);
	}
	free(B);
    free(BUFFER_A);
    free(BUFFER_C);

	MPI_Finalize(); 
	return 0;
}


void multBloques(double *A,double *B,double *C, int N, int p, int r){
  int kk, jj, i, j, k; 
  double temp;
  for(jj=0;jj<N;jj+= r){
        for(kk=0;kk<N;kk+= r){
                for(i=0;i<N/p;i++){
                        for(j = jj; j < MIN(N, jj+r); j++){
                            temp = 0.0;
                            for(k = kk; k <  MIN(N, kk+r); k++){
                                    temp += A[i*N+k] * B[k*N+j];
                            }
                            C[i*N+j] += temp;
                        }
                }
         }
   } 
}

void multBloquesOpenMP(double *A,double *B,double *C, int N, int p, int r, int t){
  int kk, jj, i, j, k; 
  double temp;
  omp_set_num_threads(t);
  #pragma omp parallel shared(A,B,C,N,p,r) private(kk,i,j,k,temp) 
  {
  	#pragma omp for schedule(static)
    for(jj=0;jj<N;jj+= r){
        for(kk=0;kk<N;kk+= r){
                for(i=0;i<N/p;i++){
                        for(j = jj; j < MIN(N, jj+r); j++){
                            temp = 0.0;
                            for(k = kk; k <  MIN(N, kk+r); k++){
                                    temp += A[i*N+k] * B[k*N+j];
                            }
                            C[i*N+j] += temp;
                        }
                }
         }
  	 } 
  }
}


void cargarMatriz(double *M, int N)
{
    int i, j, indice;
	for (i=0; i < N; i++){
		for (j=0; j < N; j++) {
            indice = i*N+j;
			M[indice] = indice;
		}
	}
}


void imprimirMatriz(double *M, int F, int C)
{
	int i, j = 0;
	for (i=0; i < F; i++) {
		printf("\n\t| ");
		for (j=0; j < C; j++)
			printf("%4f ", M[i*C+j]);
		printf("|");
	}
}

void imprimirMatrices(double *A, double *B, double * C, int N)
{
	printf("Imprimiendo Matriz A...\n");
	imprimirMatriz(A, N, N);
	printf("\n\n");
	printf("Imprimiendo Matriz B...\n");
	imprimirMatriz(B, N, N);
	printf("\n\n");
	printf("Imprimiendo Matriz C...\n");
	imprimirMatriz(C, N, N);
	printf("\n\n");
}


FILE* makeOutfile(char* basename, char *header) {
    FILE *outfile;
    if( access( basename, F_OK ) != -1 ) {
       outfile = fopen(basename, "a");
    } else {
       outfile = fopen(basename, "w");
       if(outfile != NULL)
            fprintf(outfile, "%s", header);
    }
    if (outfile == NULL)
        fprintf(stderr, "Cannot open outfile %s\n", basename);
    return outfile;
}

void guardarEjecucion(int N, int r, int t, double d){
	FILE* outfile = makeOutfile(BASENAME_CSV, HEADER_CSV);
    fprintf(outfile, "%d,%d,%d,%4f\n", N, r, t, d);
    fclose(outfile);
}


long get_integer_arg(int argc, char* argv[], int arg_index, 
        long min_val, const char* description, const char* usage_msg,
        bool print_flag, int id, void (*fun) (void) )
{
    if (arg_index >= argc) {
        if (print_flag && id == 0)
            fprintf(stderr, usage_msg, argv[0]);
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    char* end;
    long rval = strtol(argv[arg_index], &end, 10);
    if (*end != '\0') {
        if (print_flag && id == 0) {
            fprintf(stderr, "%s debe ser un número entero\n", description);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    if (rval < min_val) {
        if (print_flag && id == 0) {
            fprintf(stderr, "%s debe ser al menos %ld\n", description, min_val);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (fun != NULL) fun();
        exit(EXIT_FAILURE);
    }
    return rval;
}