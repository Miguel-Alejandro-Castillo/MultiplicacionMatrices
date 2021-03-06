/********************************************************************************************************
 *
 * Archivo: matmul-block-parallel.c
 *
 * Multiplicacion de matrices cuadradas por bloques utilizando la Interfaz de Paso de Mensajes MPI
 * + OpenMP, modelo Master Worker
 *
 * Para compilar:
 * mpicc matmul-block-parallel.c -fopenmp -O3 -o matmul-block-parallel
 *
 * Para ejecutar:
 * mpirun -np <processes> --hostfile <hostfile> <executable> <sizeMatrix> <sizeBlock> <threads>
 *******************************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <unistd.h>
#include <mpi.h>
#include <omp.h>

#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <semaphore.h>

#define USAGE_MSG "Use mpirun -np <processes> --hostfile <hostfile> %s <sizeMatrix> <sizeBlock> <threads>\n"
#define BASENAME_MPI_CSV "result_matmul-block-mpi.csv"
#define BASENAME_MPI_OPENMP_CSV "result_matmul-block-mpi+openmp.csv"
#define HEADER_CSV "processes,threads,sizeMatrix,sizeBlock,time\n"
#define min(a, b) ((a)<(b) ? (a) : (b))
#define SEM_NAME "/mysem"

void imprimirMatrices(double *A, double *B, double * C, int N);
void guardarEjecucion(int nproc, int sizeMatrix, int sizeBlock, int threads, double time);
void matmulblks(double *a, double *b, double *c, int n, int bs, int p);
void matmulblksOpenMP(double *a, double *b, double *c, int n, int bs, int p, int t);
long get_integer_arg(int argc, char* argv[], int arg_index, long min_val, const char* description, const char* usage_msg, bool print_flag, int id, void (*fun) (void) );
void imprimirMatriz(double *M, int N);
void initializeMatrix(double *S, int sizeMatrix);
void initmat(double *m, int n, int transpose);


int main(int argc, char **argv) {
  int nproc, myid, sizeMatrix, sizeBlock, threads;
  double *A, *B, *C, *BUFFER_A, *BUFFER_C, t1, t2;
  MPI_Comm comm = MPI_COMM_WORLD;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_size(comm, &nproc);
  MPI_Comm_rank(comm, &myid );

  /* Se obtienen los argumentos de linea de comandos */
  sizeMatrix = get_integer_arg(argc, argv, 1, 1, "sizeMatrix", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
  sizeBlock = get_integer_arg(argc, argv, 2, 1, "sizeBlock", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);
  threads = get_integer_arg(argc, argv, 3, 1, "threads", USAGE_MSG, true, myid, (void (*)(void))MPI_Finalize);

  /* Check command line parameters */
  if (((sizeMatrix % sizeBlock) != 0) || ((sizeMatrix % nproc) != 0))
  {
    if (myid == 0) {
      printf("%s", USAGE_MSG);
      printf("processes: number of processes (sizeMatrix should be multiple of processes)\n");
      printf("sizeBlock: block size (sizeMatrix should be multiple of sizeBlock)\n");
    }
    MPI_Finalize();
    exit(EXIT_FAILURE);
  }

  int NxN = sizeMatrix * sizeMatrix;
  BUFFER_C = (double *)malloc((NxN / nproc) * sizeof(double));
  BUFFER_A = (double *)malloc((NxN / nproc) * sizeof(double));
  B = (double *)malloc(NxN * sizeof(double));

  /* Se inicializa la submatrix buffer_C con ceros */
  initializeMatrix(BUFFER_C, NxN / nproc);

  if (myid == 0) {
    /* Se inicializan las matrices A y B */
    printf("\nInicializacion ...\n");
    A = (double *)malloc(NxN * sizeof(double));
    C = (double *)malloc(NxN * sizeof(double));
    /* Se inicializa la matrix A */
    initmat(A, sizeMatrix, 0);
    /* Se inicializa la matrix B transpuesta */
    initmat(B, sizeMatrix, 1);
    printf("Ok!!\n\n");

    printf("Multiplicacion ...\n");
    /* Inicia el conteo */
    t1 = MPI_Wtime();
  }

  /*  Se realiza la distribucion de las matrices  */
  MPI_Scatter(A, NxN /  nproc, MPI_DOUBLE, BUFFER_A, NxN / nproc, MPI_DOUBLE, 0, comm);
  MPI_Bcast(B, NxN, MPI_DOUBLE, 0, comm);

  /* Se realiza la multiplicacion, si la cantidad de threads es mayor a 1 se hace uso de OpenMP para el computo */
  if (threads == 1)
    matmulblks(BUFFER_A, B, BUFFER_C, sizeMatrix, sizeBlock, nproc);
  else
    matmulblksOpenMP(BUFFER_A, B, BUFFER_C, sizeMatrix, sizeBlock, nproc, threads);

  /* Se envian(submatrices buffer_c) y combinan los resultados en la matriz c */
  MPI_Gather(BUFFER_C, NxN / nproc, MPI_DOUBLE, C, NxN / nproc, MPI_DOUBLE, 0, comm);

  if (myid == 0) {
    /* Finaliza el conteo */
    t2 = MPI_Wtime();
    printf("Ok!!\n\n");
    if (sizeMatrix <= 16) {
      //Se imprime en pantalla las matrices a, b y c
      imprimirMatrices(A, B, C, sizeMatrix);
    }
    printf("\nDuración total de la multiplicacion de matrices %4f segundos\n\n", t2 - t1);
    /* La ejecucion se almacena en un archivo .csv */

    sem_t *sem;
    sem = sem_open(SEM_NAME, O_CREAT|O_EXCL, S_IRUSR|S_IWUSR, 1);
    if (sem ==  SEM_FAILED && errno == EEXIST) {
      //printf("semaphore  appears  to  exist  already !\n");
      sem = sem_open(SEM_NAME , 0);
    }
    sem_wait(sem);
    guardarEjecucion(nproc, sizeMatrix, sizeBlock, threads, t2 - t1);
    sem_post(sem);
    sem_close(sem);

    free(A);
    free(C);
  }
  free(B);
  free(BUFFER_A);
  free(BUFFER_C);

  MPI_Finalize();
  return 0;
}

/* Multiply square matrices, blocked version */
/* Version optimizada con multiplicacion de A * B'(transpuesta) y uso de variables temporales */
void matmulblks(double *a, double *b, double *c, int n, int bs, int p)
{
  int i0, j0, k0, i, j, k;
  int disp, dispB;
  double temp;
  for (i0 = 0; i0 < n / p; i0 += bs) {
    for (j0 = 0; j0 < n; j0 += bs) {
      for (k0 = 0; k0 < n; k0 += bs) {
        for (i = i0; i < min(i0 + bs, n / p); i++) {
          disp = i * n;
          for (j = j0; j < min(j0 + bs, n); j++) {
            temp = 0.0;
            dispB = j * n;
            for (k = k0; k < min(k0 + bs, n); k++) {
              temp += a[disp + k] * b[dispB + k];
            }
            c[disp + j] += temp;
          }
        }
      }
    }
  }

}

void matmulblksOpenMP(double *a, double *b, double *c, int n, int bs, int p, int t)
{
  int i0, j0, k0, i, j, k;
  int disp, dispB;
  double temp;
  omp_set_num_threads(t);
  #pragma omp parallel shared(c,a,b,n,bs,p) private(i0,j0,k0,i,j,k,disp,dispB,temp)
  {
    #pragma omp for schedule(static)
    for (i0 = 0; i0 < n / p; i0 += bs) {
      for (j0 = 0; j0 < n; j0 += bs) {
        for (k0 = 0; k0 < n; k0 += bs) {
          for (i = i0; i < min(i0 + bs, n / p); i++) {
            disp = i * n;
            for (j = j0; j < min(j0 + bs, n); j++) {
              temp = 0.0;
              dispB = j * n;
              for (k = k0; k < min(k0 + bs, n); k++) {
                temp += a[disp + k] * b[dispB + k];
              }
              c[disp + j] += temp;
            }
          }
        }
      }
    }
  }
}

/* Init square matrix with a specific value */
void initmat(double *m, int n, int transpose) {
  int i, j;

  if (transpose == 0) {
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        m[i * n + j] = i * n + j;
      }
    }
  } else {
    for (i = 0; i < n; i++)
    {
      for (j = 0; j < n; j++)
      {
        m[j * n + i] = i * n + j;
      }
    }
  }
}

void initializeMatrix(double * S, int sizeMatrix) {
  int i;
  for (i = 0; i < sizeMatrix ; i++) {
    S[i] = 0.0;
  };
}

void imprimirMatriz(double * M, int N)
{
  int i, j;
  printf("Contenido de la matriz:  \n ");
  for (i = 0 ; i < N; i++) {
    for (j = 0; j < N; j++)
      printf(" %.2f " , M[i * N + j]);
    printf("\n ");
  }
  printf("\n");
}

void imprimirMatrices(double * A, double * B, double * C, int N)
{
  printf("\n A: \n");
  imprimirMatriz(A, N);
  printf("\n B: \n");
  imprimirMatriz(B, N);
  printf("\n C: \n");
  imprimirMatriz(C, N);
}

FILE* makeOutfile(char* basename, char *header) {
  FILE *outfile;
  if ( access( basename, F_OK ) != -1 ) {
    outfile = fopen(basename, "a");
  } else {
    outfile = fopen(basename, "w");
    if (outfile != NULL)
    	fprintf(outfile, "%s", header);
  }
  if (outfile == NULL)
    fprintf(stderr, "Cannot open outfile %s\n", basename);
  return outfile;
}

void guardarEjecucion(int nproc, int sizeMatrix, int sizeBlock, int threads, double time) {
  FILE* outfile;
  if(threads == 1) {
	outfile = makeOutfile(BASENAME_MPI_CSV, HEADER_CSV);
  } else {
	outfile = makeOutfile(BASENAME_MPI_OPENMP_CSV, HEADER_CSV);
  }
  if(outfile != NULL){
	fprintf(outfile, "%d,%d,%d,%d,%4f\n", nproc, threads > 1 ? nproc * threads : threads, sizeMatrix, sizeBlock, time);
	fclose(outfile);
  }
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
