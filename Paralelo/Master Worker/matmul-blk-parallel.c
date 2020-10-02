/********************************************************************************************************
 *
 * Archivo: matmul-blk-parallel.c
 *
 * Multiplicacion de matrices cuadradas por bloques utilizando la Interfaz de Paso de Mensajes MPI
 * + OpenMP, modelo Master Worker
 *
 * Para compilar:
 * mpicc -fopenmp -o <executable> matmul-blk-parallel.c
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

#define USAGE_MSG "Use mpirun -np <processes> --hostfile <hostfile> %s <sizeMatrix> <sizeBlock> <threads>\n"
#define BASENAME_CSV "result_matmul-blk-parallel.csv"
#define HEADER_CSV "sizeMatrix,sizeBlock,threads,time\n"

void imprimirMatrices(double *A, double *B, double * C, int N);
void guardarEjecucion(int sizeMatrix, int sizeBlock, int threads, double time);
void mulblks(double *a, double *b, double *c, int sizeMatrix, int sizeBlock, int nproc);
void mulblksOpenMP(double *a, double *b, double *c, int sizeMatrix, int sizeBlock, int nproc, int threads);
void matmulblks(double *a, double *b, double *c, int n, int bs, int iSubMatrix, int jSubMatrix, int kSubMatrix, int sizeSubMatrix);
void initvalmat(double *mat, int n, double val, int transpose);
long get_integer_arg(int argc, char* argv[], int arg_index, long min_val, const char* description, const char* usage_msg, bool print_flag, int id, void (*fun) (void) );
void imprimirMatriz(double *M, int N);
void inicializarMatrix(double *S, int sizeMatrix);

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
  B = (double *)malloc(NxN * sizeof(double));
  BUFFER_A = (double *)malloc((NxN / nproc) * sizeof(double));
  BUFFER_C = (double *)malloc((NxN / nproc) * sizeof(double));
  /* Se inicializa la submatrix buffer_C con ceros */
  inicializarMatrix(BUFFER_C, NxN / nproc);

  if (myid == 0) {
    /* Se inicializan las matrices A y B */
    printf("\nInicializacion ...\n");
    A = (double *)malloc(NxN * sizeof(double));
    C = (double *)malloc(NxN * sizeof(double));
    initvalmat(A, sizeMatrix, 1.0, 0);
    /* Se inicializa la matrix B transpuesta */
    initvalmat(B, sizeMatrix, 1.0, 1);
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
    mulblks(BUFFER_A, B, BUFFER_C, sizeMatrix, sizeBlock, nproc);
  else
    mulblksOpenMP(BUFFER_A, B, BUFFER_C, sizeMatrix, sizeBlock, nproc, threads);

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
    printf("\nDuración total de la multiplicacion de matrices %4f segundos\n", t2 - t1);
    /* La ejecucion se almacena en un archivo .csv */
    guardarEjecucion(sizeMatrix, sizeBlock, threads, t2 - t1);
    free(A);
    free(C);
  }
  free(B);
  free(BUFFER_A);
  free(BUFFER_C);

  MPI_Finalize();
  return 0;
}

void mulblks(double *a, double *b, double *c, int sizeMatrix, int sizeBlock, int nproc) {
  int sizeSubMatrix = sizeMatrix / nproc;
  int cantSubMatrix = sizeMatrix / sizeSubMatrix;
  int auxSizeBlock = sizeBlock > sizeSubMatrix ? sizeSubMatrix : sizeBlock;
  int j, k;

  //Se multiplica una fila de cantSubMatrix submatrices de tamaño sizeSubMatrix cada una
  for ( j = 0; j < cantSubMatrix; j++) {
    for (k = 0; k < cantSubMatrix; k++) {
      //Se multiplican 2 submatrices cuadradas de a y b, el resultado se guarda en una submatrix de c
      matmulblks(a, b, c, sizeMatrix, auxSizeBlock, 0, j, k, sizeSubMatrix);
    }
  }
}

void mulblksOpenMP(double *a, double *b, double *c, int sizeMatrix, int sizeBlock, int nproc, int threads) {
  int sizeSubMatrix = sizeMatrix / nproc;
  int cantSubMatrix = sizeMatrix / sizeSubMatrix;
  int auxSizeBlock = sizeBlock > sizeSubMatrix ? sizeSubMatrix : sizeBlock;
  int j, k;


  omp_set_num_threads(threads);
  #pragma omp parallel shared(a,b,c,sizeMatrix,auxSizeBlock,sizeSubMatrix,cantSubMatrix) private(j,k)
  {
    #pragma omp for schedule(static)
    //Se multiplica una fila de cantSubMatrix submatrices de tamaño sizeSubMatrix cada una
    for ( j = 0; j < cantSubMatrix; j++) {
      for (k = 0; k < cantSubMatrix; k++) {
        //Se multiplican 2 submatrices cuadradas de a y b, el resultado se guarda en una submatrix de c
        matmulblks(a, b, c, sizeMatrix, auxSizeBlock, 0, j, k, sizeSubMatrix);
      }
    }
  }
}

/* Multiply square submatrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs, int iSubMatrix, int jSubMatrix, int kSubMatrix, int sizeSubMatrix)
{
  int i, j, k, ii, jj, kk, disp, dispA, dispB, dispC;
  int iAuxSubMatrix = iSubMatrix * sizeSubMatrix;
  int jAuxSubMatrix = jSubMatrix * sizeSubMatrix;
  int kAuxSubMatrix = kSubMatrix * sizeSubMatrix;
  int auxDispC = iAuxSubMatrix * n + jAuxSubMatrix;
  int auxDispA = iAuxSubMatrix * n + kAuxSubMatrix;
  int auxDispB = kAuxSubMatrix * n + jAuxSubMatrix;

  int numBlocks = sizeSubMatrix / bs;
  int blockElems = bs * bs;

  for (ii = 0; ii < numBlocks; ii++) {
    for (jj = 0; jj < numBlocks; jj++) {
      dispC = (ii * numBlocks + jj) * blockElems;
      for (kk = 0; kk < numBlocks; kk++) {
        dispA = (ii * numBlocks + kk) * blockElems;
        dispB = (kk * numBlocks + jj) * blockElems;
        for (i = 0; i < bs; i++) {
          for (j = 0; j < bs; j++) {
            disp = dispC + (i * bs + j) + auxDispC;
            for (k = 0; k < bs; k++) {
              c[disp] += a[dispA + (i * bs + k)  + auxDispA] * b[dispB + (j * bs + k) + auxDispB];
            }
          }
        }
      }
    }
  }
}

/* Init square matrix with a specific value */
void initvalmat(double * mat, int n, double val, int transpose)
{
  int i, j;      /* Indexes */

  if (transpose == 0) {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        mat[i * n + j] = val;
      }
    }
  } else {
    for (i = 0; i < n; i++) {
      for (j = 0; j < n; j++) {
        mat[j * n + i] = val;
      }
    }
  }
}

void inicializarMatrix(double * S, int sizeMatrix) {
  int i;
  for (i = 0; i < sizeMatrix ; i++) {
    S[i] = 0.0;
  };
}

void imprimirMatriz(double * M, int N)
{
  int i, j;
  printf("Contenido de la matriz:  \n");
  for (i = 0 ; i < N; i++) {
    for (j = 0; j < N; j++)
      printf(" %.2f " , M[i * N + j]);
    printf("\n ");
  }
  printf("\n");
}

void imprimirMatrices(double * A, double * B, double * C, int N)
{
  printf("\n  A: \n");
  imprimirMatriz(A, N);
  printf("\n  B: \n");
  imprimirMatriz(B, N);
  printf("\n  C: \n");
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

void guardarEjecucion(int sizeMatrix, int sizeBlock, int threads, double time) {
  FILE* outfile = makeOutfile(BASENAME_CSV, HEADER_CSV);
  fprintf(outfile, "%d,%d,%d,%4f\n", sizeMatrix, sizeBlock, threads, time);
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