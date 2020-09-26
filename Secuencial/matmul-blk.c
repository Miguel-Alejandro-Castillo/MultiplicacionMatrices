/* Matrix multiplication, blocked version */
/* Enzo Rucci */

#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include <unistd.h>

#define N 128
#define BS 32
#define BASENAME_CSV "result_matmul-blk-secuencial.csv"
#define HEADER_CSV "sizeMatrix,sizeBlock,time\n"

/* Init square matrix with a specific value */
void initvalmat(double *mat, int n, double val, int transpose);

/* Multiply square matrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs);

/* Time in seconds from some point in the past */
double dwalltime();

void imprimeMatriz(double *M, int n);

void guardarEjecucion(int sizeMatrix, int sizeBlock, double time);

/************** MAIN *************/
int main(int argc, char *argv[])
{
	double *a, *b, *c;
	int check = 1, n = N, bs = BS, i, j;
	double time, timetick;

	/* Check command line parameters */
	if ( (argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
	{
		printf("\nUsage: %s n bs \n  n: matrix order (nxn X nxn)\n  bs: block size (n should be multiple of bs)\n", argv[0]);
		exit(1);
	}

	/* Allocate memory */
	a = (double *) malloc(n * n * sizeof(double));
	b = (double *) malloc(n * n * sizeof(double));
	c = (double *) malloc(n * n * sizeof(double));

	/* Init matrix operands */
	initvalmat(a, n, 1.0, 0);
	initvalmat(b, n, 1.0, 1);
	/* Init matrix c, just in case */
	initvalmat(c, n, 0.0, 0);

	timetick = dwalltime();

	matmulblks(a, b, c, n, bs);

	time =  dwalltime() - timetick;

	printf("Multiplicación de matrices de %dx%d (bloque = %d) tomó %4f segundos\n", n, n, bs, time);

	if (n <= 16) {
		printf("\n  A: \n" );
		imprimeMatriz(a, n);

		printf("\n B: \n" );
		imprimeMatriz(b, n);

		printf("\n  C: \n" );
		imprimeMatriz(c, n);
	}

	// Check results
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			check = check && (c[i * n + j] == n);
		}
	}

	// Print results
	if (check)
		printf("Multiplicacion de matrices resultado correcto\n");
	else
		printf("Multiplicacion de matrices resultado erroneo\n");

	guardarEjecucion(n, bs, time);

	free(a);
	free(b);
	free(c);

	return 0;
}

/*****************************************************************/

/* Init square matrix with a specific value */
void initvalmat(double *mat, int n, double val, int transpose)
{
	int i, j;      /* Indexes */

	if (transpose == 0) {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				mat[i * n + j] = val;
			}
		}
	} else {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				mat[j * n + i] = val;
			}
		}
	}
}

/*****************************************************************/

/* Multiply square matrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs)
{
	int i, j, k, ii, jj, kk, disp, dispA, dispB, dispC;
	int numBlocks = n / bs;
	int blockElems = bs * bs;

	for (ii = 0; ii < numBlocks; ii++) {
		for (jj = 0; jj < numBlocks; jj++) {
			dispC = (ii * numBlocks + jj) * blockElems;
			for (kk = 0; kk < numBlocks; kk++) {
				dispA = (ii * numBlocks + kk) * blockElems;
				dispB = (kk * numBlocks + jj) * blockElems;
				for (i = 0; i < bs; i++) {
					for (j = 0; j < bs; j++) {
						disp = dispC + i * bs + j;
						for (k = 0; k < bs; k++) {
							c[disp] += a[dispA + i * bs + k] * b[dispB + j * bs + k];
						}
					}
				}
			}
		}
	}
}

double dwalltime()
{
	double sec;
	struct timeval tv;

	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
}


void imprimeMatriz(double *M, int n) {
	int i, j;
	printf("Contenido de la matriz:  \n ");
	for (i = 0 ; i < n; i++) {
		for (j = 0; j < n; j++)
			printf(" %.2f " , M[i * n + j]);
		printf("\n ");
	}
	printf("\n");
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

void guardarEjecucion(int sizeMatrix, int sizeBlock, double time) {
	FILE* outfile = makeOutfile(BASENAME_CSV, HEADER_CSV);
	fprintf(outfile, "%d,%d,%4f\n", sizeMatrix, sizeBlock, time);
	fclose(outfile);
}