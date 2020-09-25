/* Matrix multiplication, blocked version */
/* Enzo Rucci */

#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>

#define N 128
#define BS 32

void inicializarMatrix(double *m, int n);

/* Init square matrix */
void initvalmat(double *mat, int n, int transpose);

/* Multiply square matrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs);

/* Time in seconds from some point in the past */
double dwalltime();

void imprimeMatriz(double *M, int n);

/************** MAIN *************/
int main(int argc, char *argv[])
{
	double *a, *b, *c;
	int n = N, bs = BS, i, j, print;
	double time, timetick;

	/* Check command line parameters */
	if ( (argc != 4) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((print = atoi(argv[3])) < 0) || ((n % bs) != 0))
	{
		printf("\nUsage: %s n bs print? \n  n: matrix order (nxn X nxn)\n  bs: block size (n should be multiple of bs)\n  print?: print result \n", argv[0]);
		exit(1);
	}

	/* Allocate memory */
	a = (double *) malloc(n * n * sizeof(double));
	b = (double *) malloc(n * n * sizeof(double));
	c = (double *) malloc(n * n * sizeof(double));

	/* Init matrix operands */
	initvalmat(a, n, 0);
	initvalmat(b, n, 1);
	/* Init matrix c, just in case */
	inicializarMatrix(c, n);

	timetick = dwalltime();

	matmulblks(a, b, c, n, bs);

	time =  dwalltime() - timetick;

	printf("Multiplicación de matrices de %dx%d (bloque = %d) tomó %f segundos\n", n, n, bs, time);

	if (print == 1) {
		printf("\n\n  A: \n" );
		imprimeMatriz(a, n);

		printf("\n\n B: \n" );
		imprimeMatriz(b, n);

		printf("\n\n  C: \n" );
		imprimeMatriz(c, n);
	}

	free(a);
	free(b);
	free(c);

	return 0;
}

/*****************************************************************/

void inicializarMatrix(double *m, int n) {
	int i, j;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			m[i * n + j] = 0.0;
		}
	}
}

/* Init square matrix */
void initvalmat(double *mat, int n, int transpose)
{
	int i, j;      /* Indexes */

	if (transpose == 0) {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				mat[i * n + j] = i * n + j;
			}
		}
	} else {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				mat[j * n + i] = i * n + j;
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
				dispB = (jj * numBlocks + kk) * blockElems;
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