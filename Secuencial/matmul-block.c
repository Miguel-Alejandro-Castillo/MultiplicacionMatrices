/* Matrix multiplication, blocked version */

#include<stdio.h>
#include<stdlib.h>
#include<sys/time.h>
#include <unistd.h>

#include <fcntl.h>
#include <errno.h>
#include <sys/stat.h>
#include <semaphore.h>

#define BASENAME_CSV "result_matmul-block-secuencial.csv"
#define HEADER_CSV "sizeMatrix,sizeBlock,time\n"
#define min(a, b) ((a)<(b) ? (a) : (b))
#define SEM_NAME "/semsecuencial"

/* Init square matrix  */
void initvalmat(double *m, int n, double val, int transpose);

void initmat(double *m, int n, int transpose);

/* Multiply square matrices, blocked version */
void matmulblks(double *a, double *b, double *c, int n, int bs);

/* Time in seconds from some point in the past */
double dwalltime();

void printMatrix(double *m, int n);

void saveExecution(int sizeMatrix, int sizeBlock, double time);

/************** MAIN *************/
int main(int argc, char *argv[])
{
	double *a, *b, *c;
	int n, bs;
	double time, timetick;

	/* Check command line parameters */
	if ( (argc != 3) || ((n = atoi(argv[1])) <= 0) || ((bs = atoi(argv[2])) <= 0) || ((n % bs) != 0))
	{
		printf("\nUsage: %s n bs \n  \tn: matrix order (nxn X nxn)\n  \tbs: block size (n should be multiple of bs)\n\n", argv[0]);
		exit(1);
	}

	/* Allocate memory */
	a = (double *) malloc(n * n * sizeof(double));
	b = (double *) malloc(n * n * sizeof(double));
	c = (double *) malloc(n * n * sizeof(double));

	/* Init matrix a */
	initmat(a, n, 0);
	/* Init matrix b transpose */
	initmat(b, n, 1);
	/* Init matrix c */
	initvalmat(c, n, 0.0, 0);

	timetick = dwalltime();

	matmulblks(a, b, c, n, bs);

	time =  dwalltime() - timetick;

	printf("Multiplicación de matrices de %dx%d (bloque = %d) tomó %4f segundos\n", n, n, bs, time);

	if (n <= 16) {
		printf("\n A: \n" );
		printMatrix(a, n);

		printf("\n B: \n" );
		printMatrix(b, n);

		printf("\n C: \n" );
		printMatrix(c, n);
	}

    sem_t *sem;
    sem = sem_open(SEM_NAME, O_CREAT|O_EXCL, S_IRUSR|S_IWUSR, 1);
    if (sem ==  SEM_FAILED && errno == EEXIST) {
      //printf("semaphore  appears  to  exist  already !\n");
      sem = sem_open(SEM_NAME , 0);
    }
    sem_wait(sem);
    saveExecution(n, bs, time);
    sem_post(sem);
    sem_close(sem);

	free(a);
	free(b);
	free(c);

	return 0;
}

/*****************************************************************/

/* Multiply square matrices, blocked version */
/* Version optimizada con multiplicacion de A * B'(transpuesta) y uso de variables temporales */
void matmulblks(double *a, double *b, double *c, int n, int bs)
{
	int i0, j0, k0, i, j, k;
	int dispA, dispB, dispC;
	double temp;
	for (i0 = 0; i0 < n; i0 += bs) {
		for (j0 = 0; j0 < n; j0 += bs) {
			for (k0 = 0; k0 < n; k0 += bs) {
				for (i = i0; i < min(i0 + bs, n); i++) {
					dispC = i * n;
					dispA = i * n;
					for (j = j0; j < min(j0 + bs, n); j++) {
						temp = 0.0;
						dispB = j * n;
						for (k = k0; k < min(k0 + bs, n); k++) {
							temp += a[dispA + k] * b[dispB + k];
						}
						c[dispC + j] += temp;
					}
				}
			}
		}
	}
}

/* Init square matrix with a specific value */
void initvalmat(double *m, int n, double val, int transpose)
{
	int i, j;

	if (transpose == 0) {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				m[i * n + j] = val;
			}
		}
	} else {
		for (i = 0; i < n; i++)
		{
			for (j = 0; j < n; j++)
			{
				m[j * n + i] = val;
			}
		}
	}
}

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

void printMatrix(double *m, int n) {
	int i, j;
	printf("Contenido de la matriz:  \n ");
	for (i = 0 ; i < n; i++) {
		for (j = 0; j < n; j++)
			printf(" %.2f " , m[i * n + j]);
		printf("\n ");
	}
	printf("\n");
}

double dwalltime()
{
	double sec;
	struct timeval tv;
	gettimeofday(&tv, NULL);
	sec = tv.tv_sec + tv.tv_usec / 1000000.0;
	return sec;
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

void saveExecution(int sizeMatrix, int sizeBlock, double time) {
	FILE* outfile = makeOutfile(BASENAME_CSV, HEADER_CSV);
	fprintf(outfile, "%d,%d,%4f\n", sizeMatrix, sizeBlock, time);
	//fprintf(outfile, "\"(%d,%d)\",%4f\n", sizeMatrix, sizeBlock, time);
	fclose(outfile);
}
