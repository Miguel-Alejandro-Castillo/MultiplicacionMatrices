
/*********************************************************************************************
 * Archivo: 2D_ciclico.c
 *
 * Multiplicacion de matrices utilizando la Interfaz de Paso de Mensajes MPI
 * y el particionamiento 2-D ciclico.
 *
 * Para compilar:
 * mpicc -o <nombre executable> 2d_ciclico.c
 *
 * Para ejecutar
 * mpirun -np <procesos> <nombre ejecutable> <sizeMatrix> <sizeBlock> <print?>
 **********************************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <mpi.h>
#include <stdbool.h>

struct indices    /* Tipo de dato para representar */
{	int i, j;     /* los ndices de una matriz */
};
typedef struct indices indice;

indice getIndice(int tarea, int subnum);
void imprimirMatriz(double *M, int F, int C);
void inicializarMatrices(double *A, double *B, double *C, int N);
void cargarBuffers(double *A, double *B, double *BUFFER_A, double *BUFFER_B, int N, int r, indice indi);
void inicializarBufferC(double *BUFFER_C, int r);
void asignarBUFFEResultado(double *C, double *BUFFER_C, int N, int r, indice indi);
void multiplicarSubMatriz(double *BUFFER_A, double *BUFFER_B, double *BUFFER_C, int N, int r);
long get_integer_arg(int argc, char* argv[], int arg_index, long min_val, const char* description, const char* usage_msg, bool print_flag, void (*mpi_finalize) (void) );
void error_exit(char* msg);

int main(int argc, char *argv[])
{
	int N = 4;   /* tamanho de la matriz */
	int r = 2;   /* tamanho de la submatriz (o bloques) */
	bool print; /* 1 imprime resultado, 0 no imprime resultado */
	int subnum;	 /* numero de sub-matrices en una fila/columna de la matriz */
	int nproc;	 /* numero de nodos MPI */
	int myid; 	 /* mi propio rank */
	double *BUFFER_A;      /* buffer rxN utilizado para la distribucion de sub-matrices de A
	                       * necesarias para calcular un determinado bloque de la matriz C */
	double *BUFFER_B;	  /* buffer  Nxr utilizado para la distribucion de sub-matrices de B
	                       * necesarias para calcular un determinado bloque de la matriz C */
	double *BUFFER_C;	  /* buffer rxr utilizado para almacenar el resultado de la multiplicacion
	                       * del bloque correspondiente a C */

	int i, j, k;
	double t1, t2;

    MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);
	MPI_Comm_size(MPI_COMM_WORLD, &nproc);

	/* command-line arguments */
    char *usage_msg = "usage is %s sizeMatrix sizeBlock print? \n";
	N = get_integer_arg(argc, argv, 1, 1, "sizeMatrix", usage_msg, true, (void (*)(void))MPI_Finalize);
	r = get_integer_arg(argc, argv, 2, 1, "sizeBlock", usage_msg, true, (void (*)(void))MPI_Finalize);
    print = get_integer_arg(argc, argv, 3, 0, "print", usage_msg, true, (void (*)(void))MPI_Finalize);

	subnum = N / r;
	if (subnum * r != N ) {
		subnum++;
	}

    if (nproc > r*r)
        nproc = r*r; /* Maxima cantidad de nodos necesarios */

	if ( myid == 0){

		double *A, *B, *C;
		int id_destino = 0;

        /*Se inicializan las matrices */
		printf("Inicializacion ...\n");
		A = (double *)malloc(N*N*sizeof(double));
   		B = (double *)malloc(N*N*sizeof(double));
    	C = (double *)malloc(N*N*sizeof(double));
    	BUFFER_A = (double *)malloc(N*N*sizeof(double));
   		BUFFER_B = (double *)malloc(N*N*sizeof(double));
    	BUFFER_C = (double *)malloc(N*N*sizeof(double));
		inicializarMatrices(A, B, C, N);
		printf("Ok!!\n\n");

        /* Se realiza la multiplicacion de matrices */
		printf("Multiplicacion ...\n");
		t1 = MPI_Wtime();

		int tarea; /* Cantidad total de tareas disponibles */
		for (tarea = 0; tarea < subnum*subnum; tarea++ ){

            /* Se carga los buffers con los datos necesarios para cada terea */
			indice indi = getIndice(tarea, subnum);

			cargarBuffers(A, B, BUFFER_A, BUFFER_B, N, r, indi);

            id_destino = tarea % nproc;

            if ( id_destino != 0 ){
                /* Enviar las sub-matrices a multiplicar */
                MPI_Send(BUFFER_A, N*r, MPI_DOUBLE, id_destino, 1, MPI_COMM_WORLD);
                MPI_Send(BUFFER_B, N*r, MPI_DOUBLE, id_destino, 2, MPI_COMM_WORLD);
            }
            else{
                /* Se multiplican localmente las sub-matrices */
                multiplicarSubMatriz(BUFFER_A, BUFFER_B, BUFFER_C, N, r);
                asignarBUFFEResultado(C, BUFFER_C, N, r, indi);
            }

            if ( id_destino == nproc-1 ){
                /* Se reciben todas las tareas distribuidas hasta el momento */
                int tarea_recibida = tarea - (nproc-1);
                int proc;
                for (proc = 0; proc <= nproc-1; proc++){
                    id_destino = tarea_recibida % nproc;
                    if (id_destino != 0){
                        /* Si la tarea no corresponde al root (id=0) entonces se recibe */
                        MPI_Recv(BUFFER_C, r*r, MPI_DOUBLE, id_destino,3, MPI_COMM_WORLD, &status);
                        indi = getIndice(tarea_recibida, subnum);
                        asignarBUFFEResultado(C, BUFFER_C, N, r, indi);
                    }
                    tarea_recibida++;
                }
            }
		}

        if ( (nproc < r*r) && (((r*r) % nproc) != 0) ){
            /* Se reciben las ultimas tareas pendientes */

            int tareas_pendientes = r*r % nproc; // cuantas tareas aun estan pendientes??
            int tarea_recibida = r*r - tareas_pendientes; // a partir de que tarea??
            int x;
            for (x = 0; x < tareas_pendientes; x++){
                id_destino = tarea_recibida % nproc;
                if (id_destino != 0){
                    /* Si la tarea no corresponde al root (id=0) entonces se recibe */
                    MPI_Recv(BUFFER_C, r*r, MPI_DOUBLE, id_destino,3, MPI_COMM_WORLD, &status);
                    indice indi = getIndice(tarea_recibida, subnum);
                    asignarBUFFEResultado(C, BUFFER_C, N, r, indi);
                }
                tarea_recibida++;
            }
        }

		t2 = MPI_Wtime();
		/* Fin de la multiplicacion de matrices */
		printf("Ok!!\n\n");

		if (print) {
		    /* Si el tamanho de la matriz es pequenha entonces se imprime el resultado*/
			printf("Imprimiendo Matriz A...\n");
			imprimirMatriz(A, N, N);
			printf("ok!!\n\n");
			printf("Imprimiendo Matriz B...\n");
			imprimirMatriz(B, N, N);
			printf("ok!!\n\n");
			printf("Imprimiendo Matriz C...\n");
			imprimirMatriz(C, N, N);
			printf("ok!!\n\n");
		}
		printf("\nDuración total de la multilplicacion de matrices %4f segundos\n", t2-t1);

		free(A);
		free(B);
		free(C);
		free(BUFFER_A);
   		free(BUFFER_B);
    	free(BUFFER_C);
	}
	else{

		int tarea;
		int id_destino;
		for (tarea = 0; tarea < subnum*subnum; tarea++ ){

			indice indi = getIndice(tarea, subnum);
            id_destino = tarea % nproc;

            if ( id_destino == myid ){

            	BUFFER_A = (double *)malloc(N*N*sizeof(double));
   				BUFFER_B = (double *)malloc(N*N*sizeof(double));
    			BUFFER_C = (double *)malloc(N*N*sizeof(double));

                /* Se reciben los bloques de entrada, se multiplican y se envia el resultado */
                MPI_Recv(BUFFER_A, r*N, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);
                MPI_Recv(BUFFER_B, N*r, MPI_DOUBLE, 0, 2, MPI_COMM_WORLD, &status);
                multiplicarSubMatriz(BUFFER_A, BUFFER_B, BUFFER_C, N, r);
                MPI_Send(BUFFER_C, r*r, MPI_DOUBLE, 0, 3, MPI_COMM_WORLD);
                
                free(BUFFER_A);
   				free(BUFFER_B);
    			free(BUFFER_C);
            }
		}
	}

    MPI_Finalize();

	return 0;
}

/**
 * Esta funcion mapea una tarea a un indice de la sub-matriz
 */
indice getIndice(int tarea, int subnum)
{
	indice indi;
	indi.i = tarea / subnum;
	indi.j = tarea % subnum;
	return indi;
}

/**
 * Imprime el contenido de una matriz F(filas) x C(columnas)
 */
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

/**
 * Inicializa las matrices en forma generica
 */
void inicializarMatrices(double *A, double *B, double *C, int N)
{
    int i, j;
	for (i=0; i < N; i++){
		for (j=0; j < N; j++) {
			A[i*N+j] = i*N+j;
			B[i*N+j] = i*N+j;
			C[i*N+j]= 0;
		}
	}
}


/**
 * Rellena los buffers con los valores necesarios para calcular
 * el bloque correspondiente a la matriz resultado C
 * Parametros: A y B son las matrices de entrada
 *   indi: indica la primera posicion del bloque a calcular
 */
void cargarBuffers(double *A, double *B, double *BUFFER_A, double *BUFFER_B, int N, int r, indice indi){

    indice indiceA, indiceB;
	indiceA.i = indi.i * r;
	indiceA.j = 0;
	indiceB.i = 0;
	indiceB.j = indi.j * r;

    if ( N % r != 0){
        /* Sobran bloques de tamanhos distintos a r */
        if ( indiceA.i + r > N )
            indiceA.i = N - r;
        if ( indiceB.j + r > N )
            indiceB.j = N - r;
    }

	int i, j;
	/* Rellenamos el buufer A */
	for (i = 0; i < r ; i++){
		for (j = 0; j < N; j++){
			BUFFER_A[i*N+j] = A[ (i+indiceA.i) * N + (j+indiceA.j) ];
		}
	}

	/* Rellenamos el buufer B */
	for (i = 0; i < N ; i++){
		for (j = 0; j < r; j++){
			BUFFER_B[i*r+j] = B[(i+indiceB.i) * N + (j+indiceB.j) ];
		}
	}

}

/**
 * Inicializa el buffer C
 */
void inicializarBufferC(double *BUFFER_C, int r){

	int i, j;
	for (i = 0; i < r ; i++)
		for (j = 0; j < r; j++)
			BUFFER_C[i * r + j] = 0;
}

/**
 * Asigna los datos obtenidos en el BUFFER_C
 * a la matriz resultado C
 */
void asignarBUFFEResultado(double *C, double *BUFFER_C, int N, int r, indice indi){

     indice indiceC;
     indiceC.i = indi.i * r;
	 indiceC.j = indi.j * r;

     if ( N % r != 0){
         /* Sobran bloques de tamanhos distintos a r */
         if ( indiceC.i + r > N )
             indiceC.i = N - r;
         if ( indiceC.j + r > N )
             indiceC.j = N - r;
     }

     int i, j;
     for (i = 0; i < r; i++)
         for (j = 0; j < r; j++)
             C[(i+indiceC.i) * N + (j+indiceC.j)] = BUFFER_C[i * r + j];

}

/**
 * Realiza la multiplicacion correspondiente a un
 * determinado bloque de la matriz resultado C
 * Utiliza los BUFFER_A y BUFFER_B como parametros de entrada y
 * el resultado lo almacena en BUFFER_C
 */
void multiplicarSubMatriz(double *BUFFER_A, double *BUFFER_B, double *BUFFER_C, int N, int r){

    //inicializarBufferC(BUFFER_C, r);

    int i, j, k;
    for (i = 0; i < r; i++){
        for (j = 0; j < r; j++){
        	double temp = 0;
            for (k = 0; k < N; k++){
                temp += BUFFER_A[i * N + k] * BUFFER_B[k * r + j];
            }
            //BUFFER_C[ i * r + j] += temp;
            BUFFER_C[ i * r + j] = temp;
        }
    }
}


/* 
 * Version of get_integer_arg for MPI programs, allowing for printing
 * only in one process and execution of a mpi_finalize function before exit.
 */
long get_integer_arg(int argc, char* argv[], int arg_index, 
        long min_val, const char* description, const char* usage_msg,
        bool print_flag, void (*mpi_finalize) (void) )
{
    if (arg_index >= argc) {
        if (print_flag)
            fprintf(stderr, usage_msg, argv[0]);
        if (mpi_finalize != NULL) mpi_finalize();
        exit(EXIT_FAILURE);
    }
    char* end;
    long rval = strtol(argv[arg_index], &end, 10);
    if (*end != '\0') {
        if (print_flag) {
            fprintf(stderr, "%s debe ser un número entero\n", description);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (mpi_finalize != NULL) mpi_finalize();
        exit(EXIT_FAILURE);
    }
    if (rval < min_val) {
        if (print_flag) {
            fprintf(stderr, "%s debe ser al menos %ld\n", description, min_val);
            fprintf(stderr, usage_msg, argv[0]);
        }
        if (mpi_finalize != NULL) mpi_finalize();
        exit(EXIT_FAILURE);
    }
    return rval;
}

void error_exit(char* msg) {
    if (msg != NULL) {
        fprintf(stderr, "%s", msg);
    }
    //MPI_Finalize();
    exit(EXIT_FAILURE);
}