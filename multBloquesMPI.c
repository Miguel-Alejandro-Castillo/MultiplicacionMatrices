#include <mpi.h>


int idP;
int cantP;

int main(int argc, char *argv[]){

	//Sirve para inicializar el entorno MPI.
	MPI_Init(&argc, &argv);
	//Indica la cantidad de procesos en el comunicador.
	MPI_Comm_rank(MPI_COMM_WORLD, &idP);
	//Indica el “rank” (identificador) del proceso dentro de ese comunicador.
	MPI_Comm_size(MPI_COMM_WORLD, &cantP);

	if(idP == 0){
		
	} else {

	}
	MPI_Finalize();
	exit(0);
}