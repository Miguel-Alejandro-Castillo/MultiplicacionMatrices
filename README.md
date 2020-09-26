# Multiplicación de Matrices Cuadradas por Bloques MPI + OpenMP

* Secuencial
  * Para compilar: gcc -o multBloques multBloques.c
  * Para ejecutar: ./multBloques &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;print?&gt;
  
* Paralelo(Master Worker)
  * Para compilar: mpicc -fopenmp -o multMatricesMW multMatricesMW.c
  * Para ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; mm2d &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt;
                  <br/>  &lt;threads&gt; &nbsp;&nbsp;  =  &nbsp;&nbsp; 1
  
* MPI + OpenMP
  * Para compilar: idem caso MPI
  * Para ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; mm2d &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt; &lt;print?&gt;
                    <br/> &lt;threads&gt; &nbsp;&nbsp; > &nbsp;&nbsp; 1
                    
 
  &lt;sizeMatriz&gt; => tamaño de la matriz <br/>
  &lt;sizeBlock&gt; => tamaño del bloque <br/>
  &lt;print?&gt; => 1 imprime el resultado, 0 caso contrario <br/>
  &lt;processes&gt; => cantidad de nodos <br/>
  &lt;file_hosts&gt; => archivo hostfile <br/>
  &lt;threads&gt; =>  cantidad de threads <br/>
