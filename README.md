# Multiplicación de Matrices Cuadradas por Bloques MPI + OpenMP

* Secuencial
  * Compilar: gcc matmul-block.c -O2 -o matmul-block
  * Ejecutar: ./matmul-block &lt;sizeMatrix&gt; &lt;sizeBlock&gt;
  
* Paralelo(Master Worker)
  * MPI
     * Compilar: mpicc matmul-block-parallel.c -fopenmp -O2 -o matmul-block-parallel
     * Ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; matmul-block-parallel &lt;sizeMatrix&gt; &lt;sizeBlock&gt; 1
  
  * MPI + OpenMP
      * Compilar: idem caso MPI
      * Ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; matmul-block-parallel &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt;
                        <br/> &lt;threads&gt; &nbsp;&nbsp; > &nbsp;&nbsp; 1
                    
 
    &lt;sizeMatriz&gt; => tamaño de la matriz <br/>
    &lt;sizeBlock&gt; => tamaño del bloque <br/>
    &lt;processes&gt; => cantidad de nodos <br/>
    &lt;file_hosts&gt; => archivo hostfile <br/>
    &lt;threads&gt; =>  cantidad de threads <br/>
