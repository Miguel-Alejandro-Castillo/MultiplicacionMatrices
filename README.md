# Multiplicación de Matrices Cuadradas por Bloques MPI + OpenMP

* Secuencial
  * Para compilar: gcc -o matmul-blk matmul-blk.c
  * Para ejecutar: ./matmul-blk &lt;sizeMatrix&gt; &lt;sizeBlock&gt;
  
* Paralelo(Master Worker)
  *MPI
     * Para compilar: mpicc -fopenmp -o multMatricesMW multMatricesMW.c
     * Para ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; mm2d &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt;
                     <br/>  &lt;threads&gt; &nbsp;&nbsp;  =  &nbsp;&nbsp; 1
  
  * MPI + OpenMP
      * Para compilar: idem caso MPI
      * Para ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; mm2d &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt;
                        <br/> &lt;threads&gt; &nbsp;&nbsp; > &nbsp;&nbsp; 1
                    
 
&lt;sizeMatriz&gt; => tamaño de la matriz <br/>
&lt;sizeBlock&gt; => tamaño del bloque <br/>
&lt;processes&gt; => cantidad de nodos <br/>
&lt;file_hosts&gt; => archivo hostfile <br/>
&lt;threads&gt; =>  cantidad de threads <br/>
