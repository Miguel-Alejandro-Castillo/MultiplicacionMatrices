# Multiplicacion de Matrices Cuadrados por Bloques

* Secuencial
  * Para compilar: gcc -o multBloques multBloques.c
  * Para ejecutar: ./multBloques &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;print?&gt;
  
* MPI 
  * Para compilar: mpicc -fopenmp -o mm2d 2d_ciclico.c
  * Para ejecutar: mpirun -np &lt;processes&gt; --hostfile &lt;file_hosts&gt; mm2d &lt;sizeMatrix&gt; &lt;sizeBlock&gt; &lt;threads&gt; &lt;print?&gt;
