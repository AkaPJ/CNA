
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

/* En este ejercicio se supone un almacenamiento de los elementos de una matriz por columnas,
 * es decir, los elementos de una misma columna están almacenados en posiciones consecutivas de memoria 
 * La siguiente definición (macro) es coherente con esto, siendo m el número de filas de la matriz. 
 */
#define A(i,j) A[(j)*m+(i)]
#define Asaved(i,j) Asaved[(j)*m+(i)]

double norm2( int n, double *v );

int main( int argc, char *argv[] ) {

  MPI_Init(&argc,&argv);
  int size, rank; 
  MPI_Comm_size (MPI_COMM_WORLD, &size); 
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if( argc<3 ) {
    if( rank==0 ) printf("Usage: mpirun -np %d %s n block_size\n",size,argv[0]);
    MPI_Finalize();
    return 0;
  }
  int m, n;
  sscanf(argv[1],"%d",&m);
  int blk;
  sscanf(argv[2],"%d",&blk);
  if( blk*size > m ) {
    if( rank==0 ) printf("blk debería ser más pequeño\n");
    MPI_Finalize();
    return 0;
  }

  double *A, *Asaved;
  int blocks_per_process = (m+blk*size-1)/(blk*size);
  /* Número de columnas de la matriz que hay que reservar (m<=n) */
  n = blocks_per_process*blk*size; 

  if( rank==0 ) {
    A = (double *) malloc( m*n*sizeof(double) );
    Asaved = (double *) malloc( m*n*sizeof(double) ); /* Para comprobar el resultado */
    for( int i=0; i<m; i++ ) {
      A(i,i) = Asaved( i, i ) = ( (double) rand() ) / RAND_MAX;
      for( int j=i+1; j<m; j++ ) {
        A(i,j) = A(j,i) = Asaved(i,j) = Asaved(j,i) = ( (double) rand() ) / RAND_MAX;
      }
    }
  }

  /* Almacenamiento de destino local a cada proceso */
   /* Resercad el espacio justo para almacenar los datos de llegada */

  double *B = (double *) malloc( m*blk*sizeof(double) );
  /* Creación del tipo de dato MPI */
  MPI_Datatype blocktype;
  MPI_Type_vector( m, blk, n, MPI_DOUBLE, &blocktype );

  /* Distribución de la matriz A entre los procesos */
  MPI_Type_commit( &blocktype );
  MPI_Scatter( A, 1, blocktype, B, m*blk, MPI_DOUBLE, 0, MPI_COMM_WORLD );

  /* Recogida de la matriz desde los proceos en la matriz Asaved del root */
  MPI_Gather( B, m*blk, MPI_DOUBLE, Asaved, m*blk, MPI_DOUBLE, 0, MPI_COMM_WORLD );


  /* Comprobación del resultado */
  if( rank == 0 ) {
    for( int i=0; i<m; i++ ) {
      for( int j=0; j<n; j++ ) {
        Asaved( i, j ) -= A( i, j );
      }
    }
    printf("Error = %f\n",norm2(m*m,Asaved));
  }
  MPI_Type_free( &blocktype );
  //MPI_Type_free( &blocktype );
  free(B);
  if( rank==0 ) {
    free(Asaved);
    free(A);
  }
  MPI_Finalize();
  return 0;
}

double norm2( int n, double *v ) {
  double x = 0.0;
  for( int i=0; i<n; i++ ) {
    x += v[i]*v[i];
  }
  return sqrt(x);
}


