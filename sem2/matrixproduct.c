#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>

void print( const char *matriz, const int n, double **M ) {
  printf("%s = [\n",matriz);
  for( int i = 0; i<n; i++ ) {
    for( int j = 0; j<n; j++ ) {
      printf("%9.4f", M[i][j]);
    }
    printf("\n");
  }
  printf("]\n");
}

int main(int argc, char *argv[]){

  if( argc<2 ) {
    printf("Usage: <matrix_size>\n");
    return 1;
  }
    
  int n, i = 0, j = 0, k = 0;
  sscanf(argv[1],"%d",&n);

  double **A = (double **) malloc (n*sizeof(double*));
  double **B = (double **) malloc (n*sizeof(double*));
  double **Cseq = (double **) malloc (n*sizeof(double*));

  //Inicializar matrices 
  A[0] = (double *) malloc (n*n*sizeof(double));
  B[0] = (double *) malloc (n*n*sizeof(double));
  Cseq[0] = (double *) malloc (n*n*sizeof(double));
  for( i = 0; i<n*n; i++ ) {
    A[0][i] = (double) rand()/RAND_MAX * 2.0 - 1.0;
    B[0][i] = (double) rand()/RAND_MAX * 2.0 - 1.0;
    Cseq[0][i] = 0;
  }
  for(i=1; i<n; i++) {
    A[i] = (double *) &(A[0][i*n]);
    B[i] = (double *) &(B[0][i*n]);
    Cseq[i] = (double *) &(Cseq[0][i*n]);
  }
  //print( "A", n, A );
  //print( "B", n, B );

  double t1 = omp_get_wtime();
  for( i = 0; i < n; i++ ) {
     for( j = 0; j < n; j++ ) {
         for( k = 0; k < n; k++ ) {
            Cseq[i][j] = Cseq[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  double t2 = omp_get_wtime();
  double elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo = %.1f segundos de la multiplicación secuencial\n",elapsed);

  double **C = (double **) malloc (n*sizeof(double*));
  C[0] = (double *) calloc (n*n,sizeof(double));
  /*************************************************************************/
  /* # Código a paralelizar #                                            
   * Esto se ha de paralelizar añadiendo directivas OpenMP y 
   * sin modificar el código                                               */
  /*************************************************************************/
  for(i=1; i<n; i++) C[i] = (double *) &(C[0][i*n]);
  t1 = omp_get_wtime();
  #pragma omp parallel for private(i,j,k)
  for( i = 0; i < n; i++ ) {
     for( j = 0; j < n; j++ ) {
         for( k = 0; k < n; k++ ) {
            C[i][j] = C[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo = %.1f segundos de la multiplicación paralela\n",elapsed);

  /* Comprobación del resultado */
  /* Si la paralelización es correcta el error debe ser pequeño o cero */
  double error = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error += fabs(C[0][i]-Cseq[0][i]);
  }
  printf("Error = %.2e\n",error);
  /*************************************************************************/
  /* Código a paralelizar                                                  */
  /*************************************************************************/
   
  free(A[0]);
  free(A);
  free(B[0]);
  free(B);
  free(Cseq[0]);
  free(Cseq);
  free(C[0]);
  free(C);

   return 0;
}
