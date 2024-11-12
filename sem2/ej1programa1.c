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

  //Versión secuencial 1: ijk
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
  printf("Tiempo multiplicación secuencial 1 (ijk) = %.1f\n",elapsed);
  

  //Versión secuencial 2: ikj
  double **C1 = (double **) malloc (n*sizeof(double*));
  C1[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C1[i] = (double *) &(C1[0][i*n]);
  
  t1 = omp_get_wtime();
  for( i = 0; i < n; i++ ) {
     for( k = 0; k < n; k++ ) {
         for( j = 0; j < n; j++ ) {
            C1[i][j] = C1[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo multiplicación secuencial 2 (ikj) = %.1f\n",elapsed);

  double error1 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error1 += fabs(C1[0][i]-Cseq[0][i]);
  }
  printf("Error secuencial 2 (ikj) = %.2e\n",error1);
  
  
    //Versión secuencial 3: jik
  double **C2 = (double **) malloc (n*sizeof(double*));
  C2[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C2[i] = (double *) &(C2[0][i*n]);
  
  t1 = omp_get_wtime();
  for( j = 0; j < n; j++ ) {
     for( i = 0; i < n; i++ ) {
         for( k = 0; k < n; k++ ) {
            C2[i][j] = C2[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo multiplicación secuencial 3 (jik) = %.1f\n",elapsed);

  double error2 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error2 += fabs(C2[0][i]-Cseq[0][i]);
  }
  printf("Error secuencial 3, jik = %.2e\n",error2);


  //Versión secuencial 4: jki
  double **C3 = (double **) malloc (n*sizeof(double*));
  C3[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C3[i] = (double *) &(C3[0][i*n]);
  
  t1 = omp_get_wtime();
  for( j = 0; j < n; j++ ) {
     for( k = 0; k < n; k++ ) {
         for( i = 0; i < n; i++ ) {
            C3[i][j] = C3[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo multiplicación secuencial 4 (jki) = %.1f\n",elapsed);

  double error3 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error3 += fabs(C3[0][i]-Cseq[0][i]);
  }
  printf("Error secuencial 4 (jki) = %.2e\n",error3);
  
  
    //Versión secuencial 5: kij
  double **C4 = (double **) malloc (n*sizeof(double*));
  C4[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C4[i] = (double *) &(C4[0][i*n]);
  
  t1 = omp_get_wtime();
  for( k = 0; k < n; k++ ) {
     for( i = 0; i < n; i++ ) {
         for( j = 0; j < n; j++ ) {
            C4[i][j] = C4[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo multiplicación secuencial 5 (kij) = %.1f\n",elapsed);

  double error4 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error4 += fabs(C4[0][i]-Cseq[0][i]);
  }
  printf("Error secuencial 5 (kij) = %.2e\n",error4);
  
  
  
      //Versión secuencial 6: kji
  double **C5 = (double **) malloc (n*sizeof(double*));
  C5[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C5[i] = (double *) &(C5[0][i*n]);
  
  t1 = omp_get_wtime();
  for( k = 0; k < n; k++ ) {
     for( j = 0; j < n; j++ ) {
         for( i = 0; i < n; i++ ) {
            C5[i][j] = C5[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo multiplicación secuencial 6 (kji) = %.1f\n",elapsed);

  double error5 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error5 += fabs(C5[0][i]-Cseq[0][i]);
  }
  printf("Error secuencial 6 (kji) = %.2e\n",error5);
  
 
   
  free(A[0]);
  free(A);
  free(B[0]);
  free(B);
  free(Cseq[0]);
  free(Cseq);
  free(C1[0]);
  free(C1);
  free(C2[0]);
  free(C2);
  free(C3[0]);
  free(C3);
  free(C4[0]);
  free(C4);
  free(C5[0]);
  free(C5);
   return 0;
}
