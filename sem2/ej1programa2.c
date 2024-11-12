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

  //Versión secuencial: ijk
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
  printf("Tiempo = %.1f segundos de la multiplicación secuencial ijk.\n",elapsed);
  
  
   //Versión paralela1 ijk
  double **C = (double **) malloc (n*sizeof(double*));
  C[0] = (double *) calloc (n*n,sizeof(double));
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
  printf("Tiempo = %.1f segundos de la multiplicación paralela1 ijk. \n",elapsed);

  double error = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error += fabs(C[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela1 ikj = %.2e\n",error);
  
  

  //Versión paralela2: ikj
  double **C1 = (double **) malloc (n*sizeof(double*));
  C1[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C1[i] = (double *) &(C1[0][i*n]);
  
  t1 = omp_get_wtime();
  #pragma omp parallel for private(i,k,j)
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
  printf("Tiempo = %.1f segundos de la multiplicación paralela2 ikj. \n",elapsed);

  double error1 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error1 += fabs(C1[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela2 ikj = %.2e\n",error1);
  
  
    //Versión paralela3 jik
  double **C2 = (double **) malloc (n*sizeof(double*));
  C2[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C2[i] = (double *) &(C2[0][i*n]);
  
  t1 = omp_get_wtime();
  #pragma omp parallel for private(j,i,k)
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
  printf("Tiempo = %.1f segundos de la multiplicación paralela3 jik. \n",elapsed);

  double error2 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error2 += fabs(C2[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela3 ikj = %.2e\n",error2);


  //Versión paralela4 jki
  double **C3 = (double **) malloc (n*sizeof(double*));
  C3[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C3[i] = (double *) &(C3[0][i*n]);
  
  t1 = omp_get_wtime();
  #pragma omp parallel for private(j,k,i)
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
  printf("Tiempo = %.1f segundos de la multiplicación paralela4 jki. \n",elapsed);

  double error3 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error3 += fabs(C3[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela4 jki = %.2e\n",error3);
  
  
    //Versión paralela5 kij
  double **C4 = (double **) malloc (n*sizeof(double*));
  C4[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C4[i] = (double *) &(C4[0][i*n]);
  
  t1 = omp_get_wtime();
  
  for( k = 0; k < n; k++ ) {
     #pragma omp parallel for private(i,j)
     for( i = 0; i < n; i++ ) {
         for( j = 0; j < n; j++ ) {
            C4[i][j] = C4[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo = %.1f segundos de la multiplicación paralela5 kij. \n",elapsed);

  double error4 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error4 += fabs(C4[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela5 kij = %.2e\n",error4);
  
  
  
      //Versión paralela6 kji
  double **C5 = (double **) malloc (n*sizeof(double*));
  C5[0] = (double *) calloc (n*n,sizeof(double));
  for(i=1; i<n; i++) C5[i] = (double *) &(C5[0][i*n]);
  
  t1 = omp_get_wtime();
 
  for( k = 0; k < n; k++ ) {
     #pragma omp parallel for private(j,i)
     for( j = 0; j < n; j++ ) {
         for( i = 0; i < n; i++ ) {
            C5[i][j] = C5[i][j] + A[i][k] * B[k][j]; 
         }
     }
  }
  t2 = omp_get_wtime();
  elapsed = t2 - t1;
  //print( "C", n, C );
  printf("Tiempo = %.1f segundos de la multiplicación paralela6 kji. \n",elapsed);

  double error5 = 0.0;
  for( int i = 0; i<n*n; i++ ) {
    error5 += fabs(C5[0][i]-Cseq[0][i]);
  }
  printf("Error version paralela6 kji = %.2e\n",error5);
  
 
   
  free(A[0]);
  free(A);
  free(B[0]);
  free(B);
  free(Cseq[0]);
  free(Cseq);
  free(C[0]);
  free(C);
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
