
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpi.h>

#define A(i,j) A[(i)+(j)*n]
#define A_colectiva(i,j) A_colectiva[(i)+(j)*n]

int potencia_original( int n, double *A, double *lambda, double *q );
int potencia_colectiva( int n, double *A, double *lambda, double *q );

int main( int argc, char *argv[] ) {

  MPI_Init(&argc,&argv);
  int nprocs, rank; 
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs); 
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  if( argc<2 ) {
    if( rank==0 ) printf("Usage: mpirun -np %d %s n\n",nprocs,argv[0]);
    MPI_Finalize();
    return 0;
  }

  int n;
  sscanf(argv[1],"%d",&n);
  int cols_per_proc = n/nprocs + (n%nprocs>0?1:0);
  int total_cols = cols_per_proc*nprocs;

  double *A;
  if( rank==0 ) {
    A = (double *) malloc( n*total_cols*sizeof(double) );
    for( int i=0; i<n; i++ ) {
      A(i,i) = ( (double) rand() ) / RAND_MAX;
      for( int j=i+1; j<n; j++ ) {
        A(i,j) = A(j,i) = ( (double) rand() ) / RAND_MAX;
      }
    }
#ifdef PRINT_MATRIZ
    for( int i=0; i<n; i++ ) {
      for( int j=0; j<n; j++ ) {
        printf("%16.10f",A(i,j));
      }
      printf("\n");
    }
#endif
  }
  double *q = (double *) malloc( n*sizeof(double) );

  double lambda;
  
  //Meción del tiempo de la versión original
  double t1, t2;
  //Version original
  t1=MPI_Wtime();
  int k_original = potencia_original( n, A, &lambda, q );
  t2=MPI_Wtime();
  double tfinal = end_tori - start_tori;

  if( rank==0 ) {
    printf("Version original: lambda = %f\n",lambda);
#ifdef PRINT_VECTOR
    for(int i=0; i<n; i++ ) {
      printf("q[%d] = %f\n",i,q[i]);
    }
#endif
    printf("Version original: iteraciones = %d\n",k_original);
    printf("Version original: tiempo = %f\n",tori);

    free(q);
    free(A);
  }

  
  MPI_Barrier(MPI_COMM_WORLD);
  

  double t1_par, t2_par;
  double *A_colectiva;
   if( rank==0 ) {
   	A_colectiva = (double *) malloc( n*total_cols*sizeof(double) );
    	for( int i=0; i<n; i++ ) {
      		A_colectiva(i,i) = ( (double) rand() ) / RAND_MAX;
      		for( int j=i+1; j<n; j++ ) {
        		A_colectiva(i,j) = A_colectiva(j,i) = ( (double) rand() ) / RAND_MAX;
      		}
    }
#ifdef PRINT_MATRIZ
    for( int i=0; i<n; i++ ) {
      for( int j=0; j<n; j++ ) {
        printf("%16.10f",A_colectiva(i,j));
      }
      printf("\n");
    }
#endif
  }
  double *q_colectiva = (double *) malloc( n*sizeof(double) );
  double lambda_colectiva;
  t1_par=MPI_Wtime();
  int k_col = potencia_colectiva( n, A_colectiva, &lambda_colectiva, q_colectiva );
  t2_par=MPI_Wtime();	
  double tpar= t1_par-t2_par;
  if( rank==0 ) {
    printf("Version colectiva: lambda = %f\n",lambda_colectiva);
#ifdef PRINT_VECTOR
    for(int i=0; i<n; i++ ) {
      printf("q[%d] = %f\n",i,q_colectiva[i]);
    }
#endif
    printf("Version colectiva: iteraciones = %d\n",k_col);
    printf("Version colectiva: tiempo = %f\n",tcol);
    free(q_colectiva);
    free(A_colectiva);
  }

	
  MPI_Finalize();
  return 0;
}

int matrizvector( int n, int cols_per_proc, double *A, double *q, double *z ) {
  int nprocs, rank;
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs); 
  for( int i=0; i<n; i++ ) {
    z[i] = 0.0;
    for( int j=0; j<cols_per_proc; j++ ) {
      z[i] += A(i,j)*q[j]; 
    }
  }
  if( rank == 0 ) {
    double *r = (double *) malloc( n*sizeof(double) );
    for( int p=1; p<nprocs; p++ ) {
      MPI_Recv( r, n, MPI_DOUBLE, p, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
      for( int i=0; i<n; i++ ) z[i] += r[i];
    }
    free(r);
    for( int p=1; p<nprocs; p++ ) {
      MPI_Send( z, n, MPI_DOUBLE, p, 123, MPI_COMM_WORLD );
    }
  } else {
    MPI_Send( z, n, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD );
    MPI_Recv( z, n, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  }
  
  return 0;
}

double norm2( int n, double *v ) {
  double x = 0.0;
  for( int i=0; i<n; i++ ) {
    x += v[i]*v[i];
  }
  return sqrt(x);
}

int scal( int n, double alfa, double *v, double *w ) {
  for( int i=0; i<n; i++ ) {
    w[i] = alfa*v[i];
  }
  return 0;
}

double dot( int n, double *v, double *w ) {
  double d = 0.0;
  for( int i=0; i<n; i++ ) {
    d += v[i]*w[i];
  }
  return d;
}

int potencia_original( int n, double *A, double *lambda, double *q ) {

  int nprocs, rank;
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs); 
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  int cols_per_proc = n/nprocs + (n%nprocs>0?1:0);
  int total_cols = cols_per_proc*nprocs;

  double *Aloc = (double *) malloc( n*cols_per_proc*sizeof(double) );
  if( rank==0 ) {
    for( int p=1; p<nprocs; p++ ) {
      MPI_Send( &A[p*n*cols_per_proc], n*cols_per_proc, MPI_DOUBLE, p, 123, MPI_COMM_WORLD );
    }
    MPI_Sendrecv( A, n*cols_per_proc, MPI_DOUBLE, 0, 123, Aloc, n*cols_per_proc, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  } else {
    MPI_Recv( Aloc, n*cols_per_proc, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  }

  double *z = (double *) malloc( n*sizeof(double) );
  double alfa = 0.0;
  int k = 1;
  
  if( rank == 0 ) {
    for( int p=1; p<nprocs; p++ ) {
      MPI_Send( &A(0,0), n, MPI_DOUBLE, p, 123, MPI_COMM_WORLD );
    }
    MPI_Sendrecv( &A(0,0), n, MPI_DOUBLE, 0, 123, z, n, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  } else {
    MPI_Recv( z, n, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  }
  if( rank == 0 ) {
    for( int p=1; p<nprocs; p++ ) {
      MPI_Send( &A(0,0), 1, MPI_DOUBLE, p, 123, MPI_COMM_WORLD );
    }
    MPI_Sendrecv( &A(0,0), 1, MPI_DOUBLE, 0, 123, lambda, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  } else {
    MPI_Recv( lambda, 1, MPI_DOUBLE, 0, 123, MPI_COMM_WORLD, MPI_STATUS_IGNORE );
  }

  while( fabs(*lambda-alfa)>2.220446049250313e-16 ) {
    k++;
    alfa = *lambda;
    scal( n, 1.0/norm2(n,z), z, q );
    matrizvector( n, cols_per_proc, Aloc, &q[rank*cols_per_proc], z );
    *lambda = dot( n, q, z );
  }

  free(Aloc);
  free(z);
  return k;
}

int potencia_colectiva( int n, double *A, double *lambda, double *q ) {

  int nprocs, rank;
  MPI_Comm_size (MPI_COMM_WORLD, &nprocs); 
  MPI_Comm_rank (MPI_COMM_WORLD, &rank);
  int cols_per_proc = n/nprocs + (n%nprocs>0?1:0);
  int total_cols = cols_per_proc*nprocs;
  double local_dot;

  double *Aloc = (double *) malloc( n*cols_per_proc*sizeof(double) );
  MPI_Scatter(A,n*cols_per_proc,MPI_DOUBLE,Aloc,n*cols_per_proc,MPI_DOUBLE,0,MPI_COMM_WORLD);
 


  double *z = (double *) malloc( n*sizeof(double) );
  double alfa = 0.0;
  int k = 1;
  
  if (rank == 0) {   
      *lambda = 1.0; 
  }
  MPI_Bcast(z, n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
  MPI_Bcast(lambda, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);


  while( fabs(*lambda-alfa)>2.220446049250313e-16 ) {
    k++;
    alfa = *lambda;
    scal( n, 1.0/norm2(n,z), z, q );
    matrizvector( n, cols_per_proc, Aloc, &q[rank*cols_per_proc], z );
    *lambda=dot(n,q,z);

  }

  free(Aloc);
  free(z);
  return k;
}

