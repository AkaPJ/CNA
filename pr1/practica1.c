#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
//cada posicio consecutiva de memoria pertanyen a la mateixa columna

void generaMatriz(double *m, int n) {
    for (int i = 0; i < n; i++) {
        m[i] = (double)rand() / RAND_MAX;
    }
}

void sequential_matrix_multiplication(double* A, double* B, double* C, int n) {
    for(int i = 0; i < n i++) {
        for (int j = 0; j < n; j++) {
            C[i*n+j] = 0.0;
            for(int k = 0; k < n; k++) {
                C[i*n+j] += A[i*n+k] * B[k*n+j];
            }
        }
    }
}

double compute_error(double* sequential_C,double* parallel_C,int n) {
    double error = 0.0;
    for(int i = 0; i < n; i++) {
        error += fabs(sequential_C[i] - parallel_C[i]);
    }
    return error;
}

int main(int argc, char *argv[]) {
    int rank,size,n;
    double *A, *B, *C, *Aloc, *Bloc, *Cloc,*sequential_C;
    double start, end,time,error;
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if(argc != 2) {
        if(rank == 0) {
            fprintf(stderr,"Usage: %s <matrix size>\n", argv[0]);
        }
        MPI_Finalize();
        return -1;
    }

    n = atoi(argv[1]);

    A = (double*)malloc(n*n*sizeof(double));
    B = (double*)malloc(n*n*sizeof(double));
    C = (double*)malloc(n*n*sizeof(double));
    sequential_C = (double*)malloc(n*n*sizeof(double));

    Aloc = (double*)malloc(n*(n/size)*sizeof(double));
    Bloc = (double*)malloc(n*(n/size)*sizeof(double));
    Cloc = (double*)malloc(n*(n/size)*sizeof(double));

    if (rank==0) {
        generaMatriz(A,n);
        generaMatriz(B,n);
        sequential_matrix_multiplication(A,B,sequential_C,n);
    }

    MPI_Bcast(A,n*n,MPI_DOUBLE,0,MPI_COMM_WORLD);
    MPI_Scatter(B,n*(n/size),MPI_DOUBLE,Bloc,n*(n/size),MPI_DOUBLE,0,MPI_COMM_WORLD);

    start = MPI_Wtime();

    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n/size; j++) {
            Cloc[i*(n/size)+j] = 0.0;
            for(int k = 0; k < n; k++) {
                Cloc[i*(n/size)+j] += A[i*n+k] * Bloc[k*(n/size)+j];
            }
        }
    }
    end = MPI_Wtime();
    time = end - start;

    MPI_Gather(Cloc,n*(n/size),MPI_DOUBLE,C,n*(n/size),MPI_DOUBLE,0,MPI_COMM_WORLD);
    if(rank == 0) {
        printf("Multiplicacion paralela finalizada\n");
        error = compute_error(sequential_C,C,n);
        printf("Parallel time: %f\n",time);
        printf("Error: %f\n",error);
    }
    free(A);
    free(B);
    free(C);
    free(Aloc);
    free(Bloc);
    free(Cloc);    
    MPI_Finalize();
    
}