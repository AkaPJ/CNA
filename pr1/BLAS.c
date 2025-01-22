#include <stdlib.h>
#include <stdio.h>
#include <mpi.h>
#include <math.h>
#include <cblas.h> // Incluir la biblioteca BLAS

//cada posicio consecutiva de memoria pertanyen a la mateixa columna
void print_matrix(double* matrix, int n) {
    for(int i = 0; i < n; i++) {
        for(int j = 0; j < n; j++) {
            printf("%4.1f ", matrix[j*n + i]); // Acceso en orden de columnas
        }
        printf("\n");
    }
}

void generaMatriz(double *m, int n) {
    for (int i = 0; i < n*n; i++) {
        m[i] = rand() % 100; // Usar valores más pequeños para evitar overflow
    }
}

void sequential_matrix_multiplication(double* A, double* B, double* C, int n) {
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n, 1.0, A, n, B, n, 0.0, C, n);
}

double compute_error(double* sequential_C, double* parallel_C, int n) {
    double error = 0.0;
    for(int i = 0; i < n*n; i++) {
        error += fabs(sequential_C[i] - parallel_C[i]);
    }
    return error;
}

int main(int argc, char *argv[]) {
    int rank, size, n;
    double *A, *B, *C, *Aloc, *Bloc, *Cloc, *sequential_C;
    double start, end, time, error;

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

    if (rank == 0) {
        generaMatriz(A, n);
        generaMatriz(B, n);
        sequential_matrix_multiplication(A, B, sequential_C, n);
    }

    MPI_Bcast(A, n*n, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B, n*(n/size), MPI_DOUBLE, Bloc, n*(n/size), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    start = MPI_Wtime();

    // Multiplicación de matrices usando BLAS
    for(int i = 0; i < n; i++) {
        cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 1, n/size, n, 1.0, &A[i], n, Bloc, n, 0.0, &Cloc[i], n);
    }

    end = MPI_Wtime();
    time = end - start;

    MPI_Gather(Cloc, n*(n/size), MPI_DOUBLE, C, n*(n/size), MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if(rank == 0) {
        printf("Multiplicacion paralela finalizada\n");
        error = compute_error(sequential_C, C, n);
        printf("Parallel time: %f\n", time);
        printf("Error: %f\n", error);
        print_matrix(C, n);
        print_matrix(sequential_C, n);
    }

    free(A);
    free(B);
    free(C);
    free(Aloc);
    free(Bloc);
    free(Cloc);
    free(sequential_C);
    MPI_Finalize();

    return 0;
}