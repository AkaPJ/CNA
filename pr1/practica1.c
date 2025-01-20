#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Función para inicializar una matriz con valores aleatorios
void initialize_matrix(double* matrix, int size) {
    for (int i = 0; i < size * size; i++) {
        matrix[i] = rand() % 10;
    }
}

// Función para multiplicar matrices secuencialmente
void sequential_matrix_multiplication(double* A, double* B, double* C, int size) {
    for (int i = 0; i < size; i++) {
        for (int j = 0; j < size; j++) {
            C[i * size + j] = 0;
            for (int k = 0; k < size; k++) {
                C[i * size + j] += A[i * size + k] * B[k * size + j];
            }
        }
    }
}

// Función para multiplicar matrices en paralelo utilizando MPI
void parallel_matrix_multiplication(double* A, double* B, double* C, int size, int rank, int num_procs) {
    int rows_per_proc = size / num_procs; // Dividimos las filas entre los procesos
    double* local_A = (double*)malloc(rows_per_proc * size * sizeof(double));
    double* local_C = (double*)malloc(rows_per_proc * size * sizeof(double));

    // Distribuir filas de A entre los procesos
    MPI_Scatter(A, rows_per_proc * size, MPI_DOUBLE, local_A, rows_per_proc * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // Difundir toda la matriz B a todos los procesos
    MPI_Bcast(B, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Inicializar la matriz local de resultados
    for (int i = 0; i < rows_per_proc; i++) {
        for (int j = 0; j < size; j++) {
            local_C[i * size + j] = 0.0;
            for (int k = 0; k < size; k++) {
                local_C[i * size + j] += local_A[i * size + k] * B[k * size + j];
            }
        }
    }

    // Recolectar los resultados en la matriz C
    MPI_Gather(local_C, rows_per_proc * size, MPI_DOUBLE, C, rows_per_proc * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    free(local_A);
    free(local_C);
}

// Función para calcular el error entre dos matrices
double calculate_error(double* C_seq, double* C_par, int size) {
    double error = 0.0;
    for (int i = 0; i < size * size; i++) {
        error += fabs(C_seq[i] - C_par[i]);
    }
    return error;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int rank, num_procs;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &num_procs);

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "Usage: %s <matrix size>\n", argv[0]);
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int size = atoi(argv[1]);
    if (size % num_procs != 0) {
        if (rank == 0) {
            fprintf(stderr, "Matrix size must be divisible by the number of processes\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    double* A = (double*)malloc(size * size * sizeof(double));
    double* B = (double*)malloc(size * size * sizeof(double));
    double* C_seq = (double*)malloc(size * size * sizeof(double));
    double* C_par = (double*)malloc(size * size * sizeof(double));

    if (rank == 0) {
        initialize_matrix(A, size);
        initialize_matrix(B, size);
    }

    // Asegurarse de que todos los procesos tengan las mismas matrices A y B
    MPI_Bcast(A, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Bcast(B, size * size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    parallel_matrix_multiplication(A, B, C_par, size, rank, num_procs);

    if (rank == 0) {
        sequential_matrix_multiplication(A, B, C_seq, size);
        double error = calculate_error(C_seq, C_par, size);
        printf("Error between sequential and parallel results: %f\n", error);
    }

    free(A);
    free(B);
    free(C_seq);
    free(C_par);

    MPI_Finalize();
    return 0;
}