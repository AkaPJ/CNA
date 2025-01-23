#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

void matrix_multiply(double *A, double *B, double *C, int n) {
    for (int i = 0; i < n * n; i++) {
        C[i] = 0.0;
    }

    for (int i = 0; i < n; i++) {
        for (int k = 0; k < n; k++) {
            for (int j = 0; j < n; j++) {
                C[i + j * n] += A[i + k * n] * B[k + j * n];
            }
        }
    }
}

void initialize_matrix(double *matrix, int rows, int cols) {
    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
}

void print_matrix(const char *name, double *matrix, int rows, int cols) {
    printf("%s = [\n", name);
    for (int i = 0; i < rows; i++) {
        for (int j = 0; j < cols; j++) {
            printf("%9.4f ", matrix[i * cols + j]);
        }
        printf("\n");
    }
    printf("]\n");
}

void verify_results(double *C_parallel, double *C_sequential, int n) {
    double max_diff = 0.0;
    for (int i = 0; i < n * n; i++) {
        double diff = fabs(C_parallel[i] - C_sequential[i]);
        if (diff > max_diff) max_diff = diff;
    }
    printf("Error entre multiplicación paralela y secuencial: %e\n", max_diff);
}

void block_multiply(double *A_block, double *B_block, double *C_block, int b) {
    for (int i = 0; i < b; i++) {
        for (int k = 0; k < b; k++) {
            for (int j = 0; j < b; j++) {
                C_block[i + j * b] += A_block[i + k * b] * B_block[k + j * b];
            }
        }
    }
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);  // Get process rank
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Get total number of processes

    if (argc < 3) {
        if (rank == 0) {
            printf("Usage: %s <matrix_size> <block_size>\n", argv[0]);
        }
        MPI_Finalize();
        return 1;
    }

    int n;
    int b;
    sscanf(argv[1], "%d", &n); // Matrix size (n x n)
    sscanf(argv[2], "%d", &b); // Block size

    if (n % b != 0) {
        if (rank == 0) {
            printf("Error: n tiene que ser multiplo de b\n");
        }
        MPI_Finalize();
        return 1;
    }

    int N = n / b;

    if (size != N * N) {
        if (rank == 0)
            printf("Error: El número de procesos tiene que ser un cuadrado perfecto\n");
        MPI_Finalize();
        return -1;
    }

    MPI_Comm cart_comm;
    int dims[2];
    dims[0] = N;
    dims[1] = N;
    int periods[2];
    periods[0] = 0;
    periods[1] = 0;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 1, &cart_comm);

    int coords[2];
    MPI_Cart_coords(cart_comm, rank, 2, coords);

    double *A = NULL;
    double *B = NULL;
    double *C = NULL;
    double *C_sequential = NULL;

    if (rank == 0) {
        A = malloc(n * n * sizeof(double));
        B = malloc(n * n * sizeof(double));
        C = malloc(n * n * sizeof(double));
        C_sequential = malloc(n * n * sizeof(double));

        // Inicializar matrices
        initialize_matrix(A, n, n);
        initialize_matrix(B, n, n);

        // Imprimir matrices iniciales
        print_matrix("A", A, n, n);
        print_matrix("B", B, n, n);

        // Multiplicación secuencial
        double seq_start_time = MPI_Wtime(); // Tiempo inicial
        for (int i = 0; i < n * n; i++) {
            C_sequential[i] = 0.0;
        }
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    C_sequential[i + j * n] += A[i + k * n] * B[k + j * n];
                }
            }
        }
        double seq_end_time = MPI_Wtime(); // Tiempo final
        printf("Tiempo secuencial: %f\n", seq_end_time - seq_start_time);

        print_matrix("C_sequential", C_sequential, n, n);
    }

    // Create row communicator
    MPI_Comm row_comm;
    int remain_dims[2] = {0, 1};
    MPI_Cart_sub(cart_comm, remain_dims, &row_comm);

    // Create column communicator
    MPI_Comm column_comm;
    remain_dims[0] = 1;
    remain_dims[1] = 0;
    MPI_Cart_sub(cart_comm, remain_dims, &column_comm);

    MPI_Barrier(MPI_COMM_WORLD); // Synchronize processes
    double start_time = MPI_Wtime();

    // Scatter columns between first row processes
    double *A_block_column = NULL;
    double *B_block_column = NULL;
    if (coords[0] == 0) {
        A_block_column = malloc(n * b * sizeof(double));
        B_block_column = malloc(n * b * sizeof(double));
        MPI_Scatter(A, n * b, MPI_DOUBLE, A_block_column, n * b, MPI_DOUBLE, 0, row_comm);
        MPI_Scatter(B, n * b, MPI_DOUBLE, B_block_column, n * b, MPI_DOUBLE, 0, row_comm);
    }

    // Create block type
    MPI_Datatype block;
    MPI_Type_vector(b, b, n, MPI_DOUBLE, &block);
    MPI_Type_commit(&block);
    MPI_Datatype block_resized;
    MPI_Type_create_resized(block, 0, b * sizeof(double), &block_resized);
    MPI_Type_commit(&block_resized);

    // Scatter blocks
    double *A_block_loc = malloc(b * b * sizeof(double));
    double *B_block_loc = malloc(b * b * sizeof(double));
    MPI_Scatter(A_block_column, 1, block_resized, A_block_loc, b * b, MPI_DOUBLE, 0, column_comm);
    MPI_Scatter(B_block_column, 1, block_resized, B_block_loc, b * b, MPI_DOUBLE, 0, column_comm);

    if (coords[0] == 0) {
        free(A_block_column);
        free(B_block_column);
    }

    double *A_block = malloc(b * b * sizeof(double));
    double *B_block = malloc(b * b * sizeof(double));
    double *C_block = malloc(b * b * sizeof(double));

    for (int i = 0; i < b * b; ++i) {
        C_block[i] = 0.0;
    }

    for (int i = 0; i < N; ++i) {
        if (coords[1] == i) {
            memcpy(A_block, A_block_loc, b * b * sizeof(double));
        }
        MPI_Bcast(A_block, b * b, MPI_DOUBLE, i, row_comm);

        if (coords[0] == i) {
            memcpy(B_block, B_block_loc, b * b * sizeof(double));
        }
        MPI_Bcast(B_block, b * b, MPI_DOUBLE, i, column_comm);

        block_multiply(A_block, B_block, C_block, b);
    }

    MPI_Gather(C_block, b * b, MPI_DOUBLE, C, b * b, MPI_DOUBLE, 0, cart_comm);

    double end_time = MPI_Wtime();

    // In the root process, rearrange blocks to assemble the global matrix
    if (rank == 0) {
        printf("Time elapsed: %f seconds\n", end_time - start_time);

        double *temp = (double *)malloc(n * n * sizeof(double));
        memset(temp, 0, n * n * sizeof(double));

        // Reorder blocks correctly
        for (int pi = 0; pi < N; pi++) {
            for (int pj = 0; pj < N; pj++) {
                // Get the rank of the process in the Cartesian grid
                int coords[2] = {pj, pi};
                int block_id;
                MPI_Cart_rank(cart_comm, coords, &block_id);

                // Process each element of the block
                for (int i = 0; i < b; i++) {
                    for (int j = 0; j < b; j++) {
                        int global_i = pi * b + i;
                        int global_j = pj * b + j;
                        int block_idx = i * b + j;
                        temp[global_i * n + global_j] = C[block_id * b * b + block_idx];
                    }
                }
            }
        }

        memcpy(C, temp, n * n * sizeof(double));
        free(temp);

        print_matrix("C_parallel", C, n, n);
        verify_results(C, C_sequential, n);

        free(C_sequential);
        free(A);
        free(B);
        free(C);
    }

    free(A_block);
    free(B_block);
    free(C_block);
    free(A_block_loc);
    free(B_block_loc);

    MPI_Finalize();
    return 0;
}