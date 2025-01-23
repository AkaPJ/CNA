#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <mpi.h>

// Inicializa una matriz con valores aleatorios
void init_matrix(double *matrix, int rows, int cols, int seed_offset) {
    srand((unsigned int)time(NULL) + seed_offset);
    for (int i = 0; i < rows * cols; i++) {
        matrix[i] = rand() * 1.0 / RAND_MAX;
    }
}

// Multiplica dos bloques de matrices y acumula el resultado en C
void matmul_block(double *A, double *B, double *C, int block_size) {
    for (int i = 0; i < block_size; i++) {
        for (int j = 0; j < block_size; j++) {
            for (int k = 0; k < block_size; k++) {
                C[i * block_size + j] += A[i * block_size + k] * B[k * block_size + j];
            }
        }
    }
}

// Algoritmo SUMMA
void SUMMA(MPI_Comm comm_cart, int block_size, double *A_block, double *B_block, double *C_block, int n_blocks) {
    int coords[2];
    MPI_Cart_coords(comm_cart, MPI_Comm_rank(comm_cart, MPI_COMM_WORLD), 2, coords);
    int my_row = coords[0], my_col = coords[1];

    // Crear comunicadores por fila y columna
    MPI_Comm row_comm, col_comm;
    int remain_dims[2] = {1, 0}; // Por fila
    MPI_Cart_sub(comm_cart, remain_dims, &row_comm);
    remain_dims[0] = 0; remain_dims[1] = 1; // Por columna
    MPI_Cart_sub(comm_cart, remain_dims, &col_comm);

    // Buffers temporales para las difusiones
    double *A_temp = malloc(block_size * block_size * sizeof(double));
    double *B_temp = malloc(block_size * block_size * sizeof(double));

    for (int k = 0; k < n_blocks; k++) {
        // Difusión en filas
        if (my_col == k) {
            memcpy(A_temp, A_block, block_size * block_size * sizeof(double));
        }
        MPI_Bcast(A_temp, block_size * block_size, MPI_DOUBLE, k, row_comm);

        // Difusión en columnas
        if (my_row == k) {
            memcpy(B_temp, B_block, block_size * block_size * sizeof(double));
        }
        MPI_Bcast(B_temp, block_size * block_size, MPI_DOUBLE, k, col_comm);

        // Multiplicación local de bloques
        matmul_block(A_temp, B_temp, C_block, block_size);
    }

    free(A_temp);
    free(B_temp);
}

int main(int argc, char **argv) {
    MPI_Init(&argc, &argv);

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (argc != 2) {
        if (rank == 0) {
            fprintf(stderr, "Uso: mpirun -np <num_procesos> ./summa <N>\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int N = atoi(argv[1]); // Tamaño global de la matriz
    int n_proc = sqrt(size); // Procesos por fila/columna
    if (n_proc * n_proc != size || N % n_proc != 0) {
        if (rank == 0) {
            fprintf(stderr, "ERROR: El número de procesos debe ser un cuadrado perfecto y N debe ser divisible por sqrt(p).\n");
        }
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    // Tamaño de los bloques
    int block_size = N / n_proc;
    double *A_block = malloc(block_size * block_size * sizeof(double));
    double *B_block = malloc(block_size * block_size * sizeof(double));
    double *C_block = calloc(block_size * block_size, sizeof(double));

    double *A_global = NULL, *B_global = NULL, *C_global = NULL;
    if (rank == 0) {
        A_global = malloc(N * N * sizeof(double));
        B_global = malloc(N * N * sizeof(double));
        C_global = malloc(N * N * sizeof(double));
        init_matrix(A_global, N, N, 0);
        init_matrix(B_global, N, N, 1);
    }

    // Crear tipos derivados para bloques
    MPI_Datatype block_type, block_type_resized;
    MPI_Type_vector(block_size, block_size, N, MPI_DOUBLE, &block_type);
    MPI_Type_create_resized(block_type, 0, block_size * sizeof(double), &block_type_resized);
    MPI_Type_commit(&block_type_resized);

    // Distribuir bloques de matrices
    MPI_Scatter(A_global, 1, block_type_resized, A_block, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B_global, 1, block_type_resized, B_block, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Crear topología cartesiana
    int dims[2] = {n_proc, n_proc}, periods[2] = {0, 0};
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

    // Ejecutar SUMMA
    double start = MPI_Wtime();
    SUMMA(comm_cart, block_size, A_block, B_block, C_block, n_proc);
    double end = MPI_Wtime();

    // Recolectar resultados
    MPI_Gather(C_block, block_size * block_size, MPI_DOUBLE, C_global, block_size * block_size, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (rank == 0) {
        printf("Tiempo de ejecución: %f segundos\n", end - start);
        free(A_global);
        free(B_global);
        free(C_global);
    }

    free(A_block);
    free(B_block);
    free(C_block);
    MPI_Type_free(&block_type_resized);
    MPI_Finalize();
    return 0;
}