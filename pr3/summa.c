#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>

#define N 4  // Número de bloques
#define BLOCK_SIZE 4  // Tamaño del bloque

void matrix_multiply(double *A, double *B, double *C, int n);
void print_matrix(double *matrix, int n);

int main(int argc, char *argv[]) {
    int rank, size;
    double A[N * BLOCK_SIZE][N * BLOCK_SIZE], B_matrix[N * BLOCK_SIZE][N * BLOCK_SIZE], C[N * BLOCK_SIZE][N * BLOCK_SIZE] = {0};
    double local_A[BLOCK_SIZE][BLOCK_SIZE], local_B[BLOCK_SIZE][BLOCK_SIZE], local_C[BLOCK_SIZE][BLOCK_SIZE] = {0};
    
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Inicializar matrices A y B en el proceso 0
    if (rank == 0) {
        // Inicializar A y B con valores aleatorios
        for (int i = 0; i < N * BLOCK_SIZE; i++) {
            for (int j = 0; j < N * BLOCK_SIZE; j++) {
                A[i][j] = rand() % 10;
                B_matrix[i][j] = rand() % 10;
            }
        }
    }

    // Distribuir bloques de columnas de A y B a cada proceso
    for (int j = 0; j < N; j++) {
        MPI_Scatter(A[j * BLOCK_SIZE], BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, local_A, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);
        MPI_Scatter(B_matrix[j * BLOCK_SIZE], BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, local_B, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, 0, MPI_COMM_WORLD);

        // Calcular el producto parcial
        for (int k = 0; k < N; k++) {
            matrix_multiply(local_A, local_B, local_C, BLOCK_SIZE);
        }

        // Sumar los resultados locales en la matriz C global
        MPI_Reduce(local_C, C, BLOCK_SIZE * BLOCK_SIZE, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
    }

    // Imprimir el resultado en el proceso 0
    if (rank == 0) {
        printf("Resultado C:\n");
        print_matrix((double *)C, N * BLOCK_SIZE);
    }

    MPI_Finalize();
    return 0;
}

void matrix_multiply(double *A, double *B, double *C, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            for (int k = 0; k < n; k++) {
                C[i * n + j] += A[i * n + k] * B[k * n + j];
            }
        }
    }
}

void print_matrix(double *matrix, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%0.2f ", matrix[i * n + j]);
        }
        printf("\n");
    }
}
