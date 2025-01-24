#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include <mpi.h>

void print(const char *matriz, const int n, double **M) {
    printf("%s = [\n", matriz);
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++) {
            printf("%9.4f", M[i][j]);
        }
        printf("\n");
    }
    printf("]\n");
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv); // Inicializar MPI

    int rank, size;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank); // ID del proceso actual
    MPI_Comm_size(MPI_COMM_WORLD, &size); // Número total de procesos

    if (argc < 2) {
        if (rank == 0) printf("Usage: <matrix_size>\n");
        MPI_Finalize();
        return 1;
    }

    int n;
    sscanf(argv[1], "%d", &n);

    if (n % size != 0) {
        if (rank == 0) {
            printf("El tamaño de la matriz (%d) debe ser divisible por el número de procesos MPI (%d).\n", n, size);
        }
        MPI_Finalize();
        return 1;
    }

    // Número de filas por proceso
    int rows_per_process = n / size;
    int start = rank * rows_per_process;
    int end = start + rows_per_process;

    // Reservar memoria para las matrices
    double *A = (double *)malloc(n * sizeof(double *));
    double *B = (double *)malloc(n * sizeof(double *));
    double *Cseq = (double *)malloc(n * sizeof(double *));
    double **C = NULL;

    A[0] = (double *)malloc(n * n * sizeof(double));
    B[0] = (double *)malloc(n * n * sizeof(double));
    Cseq[0] = (double *)malloc(n * n * sizeof(double));
    for (int i = 1; i < n; i++) {
        A[i] = &(A[0][i * n]);
        B[i] = &(B[0][i * n]);
        Cseq[i] = &(Cseq[0][i * n]);
    }

    // Inicializar matrices A, B y Cseq
    for (int i = 0; i < n * n; i++) {
        A[0][i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
        B[0][i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
        Cseq[0][i] = 0;
    }

    // Versión secuencial: ikj
    if (rank == 0) {
        double t1 = omp_get_wtime();
        for (int i = 0; i < n; i++) {
            for (int k = 0; k < n; k++) {
                for (int j = 0; j < n; j++) {
                    Cseq[i][j] += A[i][k] * B[k][j];
                }
            }
        }
        double t2 = omp_get_wtime();
        double elapsed = t2 - t1;
        printf("Tiempo de la multiplicación secuencial = %.1f segundos.\n", elapsed);
    }

    // Broadcast de la matriz B a todos los procesos
    MPI_Bcast(B[0], n * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Enviar/Recibir las filas correspondientes de A a cada proceso
    double *local_A = (double *)malloc(rows_per_process * sizeof(double *));
    local_A[0] = (double *)malloc(rows_per_process * n * sizeof(double));
    for (int i = 1; i < rows_per_process; i++) {
        local_A[i] = &(local_A[0][i * n]);
    }
    MPI_Scatter(A[0], rows_per_process * n, MPI_DOUBLE, local_A[0], rows_per_process * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Reservar espacio para la parte local de C
    double *local_C = (double *)malloc(rows_per_process * sizeof(double *));
    local_C[0] = (double *)calloc(rows_per_process * n, sizeof(double));
    for (int i = 1; i < rows_per_process; i++) {
        local_C[i] = &(local_C[0][i * n]);
    }

    int i, k, j;

    // Multiplicación paralela
    double t1 = MPI_Wtime();
    #pragma omp parallel for private(i, k, j)
    for (i = 0; i < rows_per_process; i++) {
        for (k = 0; k < n; k++) {
            for (j = 0; j < n; j++) {
                local_C[i][j] += local_A[i][k] * B[k][j];
            }
        }
    }
    double t2 = MPI_Wtime();
    double elapsed = t2 - t1;

    // Reunir los resultados de C en el proceso 0
    if (rank == 0) {
        C = (double **)malloc(n * sizeof(double *));
        C[0] = (double *)malloc(n * n * sizeof(double));
        for (int i = 1; i < n; i++) {
            C[i] = &(C[0][i * n]);
        }
    }
    MPI_Gather(local_C[0], rows_per_process * n, MPI_DOUBLE, (rank == 0 ? C[0] : NULL), rows_per_process * n, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    // Imprimir resultados
    if (rank == 0) {
        printf("Tiempo de la multiplicación paralela (MPI + OpenMP): %.2f segundos.\n", elapsed);

        double error1 = 0.0;
        for (int i = 0; i < n * n; i++) {
            error1 += fabs(C[0][i] - Cseq[0][i]);
        }
        printf("Error version paralela = %.2e\n", error1);
        free(C[0]);
        free(C);
    }

    // Liberar memoria
    free(A[0]);
    free(A);
    free(B[0]);
    free(B);
    free(Cseq[0]);
    free(Cseq);
    free(local_A[0]);
    free(local_A);
    free(local_C[0]);
    free(local_C);

    MPI_Finalize();
    return 0;
}