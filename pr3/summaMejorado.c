#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

#define TOL 1e-4 // Tolerancia para la validación

int N, myrank;

// Inicializa una matriz con valores aleatorios
void init_matrix(double *matr, const int rows, const int cols, int seed_offset) {
    srand((unsigned int)time(NULL) + seed_offset);
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            matr[j * cols + i] = rand() * 1.0 / RAND_MAX;
        }
    }
}

// Imprime una matriz
void print_matrix(const int rows, const int cols, const double *matr) {
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            printf("%12.6f", matr[j * cols + i]);
        }
        printf("\n");
    }
}

// Multiplicación de matrices básica
void matmul_naive(const int m, const int n, const int k, const double *A, const double *B, double *C) {
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < k; ++i) {
            C[j * k + i] = 0.0;
            for (int l = 0; l < n; ++l) {
                C[j * k + i] += A[j * n + l] * B[l * k + i];
            }
        }
    }
}

// Implementación del algoritmo SUMMA
void SUMMA(MPI_Comm comm_cart, const int mb, double *A_loc, double *B_loc, double *C_loc) {
    int coords[2];
    MPI_Cart_coords(comm_cart, myrank, 2, coords);
    int my_col = coords[0];
    int my_row = coords[1];

    MPI_Comm row_comm, col_comm;
    int remain_dims[2] = {1, 0};
    MPI_Cart_sub(comm_cart, remain_dims, &row_comm);
    remain_dims[0] = 0; remain_dims[1] = 1;
    MPI_Cart_sub(comm_cart, remain_dims, &col_comm);

    double *A_temp = (double *)malloc(mb * mb * sizeof(double));
    double *B_temp = (double *)malloc(mb * mb * sizeof(double));

    for (int bcast_root = 0; bcast_root < N / mb; ++bcast_root) {
        if (my_col == bcast_root) memcpy(A_temp, A_loc, mb * mb * sizeof(double));
        MPI_Bcast(A_temp, mb * mb, MPI_DOUBLE, bcast_root, row_comm);

        if (my_row == bcast_root) memcpy(B_temp, B_loc, mb * mb * sizeof(double));
        MPI_Bcast(B_temp, mb * mb, MPI_DOUBLE, bcast_root, col_comm);

        matmul_naive(mb, mb, mb, A_temp, B_temp, C_loc);
    }

    free(A_temp);
    free(B_temp);
}

// Verifica el resultado contra una multiplicación secuencial
int validate_result(double *C, double *A, double *B, int n) {
    double *C_seq = (double *)calloc(n * n, sizeof(double));
    matmul_naive(n, n, n, A, B, C_seq);
    for (int i = 0; i < n * n; ++i) {
        if (fabs(C[i] - C_seq[i]) > TOL) {
            free(C_seq);
            return 0;
        }
    }
    free(C_seq);
    return 1;
}

int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);

    if (argc != 2) {
        if (myrank == 0)
            fprintf(stderr, "Uso: mpirun -np <num_procesos> ./summa <N>\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    N = atoi(argv[1]);
    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int n_proc_rows = sqrt(nprocs);

    if (n_proc_rows * n_proc_rows != nprocs || N % n_proc_rows != 0) {
        if (myrank == 0)
            fprintf(stderr, "ERROR: N debe ser divisible por el número de procesos por fila.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int dims[2] = {n_proc_rows, n_proc_rows};
    int periods[2] = {0, 0};
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, 2, dims, periods, 0, &comm_cart);

    int mb = N / n_proc_rows;
    double *A_loc = (double *)malloc(mb * mb * sizeof(double));
    double *B_loc = (double *)malloc(mb * mb * sizeof(double));
    double *C_loc = (double *)calloc(mb * mb, sizeof(double));

    double *A_glob = NULL, *B_glob = NULL, *C_glob = NULL;
    if (myrank == 0) {
        A_glob = (double *)malloc(N * N * sizeof(double));
        B_glob = (double *)malloc(N * N * sizeof(double));
        C_glob = (double *)malloc(N * N * sizeof(double));
        init_matrix(A_glob, N, N, 0);
        init_matrix(B_glob, N, N, 1);
    }

    MPI_Scatter(A_glob, mb * mb, MPI_DOUBLE, A_loc, mb * mb, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    MPI_Scatter(B_glob, mb * mb, MPI_DOUBLE, B_loc, mb * mb, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    double tstart = MPI_Wtime();
    SUMMA(comm_cart, mb, A_loc, B_loc, C_loc);
    double tend = MPI_Wtime();

    MPI_Gather(C_loc, mb * mb, MPI_DOUBLE, C_glob, mb * mb, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("Tiempo de ejecución de SUMMA: %f segundos\n", tend - tstart);
        if (validate_result(C_glob, A_glob, B_glob, N))
            printf("Resultado validado correctamente.\n");
        else
            printf("Error en la validación del resultado.\n");
        free(A_glob);
        free(B_glob);
        free(C_glob);
    }

    free(A_loc);
    free(B_loc);
    free(C_loc);
    MPI_Finalize();
    return 0;
}
