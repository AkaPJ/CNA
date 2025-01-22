#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include <mpi.h>

// Constante de tolerancia para la validación de la multiplicación
#define TOL 1e-4

// Tamaño global de la matriz cuadrada (definimos que m = n = k)
int N;
int myrank;


//crea matriz random
void init_matrix(double *matr, const int rows, const int cols) {
    srand((unsigned int)time(NULL) + myrank); // Semilla diferente por proceso
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            matr[j*cols + i] = rand() * 1.0 / RAND_MAX;
        }
    }
}

//muestra la matriz por pantalla
void print_matrix(const int rows, const int cols, const double *matr) {
    for (int j = 0; j < rows; ++j) {
        for (int i = 0; i < cols; ++i) {
            printf("%12.6f", matr[j*cols + i]);
        }
        printf("\n");
    }
    printf("\n\n\n");
}

//multiplicacion de matrices
void matmul_naive(const int m, const int n, const int k, const double *A, const double *B, double *C) {
    for (int j = 0; j < m; ++j) {
        for (int i = 0; i < k; ++i) {
            C[j*k + i] = 0.0;
            for (int l = 0; l < n; ++l) {
                C[j*k + i] += A[j*n + l] * B[l*k + i];
            }
        }
    }
}
//multiplicacion de matrices con summa
void SUMMA(MPI_Comm comm_cart, const int mb, const int nb, const int kb, double *A_loc, double *B_loc, double *C_loc) {
    int coords[2];
    MPI_Cart_coords(comm_cart, myrank, 2, coords);
    int my_col = coords[0];
    int my_row = coords[1];

    MPI_Comm row_comm, col_comm;
    int remain_dims[2] = {1, 0}; // Row communicator
    MPI_Cart_sub(comm_cart, remain_dims, &row_comm);
    remain_dims[0] = 0; remain_dims[1] = 1; // Column communicator
    MPI_Cart_sub(comm_cart, remain_dims, &col_comm);

    double *A_loc_save = (double *)calloc(mb*nb, sizeof(double));
    double *B_loc_save = (double *)calloc(nb*kb, sizeof(double));
    double *C_loc_tmp = (double *)calloc(mb*kb, sizeof(double));

    memcpy(A_loc_save, A_loc, mb*nb*sizeof(double));
    memcpy(B_loc_save, B_loc, nb*kb*sizeof(double));
    memset(C_loc, 0, mb*kb*sizeof(double));

    int nblks = N / nb;
    for (int bcast_root = 0; bcast_root < nblks; ++bcast_root) {
        if (my_col == bcast_root) memcpy(A_loc, A_loc_save, mb*nb*sizeof(double));
        MPI_Bcast(A_loc, mb*nb, MPI_DOUBLE, bcast_root, row_comm);

        if (my_row == bcast_root) memcpy(B_loc, B_loc_save, nb*kb*sizeof(double));
        MPI_Bcast(B_loc, nb*kb, MPI_DOUBLE, bcast_root, col_comm);

        matmul_naive(mb, nb, kb, A_loc, B_loc, C_loc_tmp);
        for (int i = 0; i < mb * kb; ++i) C_loc[i] += C_loc_tmp[i];
    }

    free(A_loc_save);
    free(B_loc_save);
    free(C_loc_tmp);
}
//parsea los argumentos
void parse_cmdline(int argc, char *argv[]) {
    if (argc != 2) {
        if (myrank == 0) {
            fprintf(stderr, "Uso: mpirun -np <num_procesos> ./summa <N>\n");
            MPI_Abort(MPI_COMM_WORLD, 1);
        }
    }
    N = atoi(argv[1]);
}
//main
int main(int argc, char *argv[]) {
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    parse_cmdline(argc, argv);

    int nprocs;
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
    int n_proc_rows = sqrt(nprocs);

    if (n_proc_rows * n_proc_rows != nprocs) {
        fprintf(stderr, "ERROR: El número de procesos debe ser un cuadrado perfecto.\n");
        MPI_Abort(MPI_COMM_WORLD, 1);
    }

    int ndims = 2;
    int dims[2] = {n_proc_rows, n_proc_rows};
    int periods[2] = {0, 0};
    int reorder = 0;
    MPI_Comm comm_cart;
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, periods, reorder, &comm_cart);

    int mb = N / n_proc_rows;
    double *A_loc = (double *)calloc(mb * mb, sizeof(double));
    double *B_loc = (double *)calloc(mb * mb, sizeof(double));
    double *C_loc = (double *)calloc(mb * mb, sizeof(double));

    double *C_glob = NULL;
    if (myrank == 0) C_glob = (double *)calloc(N * N, sizeof(double));

    init_matrix(A_loc, mb, mb);
    init_matrix(B_loc, mb, mb);
    print_matrix(N, N, A_loc);
    print_matrix(N, N, B_loc);

    double tstart = MPI_Wtime();
    SUMMA(comm_cart, mb, mb, mb, A_loc, B_loc, C_loc);
    double tend = MPI_Wtime();

    MPI_Gather(C_loc, mb * mb, MPI_DOUBLE, C_glob, mb * mb, MPI_DOUBLE, 0, MPI_COMM_WORLD);

    if (myrank == 0) {
        printf("La matriz resultado es:\n");
        print_matrix(N, N, C_glob);
        printf("Tiempo de ejecución de SUMMA: %f segundos\n", tend - tstart);
        free(C_glob);
    }

    free(A_loc);
    free(B_loc);
    free(C_loc);
    MPI_Finalize();
    return 0;
}
