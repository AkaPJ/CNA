#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    const int ndims = 2;
    int ierr;
    int p, rank, cart_rank;
    MPI_Comm comm2D;
    int dims[ndims], coord[ndims];
    int wrap_around[ndims];
    int reorder, nrows, ncols;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (argc < 3) {
        if (rank == 0) printf("Usage: mpirun -np %d %s nrows ncols\n", p, argv[0]);
        MPI_Finalize();
        return 0;
    }

    sscanf(argv[1], "%d", &nrows);
    sscanf(argv[2], "%d", &ncols);
    if (nrows * ncols > p) {
        if (rank == 0) printf("%d debe ser >= %d\n", p, nrows * ncols);
        MPI_Finalize();
        return 0;
    }

    dims[0] = nrows;
    dims[1] = ncols;
    wrap_around[0] = 1;
    wrap_around[1] = 1;
    reorder = 1;

    MPI_Dims_create(p, ndims, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);
    MPI_Cart_coords(comm2D, rank, ndims, coord);
    MPI_Cart_rank(comm2D, coord, &cart_rank);

    // Asignar valores fijos de las matrices A y B
    int A = 4 * coord[0] + coord[1] + 1;
    int B = 4 * coord[0] + coord[1];

    // Variables para el algoritmo de Cannon
    int C = 0;
    int temp_A, temp_B;
    MPI_Status status;

    // Desplazamiento inicial
    int src, dest;
    MPI_Cart_shift(comm2D, 1, -coord[0], &src, &dest);
    MPI_Sendrecv_replace(&A, 1, MPI_INT, dest, 0, src, 0, comm2D, &status);

    MPI_Cart_shift(comm2D, 0, -coord[1], &src, &dest);
    MPI_Sendrecv_replace(&B, 1, MPI_INT, dest, 0, src, 0, comm2D, &status);

    // Fases del algoritmo de Cannon
    for (int k = 0; k < dims[0]; k++) {
        C += A * B;

        MPI_Cart_shift(comm2D, 1, -1, &src, &dest);
        MPI_Sendrecv_replace(&A, 1, MPI_INT, dest, 0, src, 0, comm2D, &status);

        MPI_Cart_shift(comm2D, 0, -1, &src, &dest);
        MPI_Sendrecv_replace(&B, 1, MPI_INT, dest, 0, src, 0, comm2D, &status);
    }

    printf("Process %d: Cartesian coordinates = (%d, %d), C = %d\n", rank, coord[0], coord[1], C);

    MPI_Comm_free(&comm2D);
    MPI_Finalize();
    return 0;
}