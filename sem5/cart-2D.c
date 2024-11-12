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
    wrap_around[0] = 0;
    wrap_around[1] = 0;
    reorder = 1;

    MPI_Dims_create(p, ndims, dims);
    MPI_Cart_create(MPI_COMM_WORLD, ndims, dims, wrap_around, reorder, &comm2D);
    MPI_Cart_coords(comm2D, rank, ndims, coord);
    MPI_Cart_rank(comm2D, coord, &cart_rank);

    if (comm2D != MPI_COMM_NULL) {
        printf("Process %d: Cartesian coordinates = (%d, %d)\n", rank, coord[0], coord[1]);
    }

    // Crear sub-mallas (1-D) para filas y columnas
    MPI_Comm row_comm, col_comm;

    // Crear comunicador para las filas
    int remain_dims_rows[2];
    remain_dims_rows[0] = 0; // mantener la dimensión de las filas
    remain_dims_rows[1] = 1; // mantener la dimensión de las columnas
    MPI_Cart_sub(comm2D, remain_dims_rows, &row_comm);

    // Crear comunicador para las columnas
    int remain_dims_cols[2];
    remain_dims_cols[0] = 1; // mantener la dimensión de las filas
    remain_dims_cols[1] = 0; // mantener la dimensión de las columnas
    MPI_Cart_sub(comm2D, remain_dims_cols, &col_comm);

    // Imprimir coordenadas en la malla 1-D para filas y columnas
    int row_rank, col_rank;
    MPI_Comm_rank(row_comm, &row_rank);
    MPI_Comm_rank(col_comm, &col_rank);

    printf("Process %d: Cartesian coordinates = (%d, %d), Row Rank = %d, Column Rank = %d\n",
            rank, coord[0], coord[1], row_rank, col_rank);

    MPI_Comm_free(&comm2D);
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    MPI_Finalize();
    return 0;
}