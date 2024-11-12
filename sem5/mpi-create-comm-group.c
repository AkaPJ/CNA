#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int p, q; /* sqrt(p) */
    int rank, comm_rank;
    MPI_Group group_world;
    MPI_Group comm1_group, comm2_group;
    MPI_Comm comm1, comm2;
    int *comm1_procs, *comm2_procs;
    int err;
    int sum;

    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    if (p < 2) {
        if (rank == 0) printf("ERROR: The number of processes %d must be >= 2\n", p);
        MPI_Finalize();
        return 0;
    }

    // Dividir los procesos en dos grupos: pares e impares
    int comm1_size = (p + 1) / 2; // Número de procesos en el primer grupo (aproximadamente la mitad)
    int comm2_size = p / 2; // Número de procesos en el segundo grupo (aproximadamente la otra mitad)
    comm1_procs = (int *)malloc(comm1_size * sizeof(int));
    comm2_procs = (int *)malloc(comm2_size * sizeof(int));

    int i, j = 0, k = 0;
    for (i = 0; i < p; i++) {
        if (i % 2 == 0)
            comm1_procs[j++] = i;
        else
            comm2_procs[k++] = i;
    }

    // Obtener el grupo de procesos de MPI_COMM_WORLD
    MPI_Comm_group(MPI_COMM_WORLD, &group_world);

    // Crear los grupos de procesos
    MPI_Group_incl(group_world, comm1_size, comm1_procs, &comm1_group);
    MPI_Group_incl(group_world, comm2_size, comm2_procs, &comm2_group);

    // Crear los comunicadores asociados a los grupos
    MPI_Comm_create(MPI_COMM_WORLD, comm1_group, &comm1);
    MPI_Comm_create(MPI_COMM_WORLD, comm2_group, &comm2);

    // Realizar operaciones dentro de los grupos
    MPI_Comm comm;
    if (rank % 2 == 0) {
        comm = comm1;
    } else {
        comm = comm2;
    }

    if (comm != MPI_COMM_NULL) {
        MPI_Comm_rank(comm, &comm_rank);
        int value = 1;
        MPI_Bcast(&value, 1, MPI_INT, 0, comm);
        MPI_Reduce(&value, &sum, 1, MPI_INT, MPI_SUM, 0, comm);
        if (comm_rank == 0) {
            printf("Grupo %d, suma: %d\n", (rank % 2 == 0) ? 1 : 2, sum);
        }
    }

    free(comm1_procs);
    free(comm2_procs);

    MPI_Finalize();
    return 0;
}