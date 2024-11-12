#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {
    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verificación del número de procesos (simplificado para cualquier cantidad)
    if (size < 2) {
        if (rank == 0) fprintf(stderr, "Se requieren al menos 2 procesos.\n");
        MPI_Finalize();
        return 0;
    }

    // Inicialización de los vecinos (simplificada)
    int parent = (rank - 1) / 2; // En un árbol binario, el "padre" de cualquier nodo i es (i-1)/2
    int left_child = 2 * rank + 1;
    int right_child = 2 * rank + 2;

    // Mensaje inicial y contador
    int mensaje[12];
    int count = 0;
    MPI_Status status;

    // Nodo raíz inicia el mensaje
    if (rank == 0) {
        mensaje[count++] = rank; // Agrega su propio rank
    } else {
        // Recibe el mensaje del nodo padre
        MPI_Recv(mensaje, 12, MPI_INT, parent, 0, MPI_COMM_WORLD, &status);
        MPI_Get_count(&status, MPI_INT, &count);
        mensaje[count++] = rank; // Agrega su rank al mensaje
    }

    // Envía el mensaje a sus hijos, si existen en el rango de procesos
    if (left_child < size) {
        MPI_Send(mensaje, count, MPI_INT, left_child, 0, MPI_COMM_WORLD);
    }
    if (right_child < size) {
        MPI_Send(mensaje, count, MPI_INT, right_child, 0, MPI_COMM_WORLD);
    }

    // Si es un nodo hoja, imprime el camino
    if (left_child >= size && right_child >= size) {
        printf("Nodo hoja %d: Camino = ", rank);
        for (int j = 0; j < count; j++) printf("%d ", mensaje[j]);
        printf("\n");
    }

    // Finalización de MPI
    MPI_Finalize();
    return 0;
}