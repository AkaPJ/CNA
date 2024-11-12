#include <stdio.h>
#include <mpi.h>

int main(int argc, char *argv[]) {

    int rank, size;
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    // Verificación del número de procesos
    if (size != 12) {
        if (!rank) fprintf(stderr, "%s: Se debe ejecutar con 12 procesos.\n", argv[0]);
        MPI_Finalize();
        return 0;
    }

    // Definición del grafo según la figura en el PDF
    int nodes = 12;
    int index[] = {3, 5, 7, 8, 10, 12, 13, 14, 15, 16, 17, 18}; // Índice acumulado de vecinos por nodo
    int edges[] = {1, 2, 4, 0, 3, 5, 0, 6, 1, 7, 8, 4, 9, 2, 10, 11, 3, 4}; // Vecinos de cada nodo

    // Creación de la topología de grafo
    MPI_Comm MPI_COMM_GRAPH;
    MPI_Graph_create(MPI_COMM_WORLD, nodes, index, edges, 0, &MPI_COMM_GRAPH);

    // Obtención del número de vecinos y de los vecinos mismos
    int nneighbors;
    const int max_neighbors = 3;
    int neighbors[max_neighbors];
    MPI_Graph_neighbors_count(MPI_COMM_GRAPH, rank, &nneighbors);
    MPI_Graph_neighbors(MPI_COMM_GRAPH, rank, max_neighbors, neighbors);

    // Inicialización del mensaje y del contador de elementos en el mensaje
    int mensaje[3];
    int count = 0;
    MPI_Status status;
    int i;

    // Nodo raíz inicia el mensaje
    if (rank == 0) {
        mensaje[count++] = rank; // Agregamos el rank 0 al mensaje
    }

    // Comunicación con vecinos
    for (i = 0; i < nneighbors; i++) {
        if (neighbors[i] < rank) { // Nodo "padre"
            // Recibe el mensaje del nodo padre
            MPI_Recv(mensaje, 3, MPI_INT, neighbors[i], 0, MPI_COMM_GRAPH, &status);
            MPI_Get_count(&status, MPI_INT, &count); // Obtener el tamaño del mensaje recibido
        } else { // Nodo "hijo"
            // Añade su rango al mensaje antes de enviarlo a los hijos
            mensaje[count++] = rank;
            MPI_Send(mensaje, count, MPI_INT, neighbors[i], 0, MPI_COMM_GRAPH);
        }
    }

    // Si es un nodo hoja (solo tiene un vecino y es mayor que su vecino), imprime el camino
    if (nneighbors == 1 && neighbors[0] < rank) {
        printf("Nodo hoja %d: Camino = ", rank);
        for (int j = 0; j < count; j++) printf("%d ", mensaje[j]);
        printf("\n");
    }

    // Finalización de MPI
    MPI_Finalize();
    return 0;
}
