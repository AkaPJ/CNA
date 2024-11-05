#include <stdio.h>
#include <stdlib.h>
#include <omp.h>

int main(int argc, char* argv[]) {
    int threads = omp_get_num_threads();
    printf("Número de hilos: %d\n", threads);
}