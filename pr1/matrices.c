#include <stdio.h>
#include <stdlib.h>
//#include <mpi.h>

void generaMatrizCuadrada(double* m, int n) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        m[i * n + j] = (double)rand() / RAND_MAX;
    }
}

void muestraMatriz(int n, double* m) {
    for (int i = 0; i < n; i++) {
        for (int j = 0; j < n; j++)
        printf("%f ", m[i * n + j]);
        printf("\n");
    }
}

int main(int argc, char* argv[]) {
    int n = 4;
    double* m = (double*)malloc(n * n * sizeof(double));
    generaMatrizCuadrada(m, n);
    muestraMatriz(n, m);
    return 0;
}