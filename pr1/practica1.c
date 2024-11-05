#include <stdlib.h>
#include <stdio.h>

//cada posicio consecutiva de memoria pertanyen a la mateixa columna

void generaMatriz(double *m, int n) {
    for (int i = 0; i < n; i++) {
        m[i] = (double)rand() / RAND_MAX;
    }
}

void multseq(double *m, double *n, double *res, int num) {
    for (int i = 0; i < num; i++) {
        res[i] = 0;
        for (int j = 0; j < num; j++) {
            res[i] += m[j] * n[j];
        }
    }
}

void muestraMatriz(double *m, int n) {
    for (int i = 0; i < n; i++) {
        printf("%f ", m[i]);
    }
    printf("\n");
}

int main(int argc, char *argv[]) {
    int n = 4;
    double *m = (double *)malloc(n * sizeof(double));
    double *num = (double *)malloc(n * sizeof(double));
    double *res = (double *)malloc(n * sizeof(double));
    generaMatriz(m, n);
    generaMatriz(num, n);
    multseq(m, num, res, n);
    muestraMatriz(res, n);
    return 0;
}