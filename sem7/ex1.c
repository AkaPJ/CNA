#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#define A(i,j) A[(i)*k+(j)]
#define B(i,j) B[(i)*n+(j)]
#define C(i,j) C[(i)*n+(j)]

void ctimer(double *elapsed, double *ucpu, double *scpu) {
    // Implementación de la función ctimer
    // Esta función debe medir el tiempo transcurrido y el tiempo de CPU
}

void microkernel(int mc, int nc, int kc, double *Ac, double *Bc, double *Cc) {
    for (int j = 0; j < nc; j++) {
        for (int p = 0; p < kc; p++) {
            for (int i = 0; i < mc; i++) {
                Cc[i + j * mc] += Ac[i + p * mc] * Bc[p + j * kc];
            }
        }
    }
}

void macrokernel(int m, int n, int k, int mc, int nc, int kc, double *A, double *B, double *C) {
    double *Ac = (double *)malloc(mc * kc * sizeof(double));
    double *Bc = (double *)malloc(kc * nc * sizeof(double));
    double *Cc = (double *)malloc(mc * nc * sizeof(double));

    for (int jc = 0; jc < n; jc += nc) {
        for (int pc = 0; pc < k; pc += kc) {
            for (int ic = 0; ic < m; ic += mc) {
                // Copiar bloques de A y B en Ac y Bc
                for (int p = 0; p < kc; p++) {
                    for (int i = 0; i < mc; i++) {
                        Ac[i + p * mc] = A[(ic + i) * k + (pc + p)];
                    }
                }
                for (int j = 0; j < nc; j++) {
                    for (int p = 0; p < kc; p++) {
                        Bc[p + j * kc] = B[(pc + p) * n + (jc + j)];
                    }
                }
                // Inicializar Cc
                for (int j = 0; j < nc; j++) {
                    for (int i = 0; i < mc; i++) {
                        Cc[i + j * mc] = 0.0;
                    }
                }
                // Llamar al microkernel
                microkernel(mc, nc, kc, Ac, Bc, Cc);
                // Acumular resultados en C
                for (int j = 0; j < nc; j++) {
                    for (int i = 0; i < mc; i++) {
                        C[(ic + i) * n + (jc + j)] += Cc[i + j * mc];
                    }
                }
            }
        }
    }

    free(Ac);
    free(Bc);
    free(Cc);
}

int main() {
    int m, n, k, mc, nc, kc;
    printf("Pon tamaños de m, n, k: ");
    scanf("%d %d %d", &m, &n, &k);
    printf("Pon tamaños de bloque de mc, nc, kc: ");
    scanf("%d %d %d", &mc, &nc, &kc);

    double *A = (double *)malloc(m * k * sizeof(double));
    double *B = (double *)malloc(k * n * sizeof(double));
    double *Cseq = (double *)malloc(m * n * sizeof(double));
    double *C = (double *)malloc(m * n * sizeof(double));

    // Inicializar matrices
    srand(time(NULL));
    for (int i = 0; i < m * k; i++) {
        A[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
    for (int i = 0; i < k * n; i++) {
        B[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
    for (int i = 0; i < m * n; i++) {
        C[i] = Cseq[i] = 0.0;
    }

    // Multiplicación secuencial
    double elapsed, ucpu, scpu;
    ctimer(&elapsed, &ucpu, &scpu);
    for (int j = 0; j < n; j++) {
        for (int p = 0; p < k; p++) {
            for (int i = 0; i < m; i++) {
                Cseq[i * n + j] += A[i * k + p] * B[p * n + j];
            }
        }
    }
    ctimer(&elapsed, &ucpu, &scpu);
    printf("Tiempo = %.1f segundos de la multiplicación simple\n", elapsed);

    // Multiplicación con bloques
    ctimer(&elapsed, &ucpu, &scpu);
    macrokernel(m, n, k, mc, nc, kc, A, B, C);
    ctimer(&elapsed, &ucpu, &scpu);
    printf("Tiempo = %.1f segundos de la multiplicación con bloques\n", elapsed);

    // Comparar resultados
    double error = 0.0;
    for (int i = 0; i < m * n; i++) {
        error += fabs(C[i] - Cseq[i]);
    }
    printf("Error = %f\n", error);

    free(A);
    free(B);
    free(Cseq);
    free(C);

    return 0;
}