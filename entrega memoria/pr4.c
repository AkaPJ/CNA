#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <omp.h>
#include "ctimer.h"

#define Cseq(i,j) Cseq[(i)+m*(j)]
#define Copt(i,j) Copt[(i)+m*(j)]
#define Cbloq(i,j) Cbloq[(i)+m*(j)]
#define Cpar(i,j) Cpar[(i)+m*(j)]
#define A(i,j) A[(i)+m*(j)]
#define B(i,j) B[(i)+k*(j)]

void print(const char *matriz, const int m, const int n, double *M) {
    printf("%s = [\n", matriz);
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%9.4f", M[i + m * j]);
        }
        printf("\n");
    }
    printf("]\n");
}

int main(int argc, char *argv[]) {

    if (argc < 7) {
        printf("Usage: <m> <n> <k> <seq> <mc> <nc> <kc>\n");
        return 1;
    }

    int m, n, k, seq, mc, nc, kc;
    sscanf(argv[1], "%d", &m);
    sscanf(argv[2], "%d", &n);
    sscanf(argv[3], "%d", &k);
    sscanf(argv[4], "%d", &seq);
    sscanf(argv[5], "%d", &mc);
    sscanf(argv[6], "%d", &nc);
    sscanf(argv[7], "%d", &kc);

    double *A = (double *)malloc(m * k * sizeof(double));
    double *B = (double *)malloc(k * n * sizeof(double));
    double *Cseq = (double *)malloc(m * n * sizeof(double));
    double *Copt = (double *)malloc(m * n * sizeof(double));
    double *Cbloq = (double *)malloc(m * n * sizeof(double));
    double *Cpar = (double *)malloc(m * n * sizeof(double));

    // Inicializar matrices
    for (int i = 0; i < m * k; i++) {
        A[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
    for (int i = 0; i < k * n; i++) {
        B[i] = (double)rand() / RAND_MAX * 2.0 - 1.0;
    }
    for (int i = 0; i < m * n; i++) {
        Cseq[i] = Copt[i] = Cbloq[i] = Cpar[i] = 0;
    }

    if (seq) {
        /*************************/
        /* # Código de multiplicación secuencial                                */
        /*************************/
        double elapsed, ucpu, scpu;
        ctimer(&elapsed, &ucpu, &scpu);
        for (int j = 0; j < n; j++) {
            for (int p = 0; p < k; p++) {
                for (int i = 0; i < m; i++) {
                    Cseq(i, j) = Cseq(i, j) + A(i, p) * B(p, j);
                }
            }
        }
        ctimer(&elapsed, &ucpu, &scpu);
        printf("Tiempo = %.1f segundos de la multiplicación secuencial\n", elapsed);
    }

    /*************************/
    /* # Código de multiplicación optimizada                                */
    /*************************/
    double elapsed1, ucpu1, scpu1;
    ctimer(&elapsed1, &ucpu1, &scpu1);

    for (int j = 0; j < n; j++) {
        for (int p = 0; p < k; p++) {
            for (int i = 0; i < m; i++) {
                Copt(i, j) += A(i, p) * B(p, j);
            }
        }
    }

    ctimer(&elapsed1, &ucpu1, &scpu1);
    printf("Tiempo = %.1f segundos de la multiplicación optimizada\n", elapsed1);

    /*************************/
    /* # Código de multiplicación por bloques con empaquetado               */
    /*************************/
    double elapsed3, ucpu3, scpu3;
    ctimer(&elapsed3, &ucpu3, &scpu3);

    for (int jj = 0; jj < n; jj += nc) {
        for (int pp = 0; pp < k; pp += kc) {
            for (int ii = 0; ii < m; ii += mc) {
                for (int j = jj; j < jj + nc && j < n; j++) {
                    for (int p = pp; p < pp + kc && p < k; p++) {
                        for (int i = ii; i < ii + mc && i < m; i++) {
                            Cbloq(i, j) += A(i, p) * B(p, j);
                        }
                    }
                }
            }
        }
    }

    ctimer(&elapsed3, &ucpu3, &scpu3);
    printf("Tiempo = %.1f segundos de la multiplicación por bloques con empaquetado\n", elapsed3);

    /*************************/
    /* # Código de multiplicación paralela utilizando OpenMP                */
    /*************************/
    double elapsed2, ucpu2, scpu2;
    ctimer(&elapsed2, &ucpu2, &scpu2);

    #pragma omp parallel for
    for (int j = 0; j < n; j++) {
        for (int p = 0; p < k; p++) {
            for (int i = 0; i < m; i++) {
                Cpar(i, j) += A(i, p) * B(p, j);
            }
        }
    }

    ctimer(&elapsed2, &ucpu2, &scpu2);
    printf("Tiempo = %.1f segundos de la multiplicación paralela con OpenMP\n", elapsed2);

    if (seq) {
        /* Comprobación del resultado */
        double error_seq = 0.0;
        double error_opt = 0.0;
        double error_bloq = 0.0;
        double error_par = 0.0;

        for (int i = 0; i < m * n; i++) {
            error_seq += fabs(Cseq[i] - Cseq[i]);
            error_opt += fabs(Copt[i] - Cseq[i]);
            error_bloq += fabs(Cbloq[i] - Cseq[i]);
            error_par += fabs(Cpar[i] - Cseq[i]);
        }

        printf("Error (secuencial) = %.2e\n", error_seq);
        printf("Error (optimizada) = %.2e\n", error_opt);
        printf("Error (bloques) = %.2e\n", error_bloq);
        printf("Error (paralela) = %.2e\n", error_par);
    }

    free(A);
    free(B);
    free(Cseq);
    free(Copt);
    free(Cbloq);
    free(Cpar);

    return 0;
}