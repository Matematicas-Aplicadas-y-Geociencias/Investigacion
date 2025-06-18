#include <stdio.h>
#include <stdlib.h>
#include "lapacke.h"

#define N 4

void generar_matriz_y_b(int m, int n, double *A, double *b);

int main() {
    int m = N, n = N, lda = m, ldb = 1, nrhs = 1;
    int info = -1;
    double *A, *B, *X;
    int *ipiv;
    char trans = 'N';

    // Reservar Memoria
    A = (double *) malloc(m*n*sizeof(double));
    B = (double *) malloc(m*nrhs*sizeof(double));
    X = (double *) malloc(m*nrhs*sizeof(double));
    ipiv = (int *) malloc(m*sizeof(int));

    // Inicializar A y B
    generar_matriz_y_b(m, n, A, B);
    // Factorizar A
    info = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, m, n, A, lda, ipiv);
    printf("Info DGETRF = %d\n\n", info);
    // Resolver el sistema usando la matriz factorizada
    info = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, trans, n, nrhs, A, lda, ipiv, B, ldb);
    printf("Info DGETRS = %d\n\n", info);

    printf("Matriz A:\n");
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            printf("%5.2f ", A[i * n + j]);
        }
        printf("\n");
    }

    printf("Vector b:\n");
    for (int i = 0; i < n; i++) {
        printf("%5.1f\n", B[i]);
    }

    free(A);
    free(B);
    free(X);
    free(ipiv);
    return 0;
}

void generar_matriz_y_b(int m, int n, double *A, double *b) {
    for (int i = 0; i < m; i++) {
        b[i] = 0.0;
        for (int j = 0; j < n; j++) {
            if (i == j) {
                A[i * N + j] = 2.0;
            } else {
                A[i * N + j] = 1.0;
            }
            b[i] += A[i * N + j];
        }
    }
}
