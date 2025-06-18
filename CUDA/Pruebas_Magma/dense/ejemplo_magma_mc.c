#include <magma_auxiliary.h>
#include <magma_types.h>
#include <magmablas_d.h>
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include "magma_v2.h"
#include "magma_lapack.h"

#define N 20000

void generar_matriz_y_b(int m, int n, double *A, double *b);
int test_result(double *vector, int dim, double tol);

int main() {
    magma_int_t m = N, n = N, lda = m, ldb = m, nrhs = 1, info = -1;
    double *h_A, *h_B, *h_X;
    magma_int_t *ipiv;
    magmaDouble_ptr d_A, d_B;
    magma_trans_t trans = MagmaNoTrans;

    // Inicializar MAGMA
    magma_init();
    // Reservar Memoria en host
    h_A = (double *) malloc(m*n*sizeof(double));
    h_B = (double *) malloc(m*nrhs*sizeof(double));
    h_X = (double *) malloc(m*nrhs*sizeof(double));
    ipiv = (magma_int_t *) malloc(m*sizeof(magma_int_t));

    // Inicializar A y B
    generar_matriz_y_b(m, n, h_A, h_B);

    // Factorizar A
    magma_dgetrf(m, n, h_A, lda, ipiv, &info);
    printf("Info DGETRF = %d\n\n", info);

    // Reservar Memoria en device
    magma_dmalloc(&d_A, m * n);
    magma_dmalloc(&d_B, m * nrhs);

    // Copiar datos al device
    magma_dsetmatrix(m, n, h_A, lda, d_A, lda, NULL);
    magma_dsetmatrix(m, nrhs, h_B, ldb, d_B, ldb, NULL);

    // Resolver el sistema usando la matriz factorizada
    magma_dgetrs_gpu(trans, n, nrhs, d_A, lda, ipiv, d_B, ldb, &info);
    printf("Info DGETRS = %d\n\n", info);

    magma_dgetmatrix(m, n, d_A, lda, h_A, lda, NULL);
    magma_dgetmatrix(m, nrhs, d_B, ldb, h_X, ldb, NULL);

    // printf("Matriz A factorizada con LU:\n");
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         printf("%5.2f ", h_A[j * n + i]);
    //     }
    //     printf("\n");
    // }

    // printf("Vector solucion:\n");
    // for (int i = 0; i < n; i++) {
    //     printf("%5.1f\n", h_X[i]);
    // }

    test_result(h_X, n, 0.000001);

    free(h_A);
    free(h_B);
    free(h_X);
    free(ipiv);
    magma_free(d_A);
    magma_free(d_B);
    // Finalizar MAGMA
    magma_finalize();
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

int test_result(double *vector, int dim, double tol) {
    for (int i; i < dim; i++) {
        if ((1.0 - vector[i]) > tol) {
            printf("Test failure!");
            return -1;
        }
    }
    printf("Test passed!");
    return 0;
}
