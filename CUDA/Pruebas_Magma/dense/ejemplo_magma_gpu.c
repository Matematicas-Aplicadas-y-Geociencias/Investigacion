#include <math.h>
#include <stdio.h>
#include <stdbool.h>
#include <magma_v2.h>
#include <magma_lapack.h>
#include <magma_auxiliary.h>
#include <magma_types.h>
#include <magmablas_d.h>

#define N 10000

void generar_matriz_y_b(int m, int n, double *A, double *b);
int test_result(double *vector, int dim, double tol);

int main() {
    magma_int_t m = N, n = N, lda = m, ldb = m, nrhs = 1, info = -1;
    double *h_A, *h_B, *h_X;
    magma_int_t *ipiv;
    magmaDouble_ptr d_A, d_B;
    magma_trans_t trans = MagmaNoTrans;
    magma_queue_t queue;
    magma_queue_create(0, &queue);

    // Inicializar MAGMA
    magma_init();
    magma_dmalloc_cpu(&h_A, m * n);
    magma_dmalloc_cpu(&h_B, m * nrhs);
    magma_dmalloc_cpu(&h_X, m * nrhs);
    magma_imalloc_cpu(&ipiv, m);

    // Inicializar A y B
    generar_matriz_y_b(m, n, h_A, h_B);

    // Reservar Memoria en device
    magma_dmalloc(&d_A, m * n);
    magma_dmalloc(&d_B, m * nrhs);

    // Copiar datos al device
    magma_dsetmatrix(m, n, h_A, lda, d_A, lda, queue);
    magma_dsetmatrix(m, nrhs, h_B, ldb, d_B, ldb, queue);

    // Factorizar A
    magma_dgetrf_gpu(m, n, d_A, lda, ipiv, &info);
    printf("Info DGETRF = %d\n\n", info);
    // Resolver el sistema usando la matriz factorizada
    magma_dgetrs_gpu(trans, n, nrhs, d_A, lda, ipiv, d_B, ldb, &info);
    printf("Info DGETRS = %d\n\n", info);

    magma_dgetmatrix(m, n, d_A, lda, h_A, lda, queue);
    magma_dgetmatrix(m, nrhs, d_B, ldb, h_X, ldb, queue);

    // printf("Matriz A factorizada con LU:\n");
    // for (int i = 0; i < m; i++) {
    //     for (int j = 0; j < n; j++) {
    //         printf("%5.2f ", h_A[j * n + i]);
    //     }
    //     printf("\n");
    // }

    // printf("Vector solucion:\n");
    // for (int i = 0; i < n; i++) {
    //     printf("%d %5.1f\n", i, h_X[i]);
    // }

    test_result(h_X, n, 0.000001);

    magma_queue_destroy(queue);
    magma_free_cpu(h_A);
    magma_free_cpu(h_B);
    magma_free_cpu(h_X);
    magma_free_cpu(ipiv);
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
    for (int i = 0; i < dim; i++) {
        if (fabs(vector[i] - 1.0) > tol) {
            printf("Test failure at index %d: got %.8f\n", i, vector[i]);
            return -1;
        }
    }
    printf("Test passed!");
    return 0;
}
