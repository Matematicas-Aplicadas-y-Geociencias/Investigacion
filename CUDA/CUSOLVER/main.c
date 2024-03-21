#include <stdio.h>
#include <cusparse.h>
#include <cusolverSp.h>
#include <openacc.h>

#include "herramientas.h"
#include "sistema_de_prueba.h"

#define FIX_VALUE 256

const int ROWS = FIX_VALUE;
const int COLUMNS = FIX_VALUE;
const int VECTOR_SIZE = FIX_VALUE;

int main(int argc, char const *argv[])
{
    double *A = NULL;
    double *b = NULL;
    double *csrValA = NULL;
    int *csrRowPtrA = NULL;
    int *csrColIndA = NULL;
    double *x = NULL;

    int N = ROWS;
    int NNZA;
    double TOL = 1e-6;
    int REORDER = 0;
    int singularity = 4;

    A = allocate_memory_matrix(ROWS, COLUMNS);
    b = allocate_memory_vector(VECTOR_SIZE);
    x = allocate_memory_vector(VECTOR_SIZE);
    crear_matriz_de_prueba(A, ROWS, COLUMNS);
    crear_vector_de_prueba(b, VECTOR_SIZE);
    NNZA = obtener_nnzMatriz(A, ROWS, COLUMNS);

    csrValA = allocate_memory_vector(NNZA);
    csrColIndA = allocate_memory_vector_integer(NNZA);
    csrRowPtrA = allocate_memory_vector_integer(ROWS + 1);
    obtener_formato_csr(A, ROWS, COLUMNS, csrValA, csrRowPtrA, csrColIndA);

    // Variables CUDA
    cusolverSpHandle_t handle = NULL;
    cusparseMatDescr_t descrA = NULL;

    cusolverSpCreate(&handle);
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL); cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
    printf("n = %d\n", singularity);

    #pragma acc data copyin(b[:VECTOR_SIZE], csrValA[:NNZA], csrRowPtrA[:ROWS+1], csrColIndA[:NNZA]), copyout(x[:VECTOR_SIZE])
    {
        // #pragma acc kernels
        // {
            #pragma acc host_data use_device(b, csrValA, csrRowPtrA, csrColIndA, x)
            {
                // cusolverSpDcsrlsvluHost(handle, N, NNZA, descrA, csrValA, csrRowPtrA, csrColIndA, b, TOL, REORDER, x, &singularity);
                // cusolverSpDcsrlsvqr(handle, N, NNZA, descrA, csrValA, csrRowPtrA, csrColIndA, b, TOL, REORDER, x, &singularity);
                cusolverSpDcsrlsvchol(handle,N, NNZA, descrA, csrValA, csrRowPtrA, csrColIndA, b, TOL, REORDER, x, &singularity);

                // print_vector_integer(csrColIndA, NNZA);
            }
        // }
    }


    print_vector(x, VECTOR_SIZE);
    printf("n = %d\n", singularity);

    free(x);
    free(csrColIndA);
    free(csrRowPtrA);
    free(csrValA);
    free(b);
    free(A);

    return 0;
}
