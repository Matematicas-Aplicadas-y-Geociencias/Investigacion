#include <stdio.h>
#include <cusparse.h>
#include <cusolverSp.h>

#include "herramientas.h"
#include "sistema_de_prueba.h"

#define VALUE 256

const int ROWS = VALUE;
const int COLUMNS = VALUE;
const int VECTOR_SIZE = VALUE;

int main(int argc, char const *argv[])
{
    double *A = NULL;
    double *b = NULL;
    double *csrValA = NULL;
    int *csrRowPtrA = NULL;
    int *csrColIndA = NULL;

    int N = ROWS;
    int NNZA;
    double TOL = 1e-6;
    int REORDER = 0;
    double *x = NULL;
    int singularity;

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
    double *d_b = NULL;
    double *d_csrValA = NULL;
    int *d_csrRowPtrA = NULL;
    int *d_csrColIndA = NULL;
    double *d_x = NULL;

    d_b = allocateMemoryVectorDevice(VECTOR_SIZE);
    d_x = allocateMemoryVectorDevice(VECTOR_SIZE);
    d_csrValA = allocateMemoryVectorDevice(NNZA);
    d_csrColIndA = allocateMemoryVectorDeviceInteger(NNZA);
    d_csrRowPtrA = allocateMemoryVectorDeviceInteger(ROWS + 1);

    cusolverSpCreate(&handle);
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL); cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

    transfer_vector_host_to_device(b, d_b, VECTOR_SIZE);
    transfer_vector_host_to_device(x, d_x, VECTOR_SIZE);
    transfer_vector_host_to_device(csrValA, d_csrValA, NNZA);
    transfer_vector_host_to_device_integer(csrColIndA, d_csrColIndA, NNZA);
    transfer_vector_host_to_device_integer(csrRowPtrA, d_csrRowPtrA, ROWS + 1);

    // cusolverSpDcsrlsvluHost(handle, N, NNZA, descrA, csrValA, csrRowPtrA, csrColIndA, b, TOL, REORDER, x, &singularity);
    cusolverSpDcsrlsvqr(handle, N, NNZA, descrA, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_b, TOL, REORDER, d_x, &singularity);
    // cusolverSpDcsrlsvchol(handle,N, NNZA, descrA, d_csrValA, d_csrRowPtrA, d_csrColIndA, d_b, TOL, REORDER, d_x, &singularity);

    transfer_vector_device_to_host(d_x, x, VECTOR_SIZE);

    cudaDeviceSynchronize();

    printf("n = %d\n", singularity);
    print_vector(x, VECTOR_SIZE);

    cudaFree(d_csrRowPtrA);
    cudaFree(d_csrColIndA);
    cudaFree(d_csrValA);
    cudaFree(d_x);
    cudaFree(d_b);

    free(x);
    free(csrColIndA);
    free(csrRowPtrA);
    free(csrValA);
    free(b);
    free(A);

    return 0;
}
