#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>
#include <cusolverSp.h>

#include "constantes2.h"
#include "herramientas2.h"

int main(int argc, char const *argv[])
{
    // *************************************************************
    // Declaracion de variables
    int ii, jj; //, kk;
    
    // Variables de tamano de la barra
    double t, deltax;

    // Parametros fisicos del problema
    double cond_ter, start_temp, end_temp;

    // Arreglos para incognitas, termino fuente y posicion
    double *temper;
    double *xx;

    double *A, *b;
    double *resultados;
    
    // *************************************************************
    
    temper = allocate_memory_vector(ni);

    xx = allocate_memory_vector(ni);
    
    // ************************************************************

    int problem_dimension = ni - 2;
    int A_size = problem_dimension * problem_dimension;
    int b_size = problem_dimension;
    int resultados_size = problem_dimension;

    A = allocate_memory_vector(A_size);
    b = allocate_memory_vector(b_size);
    resultados = allocate_memory_vector(resultados_size);

    // ************************************************************

    /*
        * Se crea la barra 1D
    */
    t = 10.0;

    for (ii = 0; ii <= ni - 1; ii++)
    {
        double dii1 = (double)ii;
        double dni1 = (double)(ni - 1);
        xx[ii] = t * dii1 / dni1;
    }
    
    double dni1 = (double)(ni - 1);
    deltax = 1.0 / dni1;

    /*
    * Se definen los parametros fisicos del problema
    */
    cond_ter = 100.0;
    start_temp = 0;
    end_temp = 10;

    /*
    * Inicializacion de los arreglos a utilizar
    */

    initialise_array_with_value(temper, ni, 0.0);
    temper[0] = start_temp;
    temper[ni-1] = end_temp;
    initialise_array_with_value(resultados, resultados_size, 0.0);

    // *************************************************************************
    int rows_num = problem_dimension;
    int columns_num = problem_dimension;

    get_1D_heat_equation_matrix(A, rows_num, columns_num, cond_ter);
    get_vector_independent_terms(b, b_size, cond_ter, start_temp, end_temp);

    // ****************************************************************************
    double *value;
    int *col_ind;
    int *row_ptr;
    int nnz_elements;
    int value_size, col_ind_size, row_ptr_size;

    nnz_elements = get_nonzero_elements();
    value_size = nnz_elements;
    col_ind_size = value_size;
    row_ptr_size = problem_dimension;

    value = allocate_memory_vector(value_size);
    col_ind = allocate_memory_vector_int(col_ind_size);
    row_ptr = allocate_memory_vector_int(row_ptr_size);

    get_csr_schema(cond_ter, nnz_elements, value, col_ind, row_ptr);

    print_vector(value, value_size);
    print_vector_int(col_ind, col_ind_size);
    print_vector_int(row_ptr, row_ptr_size);
    // ****************************************************************************

    /*
        * Region del solver
    */
    // int N = problem_dimension;
    // int NRHS = 1;
    // int LDA = N;
    // int LDB = NRHS;
    // int INFO_TRF = 0;
    // int INFO_TRS = 0;
    // int INFO_TRI = 0;
    // int INFO = 0;
    // int *IPIV = allocate_memory_vector_int(N);

    // // INFO_TRF = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, LDA, IPIV);
    // // INFO_TRS = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, NRHS, A, LDA, IPIV, b, LDB);
    // // INFO_TRI = LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, dense_matrix, LDA, IPIV);
    // INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, A, LDA, IPIV, b, LDB);

    // // printf("Info DGETRF = %d\n", INFO_TRF);
    // // printf("Info DGETRS = %d\n", INFO_TRS);
    // // printf("Info DGETRI = %d\n", INFO_TRI);
    // printf("Info DGESV = %d\n\n", INFO);

    // get_temper_array(temper, b, start_temp, end_temp);

    // ****************************************************************************
    double TOL = 1e-5;
    int REORDER = 3;
    int singularity = 0;
    double *u;
    int u_size = b_size;

    u = allocate_memory_vector(u_size);

    cusolverSpHandle_t handle;
    cusparseMatDescr_t descrMat;
    cusolverStatus_t status;

    cusolverSpCreate(&handle);
    cusparseCreateMatDescr(&descrMat);
    cusparseSetMatType(descrMat, CUSPARSE_INDEX_BASE_ZERO);
    cusparseSetMatIndexBase(descrMat, CUSPARSE_INDEX_BASE_ZERO);

    #pragma acc data copyin(b[:b_size], value[:value_size], col_ind[:col_ind_size], row_ptr[:row_ptr_size]) copyout(u[:u_size])
    {
        status = cusolverSpDcsrlsvchol(handle, A_size, nnz_elements, descrMat, value, row_ptr, col_ind, b, TOL, REORDER, u, &singularity);
    }
    printf("status = %d\n", status);

    // ****************************************************************************

    // print_flat_matrix(A, rows_num, columns_num);
    // print_vector(b, b_size);
    // print_vector(temper, ni);
    
    // ****************************************************************************
    /*
    * Escritura de resultados
    */
    // Abrir el archivo para escribir
    FILE *file = fopen("eqHeat1D.101", "w");

    if (file != NULL) {
        // Utilizar un bucle para imprimir y guardar los datos
        for (ii = 0; ii < ni; ii++) {

            fprintf(file, "%f %f\n", xx[ii], temper[ii]);

        }
        // Cerrar el archivo después de escribir
        fclose(file);
    } else {
        printf("Error al abrir el archivo.\n");
    }

    // Liberar memoria de los punteros

    return 0;
}
