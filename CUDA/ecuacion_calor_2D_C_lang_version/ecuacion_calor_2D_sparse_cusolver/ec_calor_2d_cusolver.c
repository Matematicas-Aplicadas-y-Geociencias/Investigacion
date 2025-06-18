#include <stdio.h>
#include <stdlib.h>
#include <cusparse.h>
#include <cusolverSp.h>

#include "constantes.h"
#include "herramientas.h"

int main(int argc, char const *argv[])
{
    // ************************************************************************************

    // ****** VARIABLE DECLARATION: PRIMITIVE ******

    // LOOP ITERATION INDEX
    int ii, jj; //, kk;

    // VARIABLES FOR THE SIZE OF THE PLATE
    double a, b, delta_x, delta_y;

    // PHYSICAL PARAMETERS OF THE PROBLEM
    double cond_ter, temp_ini, temp_fin;
    double flux_aba, flux_arr; //, alpha;

    // CSR parameters
    int nonzero_elements;
    int rows_number;
    int columns_number;
    int csr_row_pointer_size;

    // CUSOLVER: HIGHT LEVEL FUNCTION PARAMETERS
    double TOLERANCE;
    int REORDER;
    int SINGULARITY;

    // ****** VARIABLE DECLARATION: POINTERS ******
    // ARRAYS FOR POSITION
    double *xx;
    double *yy;

    // ARRAYS FOR VECTOR OF INDEPENDENT TERMS AND VECTOR OF RESULTS
    double *vector_b;
    double *results;

    // ARRAYS FOR COMPRESSED SPARSE ROW STORE (CSR)
    double *csr_values;
    int *csr_column_index;
    int *csr_row_pointer;

    // ****** VARIABLE DECLARATION: POINTERS TO POINTERS ******
    // ARRAYS FOR SOURCE TERM
    double **temper;

    // BIDIMENSIONAL ARRAYS FOR NONZERO COEFFICIENTS
    double **AI, **AD, **AC;
    double **BI, **BD;

    // ****** VARIABLE DECLARATION: CUSPARSE AND CUSOLVER ******
    // HELPER FUNCTIONS PARAMETERS
    cusparseMatDescr_t descrA = NULL;
    cusolverSpHandle_t handle = NULL;

    // ERROR STATUS PARAMETERS
    cusparseStatus_t cusparse_status;
    cusolverStatus_t cusolver_status;

    // ************************************************************************************

    // 2D MESH PARAMETERS
    a = 10.0;
    b = 5.0;

    // CALCULATE DELTA VALUES
    double dmi1 = (double)(mi - 1);
    delta_x = 1.0 / dmi1;

    double dnj1 = (double)(nj - 1);
    delta_y = 1.0 / dnj1;

    // THE PHYSICAL PARAMETERS OF THE PROBLEM ARE DEFINED
    cond_ter = 100.0;
    temp_ini = 308.0;
    temp_fin = 298.0;
    flux_aba = 308.0;
    flux_arr = 308.0;
    // alpha    = 0.5; //parametro de relajacion

    // INITIALISATION OF CUSOLVER HIGHT LEVEL FUNCTION PARAMETERS
    TOLERANCE = 1e-5;
    REORDER = 3;
    SINGULARITY = 0;

    // INITIALISATION SECTION OF CSR STORAGE PARAMETERS
    nonzero_elements = get_nonzero_elements();
    rows_number = (mi - 2) * (nj - 2);
    columns_number = (mi - 2) * (nj - 2);
    csr_row_pointer_size = rows_number + 1;

    // INITIALISATION OF THE PARAMETERS OF THE HELPER FUNCTIONS CUSPARSE AND CUSOLVER
    cusolver_status = cusolverSpCreate(&handle);
    print_info_cusolver("Handle", cusolver_status);

    cusparse_status = cusparseCreateMatDescr(&descrA);
    print_info_cusparse("Matrix descriptor", cusparse_status);

    /* This parameters are by default (read documentation: https://docs.nvidia.com/ cuda/cusparse/index.html#cusparse-types-reference):
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);
    */

    // ************************************************************************************

    // MEMORY ALLOCATION SECTION

    // FOR POSITIONS
    xx = allocate_memory_vector_double(mi);
    yy = allocate_memory_vector_double(nj);

    // FOR CSR STORAGES ARRAYS
    csr_values = allocate_memory_vector_double(nonzero_elements);
    csr_column_index = allocate_memory_vector_int(nonzero_elements);
    csr_row_pointer = allocate_memory_vector_int(csr_row_pointer_size);
    vector_b = allocate_memory_vector_double(columns_number);
    results = allocate_memory_vector_double(columns_number);

    // BIDIMENSIONAL ARRAYS
    temper = allocate_memory_bidimensional_matrix_double(mi, nj);
    BI = allocate_memory_bidimensional_matrix_double(nj, mi);
    AI = allocate_memory_bidimensional_matrix_double(mi, nj);
    AD = allocate_memory_bidimensional_matrix_double(mi, nj);
    AC = allocate_memory_bidimensional_matrix_double(mi, nj);
    BD = allocate_memory_bidimensional_matrix_double(nj, mi);

    // ************************************************************************************

    // THE 2D MESH IS CREATED
    for (ii = 0; ii <= mi - 1; ii++)
    {
        double dii1 = (double)ii;
        double dmi1 = (double)(mi - 1);
        xx[ii] = a * dii1 / dmi1;
    }

    for (jj = 0; jj <= nj - 1; jj++)
    {
        double djj1 = (double)jj;
        double dnj1 = (double)(nj - 1);
        yy[jj] = b * djj1 / dnj1;
    }

    // BIDIMENSIONAL ARRAYS INITIALISATION
    initialise_bidimensional_matrix_double(AI, mi, nj, 0.0);
    initialise_bidimensional_matrix_double(AC, mi, nj, 0.0);
    initialise_bidimensional_matrix_double(AD, mi, nj, 0.0);
    initialise_bidimensional_matrix_double(BI, nj, mi, 0.0);
    initialise_bidimensional_matrix_double(BD, nj, mi, 0.0);

    double value = (temp_fin + temp_ini) / 2.0;
    initialise_bidimensional_matrix_double(temper, mi, nj, value);
    fill_boundary_conditions_temper_matrix(temper, temp_ini, temp_fin, flux_aba, flux_arr);

    // ************************************************************************************

    // ASSEMBLING ARRAYS

    for (jj = 1; jj < nj - 1; jj++)
    {
        // IN X-DIRECTION
        ensambla_tdmax(AI, AC, AD, delta_x, delta_y, cond_ter, temp_ini, temp_fin, jj);
    }

    for (ii = 1; ii < mi - 1; ii++)
    {
        // IN Y-DIRECTION
        ensambla_tdmay(BI, BD, delta_x, delta_y, cond_ter, flux_aba, flux_arr, ii);
    }

    get_vector_independent_terms(BI, AI, AD, BD, vector_b);

    get_compressed_sparse_row_storages(BI, AI, AC, AD, BD, nonzero_elements, csr_values, csr_column_index, csr_row_pointer);

    // ************************************************************************************

// PARALLEL DATA REGION OPENS
#pragma acc data copyin(vector_b[ : columns_number],              \
                        csr_values[ : nonzero_elements],          \
                        csr_column_index[ : nonzero_elements],    \
                        csr_row_pointer[ : csr_row_pointer_size]) \
    copyout(results[ : columns_number])
    {
#pragma acc host_data use_device(vector_b, csr_values, csr_row_pointer, csr_column_index, results)
        {
            // cusolverSpDcsrlsvluHost(handle, rows_number, nonzero_elements, descrA, csr_values, csr_row_pointer, csr_column_index, vector_b, TOLERANCE, REORDER, results, &SINGULARITY);
            // cusolver_status = cusolverSpDcsrlsvqr(handle, rows_number, nonzero_elements, descrA, csr_values, csr_row_pointer, csr_column_index, vector_b, TOLERANCE, REORDER, results, &SINGULARITY);
            cusolver_status = cusolverSpDcsrlsvchol(handle, rows_number, nonzero_elements, descrA, csr_values, csr_row_pointer, csr_column_index, vector_b, TOLERANCE, REORDER, results, &SINGULARITY);
        }
    }
    // PARALLEL DATA REGION IS CLOSED

    // ************************************************************************************

    if (cusolver_status == CUSOLVER_STATUS_SUCCESS)
    {
        printf("SUCCESS!");
    }
    print_info_cusolver("Cholesky Solver", cusolver_status);
    // print_info_cusolver("QR Solver", cusolver_status);
    fill_matrix_temper_with_results(temper, results);
    // write_results(xx, yy, temper); // RESULTS WRITING

    // ************************************************************************************

    // RELEASE THE MEMORY ALLOCATED FOR THE MATRIX DESCRIPTOR AND CPU-SIDE RESOURCES USED BY THE CUSOLVERSP LIBRARY

    cusolver_status = cusolverSpDestroy(handle);
    // print_info_cusolver("Release resource for handle", cusolver_status);
    cusparse_status = cusparseDestroyMatDescr(descrA);
    // print_info_cusparse("Release resource for the matrix descriptor", cusparse_status);

    // FREE UP MEMORY

    free_memory_bidimensional_matrix_double(BD, nj, mi);
    free_memory_bidimensional_matrix_double(AD, mi, nj);
    free_memory_bidimensional_matrix_double(AC, mi, nj);
    free_memory_bidimensional_matrix_double(AI, mi, nj);
    free_memory_bidimensional_matrix_double(BI, nj, mi);
    free_memory_bidimensional_matrix_double(temper, mi, nj);

    free_memory_vector_double(results);
    free_memory_vector_double(vector_b);
    free_memory_vector_int(csr_row_pointer);
    free_memory_vector_int(csr_column_index);
    free_memory_vector_double(csr_values);
    free_memory_vector_double(yy);
    free_memory_vector_double(xx);

    return 0;
}
