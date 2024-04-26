#ifndef HERRAMIENTAS_H
#define HERRAMIENTAS_H
#include <cusolverSp.h>

// Error handler
void error_handler(const char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

double *allocate_memory_vector_double(int vector_size)
{
    double *vector;
    size_t size_in_bytes;
    size_in_bytes = vector_size * sizeof(double);

    vector = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        error_handler("Allocation failure in allocate_memory_vector()");
    }
    return vector;
}

int *allocate_memory_vector_int(int vector_size)
{
    int *vector;
    size_t size_in_bytes;
    size_in_bytes = vector_size * sizeof(int);

    vector = (int *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        error_handler("Allocation failure in allocate_memory_vector()");
    }
    return vector;
}

double *allocate_memory_matrix_flat_double(int rows_number, int columns_number)
{
    double *matrix;
    size_t size_in_bytes;
    size_in_bytes = rows_number * columns_number * sizeof(double);

    matrix = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!matrix)
    {
        error_handler("Allocation failure in allocate_memory_vector()");
    }
    return matrix;
}

// Allocate memory to an matrix of double data type
double **allocate_memory_matrix_double(int rows_number, int columns_number)
{
    double **matrix;
    size_t r_size_in_bytes;
    r_size_in_bytes = rows_number * sizeof(double *);

    // Allocate pointers to rows_number
    matrix = (double **)malloc(r_size_in_bytes);
    // Check that the memory allocation was successful.
    if (!matrix)
    {
        error_handler("Allocation failure in allocate_memory_matrix()");
    }

    // Allocate pointers to columns_number
    size_t c_size_in_bytes;
    c_size_in_bytes = columns_number * sizeof(double);

    for (int ii = 0; ii < rows_number; ii++)
    {
        matrix[ii] = (double *)malloc(c_size_in_bytes);
        // Check that the memory allocation was successful.
        if (!matrix[ii])
        {
            // Liberar la memoria ya asignada para las filas antes de salir
            for (int jj = 0; jj < ii; jj++)
            {
                free(matrix[jj]);
            }
            free(matrix);
            error_handler("Allocation failure in allocate_memory_matrix()");
        }
    }

    return matrix;
}

void free_memory_vector_double(double *vector)
{
    free(vector);
}

void free_memory_vector_int(int *vector)
{
    free(vector);
}


void free_memory_matrix_flat_double(double *matrix)
{
    free(matrix);
}

void free_matrix_double(double **matrix, int rows_number, int columns_number)
{
    // Liberar la memoria al final del programa
    for (int ii = 0; ii < rows_number; ii++)
    {
        free(matrix[ii]);
    }
    free(matrix);
}

void initialise_bidimensional_array_with_value(double **matrix, const int rows_number, const int columns_number, const double value)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {

        for (size_t jj = 0; jj < columns_number; jj++)
        {
            matrix[ii][jj] = value;
        }
    }
}

void initialise_array_with_value(double *vector, const int vector_size, const double value)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        vector[ii] = value;
    }
}

int get_nonzero_elements()
{
    int nonzero_elements_number = 3 * ni - 8;
    return nonzero_elements_number;
}

void get_csr_schema(
    double cond_ter,
    int nonzero_elements,
    double *csr_values,
    int *csr_col_ind,
    int *csr_row_ptr)
{
    int kk = 0;
    int jj = 0;

    for (size_t ii = 0; ii < nonzero_elements; ii += 3)
    {
        size_t index_left = ii - 1;
        size_t index_central = ii;
        size_t index_right = ii + 1;

        if (ii != 0)
        {
            csr_values[index_left] = -1.0 * cond_ter;
            csr_col_ind[index_left] = jj - 1;
        }
        
        csr_values[index_central] = 2.0 * cond_ter;
        csr_col_ind[index_central] = jj;

        if (ii != nonzero_elements - 1)
        {
            csr_values[index_right] = -1.0 * cond_ter;
            csr_col_ind[index_right] = jj + 1;
        }
        
        jj++;

        if (ii != 0)
        {
            csr_row_ptr[kk] = index_left;
        } else {
            csr_row_ptr[kk] = index_central;
        }

        kk++;
    }

    csr_row_ptr[kk] = nonzero_elements;
    
}

void get_1D_heat_equation_matrix_dense(double *matrix, int rows_number, int columns_number, double cond_ter)
{
    int matrix_size = rows_number * columns_number;

    initialise_array_with_value(matrix, matrix_size, 0.0);

    matrix[0] = 2.0 * cond_ter;
    matrix[1] = -1.0 * cond_ter;

    for (size_t ii = 1; ii < rows_number - 1; ii++)
    {
        int index = ii * rows_number + ii;
        
        matrix[index - 1] = -1.0 * cond_ter;
        matrix[index] = 2.0 * cond_ter;
        matrix[index + 1] = -1.0 * cond_ter;

    }

    matrix[matrix_size - 2] = -1.0 * cond_ter;
    matrix[matrix_size - 1] = 2.0 * cond_ter;
    
}

void get_vector_independent_terms(double *vector, int vector_size, double cond_term, double start_temp, double end_temp)
{
    vector[0] = start_temp * cond_term;

    for (size_t ii = 1; ii < vector_size - 1; ii++)
    {
        vector[ii] = 0.0;
    }

    vector[vector_size - 1] = end_temp * cond_term;
    
}

void get_temper_vector(double *temper, double *results_vector, double start_temp, double end_temp)
{
    temper[0] = start_temp;
    for (size_t ii = 1; ii < ni - 1; ii++)
    {
        temper[ii] = results_vector[ii - 1];
    }
    temper[ni - 1] = end_temp;

}

void print_vector_int(int *vector, int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        printf("%d ", vector[ii]);
        printf("\n");
    }
    printf("\n");
}

void print_vector_double(double *vector, int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        printf("%f.4 ", vector[ii]);
        printf("\n");
    }
    printf("\n");
}

void print_matrix_flat_double(double *matrix, int rows_number, int columns_number)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {
        for (size_t jj = 0; jj < columns_number; jj++)
        {
            printf("%.4f ", matrix[ii * rows_number + jj]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_matrix_double(double **matrix, int rows_number, int columns_number)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {
        for (size_t jj = 0; jj < columns_number; jj++)
        {
            printf("%.4f ", matrix[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_formato_csr(double *csr_values, int *csr_col_ind, int *csr_row_ptr, int nonzero_elements, int csr_row_ptr_size)
{
    for (size_t ii = 0; ii < nonzero_elements; ii++)
    {
        if (ii < csr_row_ptr_size)
        {
            printf("Val[%lu] = %f | ColInd[%lu] = %d | Ptr[%lu] = %d\n", ii, csr_values[ii], ii, csr_col_ind[ii], ii, csr_row_ptr[ii]);
            continue;
        }
        printf("Val[%lu] = %f | ColInd[%lu] = %d\n", ii, csr_values[ii], ii, csr_col_ind[ii]);
    }
}

void print_info_lapacke_solver(int info)
{
    printf("\nEstado del solver de LAPACKE:\n");
    
    if (info == 0)
    {
        printf("LAPACKE SOLVER STATUS SUCCESS!...\n");
    } 
    
    if (info < 0)
    {
        printf("THE %d-TH ARGUMENT HAD AN ILLEGAL VALUE.\n", info);
    } 
    
    if (info > 0)
    {
        printf("U(%d, %d) IS EXACTLY ZERO.\n", info, info);
        printf("THE FACTORIZATION HAS BEEN COMPLETED, BUT THE FACTOR U IS EXACTLY\n");
        printf("SINGULAR, SO THE SOLUTION COULD NOT BE COMPUTED.\n");
    }
}

char *cuda_solver_get_status(cusolverStatus_t status)
{

    switch (status) {
        case CUSOLVER_STATUS_SUCCESS:
            return "CUSOLVER_STATUS_SUCCESS";
        case CUSOLVER_STATUS_NOT_INITIALIZED:
            return "CUSOLVER_STATUS_NOT_INITIALIZED";
        case CUSOLVER_STATUS_ALLOC_FAILED:
            return "CUSOLVER_STATUS_ALLOC_FAILED";
        case CUSOLVER_STATUS_INVALID_VALUE:
            return "CUSOLVER_STATUS_INVALID_VALUE";
        case CUSOLVER_STATUS_ARCH_MISMATCH:
            return "CUSOLVER_STATUS_ARCH_MISMATCH";
        case CUSOLVER_STATUS_EXECUTION_FAILED:
            return "CUSOLVER_STATUS_EXECUTION_FAILED";
        case CUSOLVER_STATUS_INTERNAL_ERROR:
            return "CUSOLVER_STATUS_INTERNAL_ERROR";
        case CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_MATRIX_TYPE_NOT_SUPPORTED";
        case CUSOLVER_STATUS_NOT_SUPPORTED:
            return "CUSOLVER_STATUS_NOT_SUPPORTED";
  }
  return "NULL";
}

void print_info_cuda_solver(cusolverStatus_t status)
{
    char *info;

    info = cuda_solver_get_status(status);

    printf("\nEstado del solver de CUDA:\n");
    printf("%s...\n", info);

}
#endif