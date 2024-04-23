#ifndef HERRAMIENTAS2_H
#define HERRAMIENTAS2_H

// Error handler
void error_handler(const char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

// Allocate memory to an vector of double data type
double *allocate_memory_vector(int vector_size)
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

// Allocate memory to an vector of int data type
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

// Allocate memory to an matrix of double data type
double **allocate_memory_matrix(int rows, int columns)
{
    double **matrix;
    size_t r_size_in_bytes;
    r_size_in_bytes = rows * sizeof(double *);

    // Allocate pointers to rows
    matrix = (double **)malloc(r_size_in_bytes);
    // Check that the memory allocation was successful.
    if (!matrix)
    {
        error_handler("Allocation failure in allocate_memory_matrix()");
    }

    // Allocate pointers to columns
    size_t c_size_in_bytes;
    c_size_in_bytes = columns * sizeof(double);

    for (int i = 0; i < rows; i++)
    {
        matrix[i] = (double *)malloc(c_size_in_bytes);
        // Check that the memory allocation was successful.
        if (!matrix[i])
        {
            // Liberar la memoria ya asignada para las filas antes de salir
            for (int j = 0; j < i; j++)
            {
                free(matrix[j]);
            }
            free(matrix);
            error_handler("Allocation failure in allocate_memory_matrix()");
        }
    }

    return matrix;
}

void free_matrix(double **matrix, int rows, int columns)
{
    // Liberar la memoria al final del programa
    for (int i = 0; i < rows; i++)
    {
        free(matrix[i]);
    }
    free(matrix);
}

void initialise_bidimensional_array_with_value(double **matrix, const int rows_num, const int columns_num, const double value)
{
    for (size_t i = 0; i < rows_num; i++)
    {

        for (size_t j = 0; j < columns_num; j++)
        {
            matrix[i][j] = value;
        }
    }
}

void initialise_array_with_value(double *vector, const int array_size, const double value)
{
    for (size_t i = 0; i < array_size; i++)
    {
        vector[i] = value;
    }
}

int get_nonzero_elements()
{
    int nonzero_elements_number = 3 * ni - 8;
    return nonzero_elements_number;
}

int get_column_index(int ii)
{
    int indice_columna;
    indice_columna = ii - 1;
    return indice_columna;
}

void get_csr_schema(
    double cond_ter,
    int nonzero_elements,
    double *csr_val,
    int *csr_col_ind,
    int *csr_row_ptr)
{
    int kk = 0;
    int jj = 0;
    // size_t stop = nonzero_elements + 1;

    for (size_t ii = 0; ii < nonzero_elements; ii += 3)
    {
        size_t index_left = ii - 1;
        size_t index_central = ii;
        size_t index_right = ii + 1;

        if (ii != 0)
        {
            csr_val[index_left] = -1.0 * cond_ter;
            csr_col_ind[index_left] = jj - 1;
        }
        
        csr_val[index_central] = 2.0 * cond_ter;
        csr_col_ind[index_central] = jj;

        if (ii != nonzero_elements - 1)
        {
            csr_val[index_right] = -1.0 * cond_ter;
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
    
}

void get_1D_heat_equation_matrix(double *matrix, int rows_num, int columns_num, double cond_ter)
{
    int matrix_size = rows_num * columns_num;

    initialise_array_with_value(matrix, matrix_size, 0.0);

    matrix[0] = 2.0 * cond_ter;
    matrix[1] = -1.0 * cond_ter;

    for (size_t i = 1; i < rows_num - 1; i++)
    {
        int index = i * rows_num + i;
        
        matrix[index - 1] = -1.0 * cond_ter;
        matrix[index] = 2.0 * cond_ter;
        matrix[index + 1] = -1.0 * cond_ter;

    }

    matrix[matrix_size - 2] = -1.0 * cond_ter;
    matrix[matrix_size - 1] = 2.0 * cond_ter;
    
}

// b: vector de terminos independientes b
void get_vector_independent_terms(double *b, int vector_size, double cond_term, double start_temp, double end_temp)
{
    b[0] = start_temp * cond_term;

    for (size_t i = 1; i < vector_size-1; i++)
    {
        b[i] = 0.0;
    }

    b[vector_size - 1] = end_temp * cond_term;
    
}

void get_temper_array(double *temper, double *outcome_vector, double start_temp, double end_temp)
{
    temper[0] = start_temp;
    for (size_t i = 1; i < ni-1; i++)
    {
        temper[i] = outcome_vector[i - 1];
    }
    temper[ni - 1] = end_temp;

}

// Print to console a Matrix
void print_matrix(double **matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("%f ", matrix[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

// Print to console a Vector
void print_vector(double *vector, int vector_size)
{
    for (int i = 0; i < vector_size; ++i)
    {
        printf("%f ", vector[i]);
        printf("\n");
    }
    printf("\n");
}

// Print to console a Vector
void print_vector_int(int *vector, int vector_size)
{
    for (int i = 0; i < vector_size; ++i)
    {
        printf("%d ", vector[i]);
        printf("\n");
    }
    printf("\n");
}

void print_flat_matrix(double *matrix, int rows_num, int columns_num)
{
    for (size_t i = 0; i < rows_num; i++)
    {
        for (size_t j = 0; j < columns_num; j++)
        {
            printf("%.7f ", matrix[i * columns_num + j]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_formato_csr(double *csrVal, int *csrIndCol, int *csrPtr, int elementos_no_cero, int tamanio_ptr)
{
    for (int ii = 0; ii < elementos_no_cero; ii++)
    {
        if (ii < tamanio_ptr)
        {
            printf("Val[%d] = %f | ColInd[%d] = %d | Ptr[%d] = %d\n", ii, csrVal[ii], ii, csrIndCol[ii], ii, csrPtr[ii]);
            continue;
        }
        printf("Val[%d] = %f | ColInd[%d] = %d\n", ii, csrVal[ii], ii, csrIndCol[ii]);
    }
}

#endif