#ifndef HERRAMIENTAS_H
#define HERRAMIENTAS_H

// ************************************************************************************

void clang_error_handler(const char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

char *cusolver_get_error_status(cusolverStatus_t status)
{
    switch (status)
    {
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
    return "UNRECOGNIZED ERROR CODE";
}

void print_info_cusolver(char *event_to_be_reported, cusolverStatus_t status)
{
    char *info;

    if (status != CUSOLVER_STATUS_SUCCESS)
    {
        info = cusolver_get_error_status(status);
        printf("cusolver error - %s: ", event_to_be_reported);
        printf("%s...\n\n", info);
    }
    else
    {
        printf("%s: SUCCESS ", event_to_be_reported);
        printf("\n\n");
    }
}

void print_info_cusparse(char *event_to_be_reported, cusparseStatus_t status)
{
    const char *info;

    if (status != CUSPARSE_STATUS_SUCCESS)
    {
        info = cusparseGetErrorName(status);
        printf("cusparse error - %s: ", event_to_be_reported);
        printf("%s...\n\n", info);
    }
    else
    {
        printf("%s: SUCCESS", event_to_be_reported);
        printf("\n\n");
    }
}

// ************************************************************************************

double *allocate_memory_vector_double(const int vector_size)
{
    double *vector;
    size_t size_in_bytes;

    size_in_bytes = vector_size * sizeof(double);
    vector = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        clang_error_handler("Allocation failure in allocate_memory_vector_double()");
    }
    return vector;
}

int *allocate_memory_vector_int(const int vector_size)
{
    int *vector;
    size_t size_in_bytes;

    size_in_bytes = vector_size * sizeof(int);
    vector = (int *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        clang_error_handler("Allocation failure in allocate_memory_vector_int()");
    }
    return vector;
}

// ************************************************************************************

int *allocate_memory_unidimensional_matrix_double(const int rows_number, const int columns_number)
{
    int *matrix;
    size_t size_in_bytes;

    size_in_bytes = rows_number * columns_number * sizeof(int);
    matrix = (int *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!matrix)
    {
        clang_error_handler("Allocation failure in allocate_memory_vector_int()");
    }
    return matrix;
}

double **allocate_memory_bidimensional_matrix_double(const int rows_number, const int columns_number)
{
    double **matrix;
    size_t row_size_in_bytes;

    row_size_in_bytes = rows_number * sizeof(double *);
    // Allocate pointers to rows_number
    matrix = (double **)malloc(row_size_in_bytes);

    // Check that the memory allocation was successful.
    if (!matrix)
    {
        clang_error_handler("Allocation failure in allocate_bidimensional_memory_matrix_double() - rows");
    }

    // Allocate pointers to columns_number
    size_t column_size_in_bytes;
    column_size_in_bytes = columns_number * sizeof(double);

    for (size_t ii = 0; ii < rows_number; ii++)
    {
        matrix[ii] = (double *)malloc(column_size_in_bytes);
        // Check that the memory allocation was successful.
        if (!matrix[ii])
        {
            // Liberar la memoria ya asignada para las rows_number antes de salir
            for (size_t jj = 0; jj < ii; jj++)
            {
                free(matrix[jj]);
            }

            free(matrix);
            clang_error_handler("Allocation failure in allocate_bidimensional_memory_matrix_double() - columns");
        }
    }

    return matrix;
}

// ************************************************************************************

void free_memory_vector_double(double *vector)
{
    free(vector);
}

void free_memory_vector_int(int *vector)
{
    free(vector);
}

void free_memory_unidimensional_matrix_double(double *matrix)
{
    free(matrix);
}

void free_memory_bidimensional_matrix_double(double **matrix, const int rows_number, const int columns_number)
{
    // Freeing up memory at the end of the programme
    for (size_t ii = 0; ii < rows_number; ii++)
    {
        free(matrix[ii]);
    }
    free(matrix);
}

// ************************************************************************************

void initialise_vector_double(double *vector, const int vector_size, const double initial_value)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        vector[ii] = initial_value;
    }
}

void initialise_unidimensional_matrix_double(double *matrix, const int rows_number, const int columns_number, const double initial_value)
{
    int matrix_size = rows_number * columns_number;
    for (size_t ii = 0; ii < matrix_size; ii++)
    {
        matrix[ii] = initial_value;
    }
}

void initialise_bidimensional_matrix_double(double **matriz, const int rows_number, const int columns_number, const double initial_value)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {

        for (size_t jj = 0; jj < columns_number; jj++)
        {
            matriz[ii][jj] = initial_value;
        }
    }
}

// ************************************************************************************
int get_nonzero_elements()
{

    int nonzero_elements = 12 + (mi + nj - 8) * 8 + (mi - 4) * (nj - 4) * 5;
    return nonzero_elements;
}

int get_column_index(const int ii, const int jj)
{
    int column_index;
    // column_index = (jj - 1) * (mi - 2) + ii - (mi - 1) - 1;
    // column_index = (jj - 1) * (mi - 2) + ii - (mi - 1); // for index-one
    column_index = jj * (mi - 2) + ii - (mi - 1); // for index-zero
    return column_index;
}

void get_compressed_sparse_row_storages(double **BI, double **AI, double **AC, double **AD, double **BD, int nonzero_elements, double *csr_values, int *csr_column_index, int *csr_row_pointer)
{
    int kk = -1;      // csr_values index
    int tt = 0;       // csr_row_pointer index
    int counter = -1; // nonzero elements counter
    csr_row_pointer[tt] = 0;

    for (size_t jj = 1; jj < nj - 1; jj++)
    {
        for (size_t ii = 1; ii < mi - 1; ii++)
        {
            if (jj >= 2 && jj < nj - 1)
            {
                kk++;
                csr_values[kk] = BI[jj - 1][ii];
                csr_column_index[kk] = get_column_index(ii, jj - 1);
                counter++;
            }

            if (ii >= 2 && ii < mi - 1)
            {
                kk++;
                csr_values[kk] = AI[ii - 1][jj];
                csr_column_index[kk] = get_column_index(ii - 1, jj);
                counter++;
            }

            kk++;

            if (kk >= nonzero_elements)
                break;

            csr_values[kk] = AC[ii][jj];
            csr_column_index[kk] = get_column_index(ii, jj);
            counter++;

            if (ii >= 1 && ii < mi - 2)
            {
                kk++;
                csr_values[kk] = AD[ii + 1][jj];
                csr_column_index[kk] = get_column_index(ii + 1, jj);
                counter++;
            }

            if (jj >= 1 && jj < nj - 2)
            {
                kk++;
                csr_values[kk] = BD[jj + 1][ii];
                csr_column_index[kk] = get_column_index(ii, jj + 1);
                counter++;
            }

            tt++;
            csr_row_pointer[tt] = counter + 1;
        }
    }
}

void get_vector_independent_terms(double **BI, double **AI, double **AD, double **BD, double *b)
{
    int kk = 0;
    for (size_t jj = 1; jj < nj - 1; jj++)
    {
        for (size_t ii = 1; ii < mi - 1; ii++)
        {
            if (jj == 1)
            {
                b[kk] += -BI[jj - 1][ii];
            }

            if (ii == 1)
            {
                b[kk] += -AI[ii - 1][jj];
            }

            if (ii == mi - 2)
            {
                b[kk] += -AD[ii + 1][jj];
            }

            if (jj == nj - 2)
            {
                b[kk] += -BD[jj + 1][ii];
            }

            kk++;
        }
    }
}

// ************************************************************************************

void fill_boundary_conditions_temper_matrix(double **matrix, double temp_ini, double temp_fin, double flux_aba, double flux_arr)
{
    for (size_t jj = 1; jj < nj - 1; jj++)
    {
        // BI
        matrix[0][jj] = temp_ini;
        matrix[mi - 1][jj] = temp_fin;
    }
    for (size_t ii = 1; ii < mi - 1; ii++)
    {
        matrix[ii][0] = flux_aba;
        matrix[ii][nj - 1] = flux_arr;
    }
}

void fill_matrix_temper_with_results(double **matriz, double *vector)
{
    size_t kk = 0; // vector index

    for (size_t jj = 1; jj < nj - 1; jj++)
    {
        for (size_t ii = 1; ii < mi - 1; ii++)
        {
            matriz[ii][jj] = vector[kk];
            kk++;
        }
    }
}

// ************************************************************************************

void print_vector_double(double *vector, const int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ++ii)
    {
        printf("%4.4f ", vector[ii]);
        printf("\n");
    }
    printf("\n");
}

void print_vector_int(int *vector, const int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ++ii)
    {
        printf("%d ", vector[ii]);
        printf("\n");
    }
    printf("\n");
}

void print_unidimensional_matrix_double(double *matrix, const int rows_number, const int columns_number)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {
        for (size_t jj = 0; jj < columns_number; jj++)
        {
            size_t index = ii * rows_number + jj;
            printf("%4.4f ", matrix[index]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_bidimensional_matrix_double(double **matrix, const int rows_number, const int columns_number)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {
        for (size_t jj = 0; jj < columns_number; jj++)
        {
            printf("%4.4f ", matrix[ii][jj]);
        }
        printf("\n");
    }
    printf("\n");
}

void print_csr_format_stores(double *csr_values, int *csr_column_index, int *csr_row_pointer, int nonzero_elements, const int csr_row_pointer_size)
{
    for (int ii = 0; ii < nonzero_elements; ii++)
    {
        if (ii < csr_row_pointer_size)
        {
            printf("Val[%d] = %f | ColInd[%d] = %d | Ptr[%d] = %d\n", ii, csr_values[ii], ii, csr_column_index[ii], ii, csr_row_pointer[ii]);
            continue;
        }
        printf("Val[%d] = %f | ColInd[%d] = %d\n", ii, csr_values[ii], ii, csr_column_index[ii]);
    }
}

// ************************************************************************************

void write_results(double *x_coordinates, double *y_coordinates, double **temper)
{
    // Abrir el archivo para escribir
    FILE *file = fopen("cusolver.101", "w");

    if (file != NULL)
    {

        // Utilizar un bucle para imprimir y guardar los datos
        for (size_t ii = 0; ii < mi; ii++)
        {
            for (size_t jj = 0; jj < nj; jj++)
            {
                fprintf(file, "%f %f %f\n", x_coordinates[ii], y_coordinates[jj], temper[ii][jj]);
            }

            fprintf(file, "\n");
        }

        // Cerrar el archivo después de escribir
        fclose(file);
    }
    else
    {

        printf("Error al abrir el archivo.\n");
    }
}

#endif