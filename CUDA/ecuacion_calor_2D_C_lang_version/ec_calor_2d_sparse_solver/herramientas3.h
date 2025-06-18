#ifndef HERRAMIENTAS3_H
#define HERRAMIENTAS3_H

// **************************************************************************************************
// **************************************************************************************************

void error_handler(const char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

// **************************************************************************************************
// **************************************************************************************************

double *allocate_memory_vector_double(const int vector_size)
{
    double *vector;
    size_t size_in_bytes;

    size_in_bytes = vector_size * sizeof(double);
    vector = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        error_handler("Allocation failure in allocate_memory_vector_double()");
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
        error_handler("Allocation failure in allocate_memory_vector_int()");
    }
    return vector;
}

double *allocate_memory_unidimensional_matrix_double(const int rows_number, const int columns_number)
{
    double *matrix;
    size_t size_in_bytes;

    size_in_bytes = rows_number * columns_number * sizeof(double);
    matrix = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!matrix)
    {
        error_handler("Allocation failure in allocate_memory_matrix_unidimensional_double()");
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
        error_handler("Allocation failure in allocate_bidimensional_memory_matrix_double() - rows");
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
            // Liberar la memoria ya asignada para las filas antes de salir
            for (size_t jj = 0; jj < ii; jj++)
            {
                free(matrix[jj]);
            }

            free(matrix);
            error_handler("Allocation failure in allocate_bidimensional_memory_matrix_double() - columns");
        }
    }

    return matrix;
}

// **************************************************************************************************
// **************************************************************************************************

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

// **************************************************************************************************
// **************************************************************************************************

void initialise_bidimensional_matrix_double(double **matrix, const int rows_number, const int columns_number, const double initial_value)
{
    for (size_t ii = 0; ii < rows_number; ii++)
    {

        for (size_t jj = 0; jj < columns_number; jj++)
        {
            matrix[ii][jj] = initial_value;
        }
    }
}

void initialise_vector_double(double *vector, const int vector_size, const double initial_value)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        vector[ii] = initial_value;
    }
}

// **************************************************************************************************
// **************************************************************************************************

int get_nonzero_elements()
{
    int nonzero_elements = 12 + (mi + nj - 8) * 8 + (mi - 4) * (nj - 4) * 5;
    return nonzero_elements;
}

int get_column_index(const int ii, const int jj)
{
    int column_index;
    column_index = (jj - 1) * (mi - 2) + ii - (mi - 1) - 1;
    // column_index = (jj - 1) * (mi - 2) + ii - (mi - 1); // for index-one
    // column_index = jj * (mi - 2) + ii - (mi - 1); // for index-zero
    return column_index;
}

int gcindex(const int ii, const int jj)
{
    int column_index;
    // column_index = (jj - 2) * (mi - 2) + ii - 1; for index_one
    column_index = (jj - 1) * (mi - 2) + ii - 1; // for index_zero
    return column_index;
}


void get_csr_format_stores(double **BI, double **AI, double **AC, double **AD, double **BD, const int nonzero_elements, double *csr_values, int *csr_column_index, int *csr_row_pointer)
{
    int kk = -1;
    int tt = 0;
    int counter = -1;
    csr_row_pointer[tt] = 0;

    for (int jj = 1; jj < nj - 1; jj++)
    {
        for (int ii = 1; ii < mi - 1; ii++)
        {
            if (jj >= 2 && jj < nj - 1)
            {
                kk++;
                csr_values[kk] = BI[jj - 1][ii];
                csr_column_index[kk] = get_column_index(ii + 1, jj);
                counter++;
            }

            if (ii >= 2 && ii < mi - 1)
            {
                kk++;
                csr_values[kk] = AI[ii - 1][jj];
                csr_column_index[kk] = get_column_index(ii, jj + 1);
                counter++;
            }

            kk++;

            if (kk >= nonzero_elements)
                break;
            
            csr_values[kk] = AC[ii][jj];
            csr_column_index[kk] = get_column_index(ii + 1, jj + 1);
            counter++;

            if (ii >= 1 && ii < mi - 2)
            {
                kk++;
                csr_values[kk] = AD[ii + 1][jj];
                csr_column_index[kk] = get_column_index(ii + 2, jj + 1);
                counter++;
            }

            if (jj >= 1 && jj < nj - 2)
            {
                kk++;
                csr_values[kk] = BD[jj + 1][ii];
                csr_column_index[kk] = get_column_index(ii + 1, jj + 2);
                counter++;
            }

            tt++;
            csr_row_pointer[tt] = counter + 1;
        }

    }
    
}

void gcsrfs(double *values, int *row_pointer, int *column_index, int nonzero_elements, double cond_ter, double deltax, double deltay)
{    
    int kk = -1;
    int tt = 0;
    int counter = -1;
    row_pointer[tt] = 0;

    for (size_t jj = 1; jj < nj - 1; jj++)
    {
        for (size_t ii = 1; ii < mi - 1; ii++)
        {
            // Aparecen las BI
            if (jj >= 2 && jj < nj - 1)
            {
                kk++;
                values[kk] = -1.0 * cond_ter / (deltay * deltay);
                column_index[kk] = gcindex(ii, jj - 1);
                counter++;
            }

            // AI
            if (ii >= 2 && ii < mi - 1)
            {
                kk++;
                values[kk] = -1 * cond_ter / (deltax * deltax);
                column_index[kk] = gcindex(ii - 1, jj);
                counter++;
            }

            kk++;

            if (kk >= nonzero_elements)
                break;
            
            // AC
            values[kk] = 2 * cond_ter * (1 / (deltax * deltax) + 1.0 / (deltay * deltay));
            column_index[kk] = gcindex(ii, jj);
            counter++;

            // AD
            if (ii >= 1 && ii < mi - 2)
            {
                kk++;
                values[kk] = -1.0 * cond_ter / (deltax * deltax);
                column_index[kk] = gcindex(ii + 1, jj);
                counter++;
            }

            // BD
            if (jj >= 1 && jj < nj - 2)
            {
                kk++;
                values[kk] = -1.0 * cond_ter / (deltay * deltay);
                column_index[kk] = gcindex(ii, jj + 1);
                counter++;
            }

            tt++;
            row_pointer[tt] = counter + 1;
        }

    }

}

void get_vector_independent_terms(double **BI, double **AI, double **AD, double **BD, double *b)
{
    int kk = 0;
    for (int jj = 1; jj < nj-1; jj++)
    {
        for (int ii = 1; ii < mi-1; ii++)
        {
            if (jj == 1)
            {
                b[kk] += -BI[jj-1][ii];
            }

            if (ii == 1)
            {
                b[kk] += -AI[ii-1][jj];
            }

            if (ii == mi-2)
            {
                b[kk] += -AD[ii+1][jj];
            }

            if (jj == nj-2)
            {
                b[kk] += -BD[jj+1][ii];
            }
            
            kk++;
        }
        
    }
    
}

void gvit(double *b, double cond_ter, double deltax, double deltay)
{
    int kk = 0;
    for (int jj = 1; jj < nj-1; jj++)
    {
        for (int ii = 1; ii < mi-1; ii++)
        {
            // BI
            if (jj == 1)
            {
                b[kk] += 0.0 * cond_ter/ (deltay * deltay);
            }

            // AI
            if (ii == 1)
            {
                b[kk] += 0.0* cond_ter / (deltax * deltax);
            }

            // AD
            if (ii == mi-2)
            {
                b[kk] += 0.0* cond_ter / (deltax * deltax);
            }

            // BD
            if (jj == nj-2)
            {
                b[kk] += 3.0* cond_ter / (deltay * deltay);
            }
            
            kk++;
        }
        
    }

}

// **************************************************************************************************
// **************************************************************************************************

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

// **************************************************************************************************
// **************************************************************************************************

void fill_boundary_conditions_temper_matrix(double **matrix)
{
    for (size_t jj = 1; jj < nj-1; jj++)
    {
        // BI
        matrix[0][jj] = 0.0;
        matrix[mi-1][jj] = 0.0;
    }
    for (size_t ii = 1; ii < mi-1; ii++)
    {
        matrix[ii][0] = 0.0;
        matrix[ii][nj-1] = 3.0;
    }
}

void completar_matriz_temper(double **matrix, double *vector, int tam_matriz_completa)
{
    size_t vector_index = 0;
    
    for (size_t jj = 1; jj < nj-1; jj++)
    {
        for (size_t ii = 1; ii < mi-1; ii++)
        {
            matrix[ii][jj] = vector[vector_index];
            vector_index++;
            // matrix[ii][jj] = vector[(ii) + 3 * (ii-1) - 1];
        }
    }
    
}

void csr_to_dense(double *csr_val, int *csr_column_index, int *csr_row_pointer, double *dense_matrix, int matrix_size)
{
    // Iterate through each current_row in the CSR format
  for (size_t ii = 0; ii < matrix_size; ii++) {
    // Loop through the non-zero elements in the current current_row
    for (size_t jj = csr_row_pointer[ii]; jj < csr_row_pointer[ii + 1]; jj++) {
      int col_index = csr_column_index[jj];
      double value = csr_val[jj];
      dense_matrix[ii * matrix_size + col_index] = value; // Set the corresponding element in the dense matrix
    }
  }
}

void comparator(int *vector_1, int *vector_2, int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        if (vector_1[ii] != vector_2[ii])
        {
            printf("Los vectores son distintos en la posicion %zu...\n", ii);
            printf("Programa abortado...\n");
            break;
        }
        
    }
    printf("Los vectores son iguales...\n");
    
}

void comparator_double(double *vector_1, double *vector_2, int vector_size)
{
    for (size_t ii = 0; ii < vector_size; ii++)
    {
        if (vector_1[ii] != vector_2[ii])
        {
            printf("Los vectores son distintos en la posicion %zu...\n", ii);
            printf("Programa abortado...\n");
            break;
        }
        
    }
    printf("Los vectores son iguales...\n");
    
}

#endif