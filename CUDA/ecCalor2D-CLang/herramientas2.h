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

void inicializar_matriz(double **matriz, const int filas, const int columnas, const double valor_inicial)
{
    for (int i = 0; i < filas; i++)
    {

        for (int j = 0; j < columnas; j++)
        {
            matriz[i][j] = valor_inicial;
        }
    }
}

void inicializar_vector(double *vector, const int tamanio_vector, const double valor_inicial)
{
    for (int i = 0; i < tamanio_vector; i++)
    {
        vector[i] = valor_inicial;
    }
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