#ifndef HERRAMIENTAS_H
#define HERRAMIENTAS_H

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
            for (int j = 0; j < i; j++) {
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
    for (int i = 0; i < rows; i++) {
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

#endif