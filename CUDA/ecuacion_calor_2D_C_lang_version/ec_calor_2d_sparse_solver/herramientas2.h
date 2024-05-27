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

int obtener_total_elementos_no_cero()
{
    int elementos_no_cero = 12 + (mi + nj - 8) * 8 + (mi - 4) * (nj - 4) * 5;
    return elementos_no_cero;
}

int obtener_indice_columna(int ii, int jj)
{
    int indice_columna;
    indice_columna = (jj - 1) * (mi - 2) + ii - (mi - 1) - 1;
    return indice_columna;
}

void obtener_formato_csr(
    double **BI,
    double **AI,
    double **AC,
    double **AD,
    double **BD,
    int elementos_no_cero,
    double *csrVal,
    int *csrIndCol,
    int *csrPtr)
{
    int kk = -1;
    int tt = 0;
    int counter = -1;
    csrPtr[tt] = 0;

    for (int jj = 1; jj < nj - 1; jj++)
    {
        for (int ii = 1; ii < mi - 1; ii++)
        {
            if (jj >= 2 && jj < nj - 1)
            {
                kk++;
                csrVal[kk] = BI[jj - 1][ii];
                csrIndCol[kk] = obtener_indice_columna(ii + 1, jj);
                counter++;
            }

            if (ii >= 2 && ii < mi - 1)
            {
                kk++;
                csrVal[kk] = AI[ii - 1][jj];
                csrIndCol[kk] = obtener_indice_columna(ii, jj + 1);
                counter++;
            }

            kk++;

            if (kk >= elementos_no_cero)
                break;
            
            csrVal[kk] = AC[ii][jj];
            csrIndCol[kk] = obtener_indice_columna(ii + 1, jj + 1);
            counter++;

            if (ii >= 1 && ii < mi - 2)
            {
                kk++;
                csrVal[kk] = AD[ii + 1][jj];
                csrIndCol[kk] = obtener_indice_columna(ii + 2, jj + 1);
                counter++;
            }

            if (jj >= 1 && jj < nj - 2)
            {
                kk++;
                csrVal[kk] = BD[jj + 1][ii];
                csrIndCol[kk] = obtener_indice_columna(ii + 1, jj + 2);
                counter++;
            }

            tt++;
            csrPtr[tt] = counter + 1;
        }

    }
    
}

// b: vector de terminos independientes b
void obtener_vector_terminos_independientes(
    double **BI,
    double **AI,
    double **AD,
    double **BD,
    double *b)
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

void print_flat_matrix(double *matrix, int matrix_size)
{
    for (int i = 0; i < matrix_size; i++)
    {
        for (int j = 0; j < matrix_size; j++)
        {
            printf("%.7f ", matrix[i * matrix_size + j]);
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

void llenar_matriz_temper(double **matriz)
{
    for (int j = 1; j < nj-1; j++)
    {
        matriz[0][j] = 0.0;
        matriz[mi-1][j] = 0.0;
    }
    for (int i = 1; i < mi-1; i++)
    {
        matriz[i][0] = 0.0;
        matriz[i][nj-1] = 3.0;
    }
}

void completar_matriz_temper(double **matriz, double *vector, int tam_matriz_completa)
{
    for (int j = 1; j < mi-1; j++)
    {
        for (int i = 1; i < nj-1; i++)
        {
            matriz[j][i] = vector[(i) + 3 * (i-1) - 1];
        }
    }
    
}

void csr_to_dense(
    double *csr_val,
    int *csr_col_ind,
    int *csr_row_ptr,
    double *dense_matrix,
    int matrix_size
)
{
    // Iterate through each row in the CSR format
  for (int i = 0; i < matrix_size; i++) {
    // Loop through the non-zero elements in the current row
    for (int j = csr_row_ptr[i]; j < csr_row_ptr[i + 1]; j++) {
      int col_index = csr_col_ind[j];
      double value = csr_val[j];
      dense_matrix[i * matrix_size + col_index] = value; // Set the corresponding element in the dense matrix
    }
  }
}

#endif