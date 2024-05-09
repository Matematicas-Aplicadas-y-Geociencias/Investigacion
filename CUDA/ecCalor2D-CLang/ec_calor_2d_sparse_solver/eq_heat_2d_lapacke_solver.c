#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#include "constantes3.h"
#include "herramientas3.h"

int main(int argc, char const *argv[])
{
    // *************************************************************
    // Declaracion de variables
    int ii, jj; //, kk;
    
    // Variables de tamano de la placa
    double a, b, deltax, deltay;

    // Parametros fisicos del problema
    double cond_ter, temp_ini, temp_fin;
    double flux_aba, flux_arr; //, alpha;

    // Arreglos para incognitas, termino fuente y posicion
    double **temper, **temp_ant;
    double *xx;
    double *yy;

    //Matriz a invertir
    double **AI, **AD, **AC;
    double **BI, **BD, **BC;

    // CSR Sparse Format
    double *csr_values;
    int *csr_column_index;
    int *csr_row_pointer;
    double *bx;
    double *results;

    int nonzero_elements;
    int rows_number;
    int columns_number;
    int csr_row_pointer_size;

    double TOL = 1e-5;
    int REORDER = 3;
    int singularity = 0;

    // *************************************************************
    double *values;
    int *column_index;
    int *row_pointer;
    double *brh;

    // *************************************************************
    
    temper = allocate_memory_bidimensional_matrix_double(mi, nj);
    temp_ant = allocate_memory_bidimensional_matrix_double(mi, nj);

    xx = allocate_memory_vector_double(mi);
    yy = allocate_memory_vector_double(nj);

    
    AI = allocate_memory_bidimensional_matrix_double(mi, nj);
    AD = allocate_memory_bidimensional_matrix_double(mi, nj);
    AC = allocate_memory_bidimensional_matrix_double(mi, nj);

    BI = allocate_memory_bidimensional_matrix_double(nj, mi);
    BD = allocate_memory_bidimensional_matrix_double(nj, mi);
    BC = allocate_memory_bidimensional_matrix_double(nj, mi);

    // ************************************************************

    nonzero_elements = get_nonzero_elements();
    rows_number = (mi-2) * (nj-2);
    columns_number = (mi-2) * (nj-2);
    csr_row_pointer_size = rows_number + 1;

    csr_values          = allocate_memory_vector_double(nonzero_elements);
    csr_column_index    = allocate_memory_vector_int(nonzero_elements);
    csr_row_pointer     = allocate_memory_vector_int(csr_row_pointer_size);
    bx                  = allocate_memory_vector_double(columns_number);
    brh                 = allocate_memory_vector_double(columns_number);
    results             = allocate_memory_vector_double(columns_number);

    values = allocate_memory_vector_double(nonzero_elements);
    column_index = allocate_memory_vector_int(nonzero_elements);
    row_pointer = allocate_memory_vector_int(csr_row_pointer_size);

    // ************************************************************

    /*
        * Se crea la malla 2D
    */
    a = 10.0;
    b = 5.0;

    for (ii = 0; ii <= mi-1; ii++)
    {
        double dii1 = (double)ii;
        double dmi1 = (double)(mi - 1);
        xx[ii] = a * dii1 / dmi1;
    }

    for (jj = 0; jj <= nj-1; jj++)
    {
        double djj1 = (double)jj;
        double dnj1 = (double)(nj -1);
        yy[jj] = b * djj1 / dnj1;
    }
    
    double dmi1 = (double)(mi - 1);
    deltax = 1.0 / dmi1;
    double dnj1 = (double)(nj - 1);
    deltay = 1.0 / dnj1;

    /*
    * Se definen los parametros fisicos del problema
    */
    cond_ter = 100.0;
    temp_ini = 0.0;
    temp_fin = 0.0;
    flux_aba = 5.0;
    flux_arr = 5.0;
    // alpha    = 0.5; //parametro de relajacion

    /*
    * Inicializacion de los arreglos a utilizar
    */
    initialise_bidimensional_matrix_with_double(AI, mi, nj, 0.0);
    initialise_bidimensional_matrix_with_double(AC, mi, nj, 0.0);
    initialise_bidimensional_matrix_with_double(AD, mi, nj, 0.0);
    initialise_bidimensional_matrix_with_double(BI, nj, mi, 0.0);
    initialise_bidimensional_matrix_with_double(BC, nj, mi, 0.0);
    initialise_bidimensional_matrix_with_double(BD, nj, mi, 0.0);
    
    double valor = (temp_fin+temp_ini)/2.0;
    initialise_bidimensional_matrix_with_double(temper, mi, nj, valor);
    initialise_bidimensional_matrix_with_double(temp_ant, mi, nj, valor);

    // *************************************************************************
    for (jj = 1; jj < nj-1; jj++)
    {
        /*
        * Ensamblando matrices en direccion x
        */
        ensambla_tdmax(AI,AC,AD,deltax,deltay,cond_ter,temp_ini,temp_fin,jj);
    }

    for (ii = 1; ii < mi-1; ii++)
    {
        /*
        * Ensamblamos matrices tridiagonales en direccion y
        */
        ensambla_tdmay(BI,BC,BD,deltax,deltay,cond_ter,flux_aba,flux_arr,ii);
    }

    // ****************************************************************************

    get_vector_independent_terms(BI,AI,AD,BD,bx);
    gvit(brh, cond_ter, deltax, deltay);
    get_csr_format_stores(BI, AI, AC, AD, BD, nonzero_elements, csr_values, csr_column_index, csr_row_pointer);
    gcsrfs(values, row_pointer, column_index, nonzero_elements, cond_ter, deltax, deltay);

    // *************************************************************************
    // print_matrix(BI, nj, mi);
//     print_matrix(AI, mi, nj);
    // print_matrix(AC, mi, nj);
    // print_matrix(AD, mi, nj);
    // print_matrix(BD, nj, mi);
    // ****************************************************************************
    // CSR to Dense Format
    double *dense_matrix;
    int matrix_size = rows_number * columns_number;

    dense_matrix = allocate_memory_unidimensional_matrix_double(rows_number, columns_number);
    initialise_vector_with_double(dense_matrix, matrix_size, 0.0);
    // csr_to_dense(csr_values, csr_column_index, csr_row_pointer, dense_matrix, rows_number);
    csr_to_dense(values, column_index, row_pointer, dense_matrix, rows_number);
    // print_flat_matrix(dense_matrix, rows_number);
    // ****************************************************************************
    /* Region del solver    */
    int N = rows_number;
    int NRHS = 1;
    int LDA = N;
    int LDB = NRHS;
    int INFO = 0;
    int *IPIV = allocate_memory_vector_int(N);

    // print_vector(bx, rows_number);
    // INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, dense_matrix, LDA, IPIV, bx, LDB);
    INFO = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, dense_matrix, LDA, IPIV, brh, LDB);

    printf("Info DGESV = %d\n\n", INFO);

    // print_vector(bx, N);

    // ****************************************************************************
    
    print_vector_double(brh, rows_number);
    // print_vector(results, rows_number);
    // print_formato_csr(csr_values, csr_column_index, csr_row_pointer, nonzero_elements, csr_row_pointer_size);
    fill_boundary_conditions_temper_matrix(temper);
    // completar_matriz_temper(temper, bx, rows_number);
    // TODO: Revisar y corregir la manera en como se guardan los valores del resultado obtenido con el solver en la matriz temper y también como se guardan en el archivo de texto.
    completar_matriz_temper(temper, brh, rows_number);
    // print_matrix(temper, mi, nj);
    
    // ****************************************************************************
    // * Escritura de resultados
    // Abrir el archivo para escribir
    FILE *file = fopen("LS.101", "w");

    if (file != NULL) {
        // Utilizar un bucle para imprimir y guardar los datos
        for (ii = 0; ii < mi; ii++) {
            for (jj = 0; jj < nj; jj++)
            {
                fprintf(file, "%f %f %f\n", xx[ii], yy[jj], temper[ii][jj]);
            }
            fprintf(file, "\n");
        }
        // Cerrar el archivo después de escribir
        fclose(file);
    } else {
        printf("Error al abrir el archivo.\n");
    }

    // Liberar memoria de los punteros
    free_memory_vector_int(csr_row_pointer);
    free_memory_vector_int(csr_column_index);
    free_memory_vector_double(csr_values);

    free_memory_bidimensional_matrix_double(BC,nj,mi);
    free_memory_bidimensional_matrix_double(BD,nj,mi);
    free_memory_bidimensional_matrix_double(BI,nj,mi);

    free_memory_bidimensional_matrix_double(AC,mi,nj);
    free_memory_bidimensional_matrix_double(AD,mi,nj);
    free_memory_bidimensional_matrix_double(AI,mi,nj);

    free_memory_vector_double(yy);
    free_memory_vector_double(xx);

    free_memory_bidimensional_matrix_double(temp_ant, mi, nj);
    free_memory_bidimensional_matrix_double(temper, mi, nj);

    return 0;
}
