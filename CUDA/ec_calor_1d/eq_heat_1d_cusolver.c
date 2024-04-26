#include <stdio.h>
#include <stdlib.h>
#include <cusparse.h>
#include <cusolverSp.h>

#include "constantes.h"
#include "herramientas.h"

int main(int argc, char const *argv[])
{
    // *************************************************************
    // * Declaracion e inicializacion de variables

    // Variables de tamano de la barra
    int points_number = ni;
    double bars_size, deltax;

    // Parametros fisicos del problema
    double cond_ter, start_temp, end_temp;

    // Variables que definen el sistema de ecuaciones lineales a resolver
    int rows_number = points_number - 2;
    int columns_number = points_number - 2;

    // Arreglos para incognitas, termino fuente y posicion
    double *temper;
    double *xx;
    int temper_size = points_number;
    int xx_size     = points_number;

    // Vector de terminos independientes y vector de resultados
    double *b, *results;
    int b_size          = columns_number;
    int results_size    = columns_number;

    // Variables y vectores necesarios para obtener el formato CSR de la matriz del sistema a resolver
    double *csr_values;
    int *csr_col_ind, *csr_row_ptr;

    int nonzero_elements    = get_nonzero_elements();
    int csr_values_size     = nonzero_elements;
    int csr_col_ind_size    = csr_values_size;
    int csr_row_ptr_size    = rows_number + 1;

    // Variables necesarias para utilizar la libreria de cusolver 
    double TOL;
    int REORDER, SINGULARITY;

    cusolverSpHandle_t handle = NULL;
    cusparseMatDescr_t descrMat = NULL;
    cusolverStatus_t status;

    // *************************************************************
    // Asignacion de memoria de los vectores

    temper  = allocate_memory_vector_double(temper_size);
    xx      = allocate_memory_vector_double(xx_size);
    b       = allocate_memory_vector_double(b_size);
    results = allocate_memory_vector_double(results_size);

    csr_values  = allocate_memory_vector_double(csr_values_size);
    csr_col_ind = allocate_memory_vector_int(csr_col_ind_size);
    csr_row_ptr = allocate_memory_vector_int(csr_row_ptr_size);

    // ************************************************************

    // * Se crea la barra 1D
    bars_size = 10.0;

    for (size_t ii = 0; ii <= ni - 1; ii++)
    {
        double dii1 = (double)ii;
        double dni1 = (double)(ni - 1);
        xx[ii] = bars_size * dii1 / dni1;
    }
    
    double dni1 = (double)(ni - 1);
    deltax = 1.0 / dni1;

    // * Se definen los parametros fisicos del problema

    cond_ter    = 100.0;
    start_temp  = 1.0;
    end_temp    = 10.0;

    // * Inicializacion de los vectores a utilizar

    initialise_array_with_value(temper, temper_size, 0.0);
    temper[0] = start_temp;
    temper[temper_size-1] = end_temp;
    initialise_array_with_value(results, results_size, 0.0);
    get_vector_independent_terms(b, b_size, cond_ter, start_temp, end_temp);
    get_csr_schema(cond_ter, nonzero_elements, csr_values, csr_col_ind, csr_row_ptr);

    // ****************************************************************************

    TOL = 1e-5;
    REORDER = 0;
    SINGULARITY = 0;

    cusolverSpCreate(&handle);
    cusparseCreateMatDescr(&descrMat);
    // cusparseSetMatType(descrMat, CUSPARSE_MATRIX_TYPE_GENERAL);
    // cusparseSetMatIndexBase(descrMat, CUSPARSE_INDEX_BASE_ZERO);

    #pragma acc data copyin(b[:b_size], csr_values[:csr_values_size], csr_col_ind[:csr_col_ind_size], csr_row_ptr[:csr_row_ptr_size]) copyout(results[:results_size])
    {
        #pragma acc host_data use_device(b, csr_values, csr_row_ptr, csr_col_ind, results) 
        {

            // status = cusolverSpDcsrlsvchol(handle, rows_number, nonzero_elements, descrMat, csr_values, csr_row_ptr, csr_col_ind, b, TOL, REORDER, results, &SINGULARITY);
            status = cusolverSpDcsrlsvchol(handle, rows_number, nonzero_elements, descrMat, csr_values, csr_row_ptr, csr_col_ind, b, TOL, REORDER, results, &SINGULARITY);

        }
    }

    print_info_cuda_solver(status);

    get_temper_vector(temper, results, start_temp, end_temp);

    // ****************************************************************************
    // * Escritura de results

    // Abrir el archivo para escribir
    FILE *file = fopen("eqHeat1D.101", "w");

    if (file != NULL) {
        // Utilizar un bucle para imprimir y guardar los datos
        for (size_t ii = 0; ii < ni; ii++) {

            fprintf(file, "%f %f\n", xx[ii], temper[ii]);

        }
        // Cerrar el archivo después de escribir
        fclose(file);
    } else {
        printf("Error al abrir el archivo.\n");
    }

    // Liberar memoria de los punteros
    free_memory_vector_int(csr_row_ptr);
    free_memory_vector_int(csr_col_ind);
    free_memory_vector_double(csr_values);
    free_memory_vector_double(results);
    free_memory_vector_double(b);
    free_memory_vector_double(xx);
    free_memory_vector_double(temper);

    return 0;
}