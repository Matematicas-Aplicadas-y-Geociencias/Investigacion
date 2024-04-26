#include <stdio.h>
#include <stdlib.h>
#include <lapacke.h>

#include "constantes2.h"
#include "herramientas2.h"

int main(int argc, char const *argv[])
{
    // *************************************************************
    // * Declaracion e inicializacion de variables comunes

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
    int xx_size = points_number;

    // Arreglos: matriz del sistema (densa), vector de terminos independientes y vector de resultados
    double *A, *b, *results;
    int A_size = rows_number * columns_number;
    int b_size = columns_number;
    int results_size = columns_number;
    
    // *************************************************************
    // Asignacion de memoria de los arreglos comunes

    temper = allocate_memory_vector(temper_size);
    xx = allocate_memory_vector(xx_size);
    A = allocate_memory_vector(A_size);
    b = allocate_memory_vector(b_size);
    results = allocate_memory_vector(results_size);

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

    cond_ter = 100.0;
    start_temp = 1.0;
    end_temp = 10.0;

    // * Inicializacion de los arreglos a utilizar

    initialise_array_with_value(temper, temper_size, 0.0);
    temper[0] = start_temp;
    temper[temper_size-1] = end_temp;
    initialise_array_with_value(results, results_size, 0.0);
    get_vector_independent_terms(b, b_size, cond_ter, start_temp, end_temp);

    // ****************************************************************************
    // * Variables Y arreglos para usar la libreria LAPACKE

    int N = rows_number;
    int NRHS = 1;
    int LDA = N;
    int LDB = NRHS;

    // int INFO_DGETRF = 0;
    // int INFO_DGETRS = 0;
    // int INFO_DGETRI = 0;
    int INFO_DGESV = 0;

    int *IPIV;
    IPIV = allocate_memory_vector_int(N);

    // Rutinas para usar y ejecutar funciones de la libreria LAPACKE
    get_1D_heat_equation_matrix(A, rows_number, columns_number, cond_ter); 
    
    // INFO_DGETRF = LAPACKE_dgetrf(LAPACK_ROW_MAJOR, N, N, A, LDA, IPIV);
    // INFO_DGETRS = LAPACKE_dgetrs(LAPACK_ROW_MAJOR, 'N', N, NRHS, A, LDA, IPIV, b, LDB);
    // INFO_DGETRI = LAPACKE_dgetri(LAPACK_ROW_MAJOR, N, A, LDA, IPIV);
    INFO_DGESV = LAPACKE_dgesv(LAPACK_ROW_MAJOR, N, NRHS, A, LDA, IPIV, b, LDB);

    // printf("Info DGETRF = %d\n", INFO_DGETRF);
    // printf("Info DGETRS = %d\n", INFO_DGETRS);
    // printf("Info DGETRI = %d\n", INFO_DGETRI);
    printf("Info DGESV = %d\n\n", INFO_DGESV);

    printf("Informacion de la rutina DGESV:\n");
    if (INFO_DGESV == 0)
    {
        printf("¡Solucion encontrada con exito!");
    } 

    if (INFO_DGESV < 0)
    {
        printf("El %d-esimo argumento tiene un valor ilegal.", INFO_DGESV);
    }
    if (INFO_DGESV > 0)
    {
        printf("U(%d, %d) is exactly zero.", INFO_DGESV, INFO_DGESV);
    }
    
    {
        
    }
    
    
    


    // ****************************************************************************
    // * Escritura de results

    // Vaciar resultado obtenido del solver en temper
    get_temper_array(temper, b, start_temp, end_temp);

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

    return 0;
}
