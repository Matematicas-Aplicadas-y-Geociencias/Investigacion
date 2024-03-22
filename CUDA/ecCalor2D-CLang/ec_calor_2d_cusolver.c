#include <stdio.h>
#include <stdlib.h>
#include <cusparse.h>
#include <cusolverSp.h>

#include "constantes2.h"
#include "herramientas2.h"

int main(int argc, char const *argv[])
{
    // *************************************************************
    // Declaracion de variables
    int ii, jj, kk;
    
    // Variables de tamano de la placa
    double a, b, deltax, deltay;

    // Parametros fisicos del problema
    double cond_ter, temp_ini, temp_fin;
    double flux_aba, flux_arr, alpha;

    // Arreglos para incognitas, termino fuente y posicion
    double **temper, **temp_ant;
    double *xx;
    double *yy;

    //Matriz a invertir
    double **AI, **AD, **AC;
    double **BI, **BD, **BC;

    // CSR Sparse Format
    double *csr_valores;
    int *csr_col_ind;
    int *csr_ptr;
    double *resultados;

    int numero_elementos_no_cero;
    int tamanio_matriz_completa;
    int tamanio_csr_ptr;

    double TOL = 1e-6;
    int REORDER = 0;
    int singularity = 0;

    cusolverSpHandle_t handle = NULL;
    cusparseMatDescr_t descrA = NULL;
    
    // *************************************************************
    
    temper = allocate_memory_matrix(mi, nj);
    temp_ant = allocate_memory_matrix(mi, nj);

    xx = allocate_memory_vector(mi);
    yy = allocate_memory_vector(nj);

    
    AI = allocate_memory_matrix(mi, nj);
    AD = allocate_memory_matrix(mi, nj);
    AC = allocate_memory_matrix(mi, nj);

    BI = allocate_memory_matrix(nj, mi);
    BD = allocate_memory_matrix(nj, mi);
    BC = allocate_memory_matrix(nj, mi);

    // ************************************************************

    numero_elementos_no_cero = obtener_total_elementos_no_cero();
    tamanio_matriz_completa = (mi-2) * (nj-2);
    tamanio_csr_ptr = tamanio_matriz_completa + 1;

    csr_valores = allocate_memory_vector(numero_elementos_no_cero);
    csr_col_ind = allocate_memory_vector_int(numero_elementos_no_cero);
    csr_ptr = allocate_memory_vector_int(tamanio_csr_ptr);
    resultados = allocate_memory_vector(tamanio_matriz_completa);

    // ************************************************************
    
    cusolverSpCreate(&handle);
    cusparseCreateMatDescr(&descrA);
    cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL); 
    cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

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
    temp_ini = 308.0;
    temp_fin = 298.0;
    flux_aba = 5.0;
    flux_arr = 5.0;
    alpha    = 0.5; //parametro de relajacion

    /*
    * Inicializacion de los arreglos a utilizar
    */
    inicializar_matriz(AI, mi, nj, 0.0);
    inicializar_matriz(AC, mi, nj, 0.0);
    inicializar_matriz(AD, mi, nj, 0.0);
    inicializar_matriz(BI, nj, mi, 0.0);
    inicializar_matriz(BC, nj, mi, 0.0);
    inicializar_matriz(BD, nj, mi, 0.0);
    
    double valor = (temp_fin+temp_ini)/2.0;
    inicializar_matriz(temper, mi, nj, valor);
    inicializar_matriz(temp_ant, mi, nj, valor);

    /*
    * Abrimos la region de datos paralela
    */
    #pragma acc data copy(temper[:mi][:nj]) \
    copyin(deltax,deltay,temp_ant[:mi][:nj],cond_ter,temp_ini,\
    temp_fin,flux_aba,flux_arr,alpha) \
    create(AI[:mi][:nj],AC[:mi][:nj],AD[:mi][:nj],\
    BI[:nj][:mi],BC[:nj][:mi],BD[:nj][:mi])
    {
        /*
        * Inicia el ciclo que recorre la coordenada y resolviendo problemas 1D en la direccion de x
        */
        #pragma acc parallel loop
        for (jj = 1; jj < nj-1; jj++)
        {
            /*
            * Ensamblando matrices en direccion x
            */
            ensambla_tdmax(AI,AC,AD,deltax,deltay,cond_ter,temp_ini,temp_fin,jj);
        }

        /*
        * Inicia el ciclo que recorre la coordenada x resolviendo problemas 1D en la direccion de y
        */
        #pragma acc parallel loop
        for (ii = 1; ii < mi-1; ii++)
        {
            /*
            * Ensamblamos matrices tridiagonales en direccion y
            */
            ensambla_tdmay(BI,BC,BD,deltax,deltay,cond_ter,flux_aba,flux_arr,ii);
        }
            
    /*
    * Cerramos la region de datos paralela
    */
    }

    // ****************************************************************************
    
    obtener_vector_terminos_independientes(BI,AI,AD,BD,resultados);
    obtener_formato_csr(BI, AI, AC, AD, BD, numero_elementos_no_cero, csr_valores, csr_col_ind, csr_ptr);
    // print_vector(resultados, tamanio_matriz_completa);
    // print_formato_csr(csr_valores, csr_col_ind, csr_ptr, numero_elementos_no_cero, tamanio_csr_ptr);
    
    // ****************************************************************************
    /*
    * Escritura de resultados
    */
    // Abrir el archivo para escribir
    // FILE *file = fopen("clang.101", "w");

    // if (file != NULL) {
    //     // Utilizar un bucle para imprimir y guardar los datos
    //     for (ii = 0; ii < mi; ii++) {
    //         for (jj = 0; jj < nj; jj++)
    //         {
    //             fprintf(file, "%f %f %f\n", xx[ii], yy[jj], temper[ii][jj]);
    //         }
    //         fprintf(file, "\n");
    //     }
    //     // Cerrar el archivo después de escribir
    //     fclose(file);
    // } else {
    //     printf("Error al abrir el archivo.\n");
    // }
    

    // Liberar memoria de los punteros
    free(csr_ptr);
    free(csr_col_ind);
    free(csr_valores);

    free_matrix(BC,nj,mi);
    free_matrix(BD,nj,mi);
    free_matrix(BI,nj,mi);

    free_matrix(AC,mi,nj);
    free_matrix(AD,mi,nj);
    free_matrix(AI,mi,nj);

    free(yy);
    free(xx);

    free_matrix(temp_ant, mi, nj);
    free_matrix(temper, mi, nj);

    return 0;
}
