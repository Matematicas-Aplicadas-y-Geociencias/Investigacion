#include <stdio.h>
#include <stdlib.h>

#include "constantes.h"
#include "herramientas.h"

int main(int argc, char const *argv[])
{
    int ii, jj, kk;
    
    // Variables de tamano de la placa
    double a, b, deltax, deltay;

    // Parametros fisicos del problema
    double cond_ter, temp_ini, temp_fin;
    double flux_aba, flux_arr, alpha;

    // Arreglos para incognitas, termino fuente y posicion
    double temper[mi][nj], temp_ant[mi][nj];
    double sourcex[mi], xx[mi];
    double sourcey[nj], yy[nj];
    double resultx[mi][nj], tempx[mi][nj];
    double resulty[nj][mi], tempy[nj][mi];

    //Matriz a invertir
    double AI[mi][nj], AD[mi][nj], AC[mi][nj];
    double BI[nj][mi], BD[nj][mi], BC[nj][mi];
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

    inicializar_matriz(resultx, mi, nj, 0.0);
    inicializar_matriz(resulty, nj, mi, 0.0);
    inicializar_matriz(tempx, mi, nj, 0.0);
    inicializar_matriz(tempy, nj, mi, 0.0);
    /*
    * Bucle de pseudotiempo
    */
    for (kk = 0; kk < 30; kk++)
    {
        /*
        * Inicia el ciclo que recorre la coordenada y resolviendo problemas 1D en la direccion de x
        */
       for (jj = 1; jj < nj-1; jj++)
       {
            /*
            * Ensamblando matrices en direccion x
            */
            ensambla_tdmax(AI,AC,AD,resultx,deltax,deltay,temp_ant,cond_ter,temp_ini,temp_fin,jj);
       }
       /*
       * Llamamos al resolvedor
       */
        for (jj = 1; jj < nj-1; jj++)
        {
            // TODO: Crear funcion para tomar columnas de una matriz como vectores
            // * Llamar a la funcion tri
            for (ii = 0; ii < mi; ii++)
            {
                temper[ii][jj] = resultx[ii][jj];
            }
            
        }
        /*
        * Inicia el ciclo que recorre la coordenada x resolviendo problemas 1D en la direccion de y
        */
        for (ii = 1; ii < mi-1; ii++)
        {
            /*
            * Ensamblamos matrices tridiagonales en direccion y
            */
           ensambla_tdmay(BI,BC,BD,resulty,deltax,deltay,temper,cond_ter,flux_aba,flux_arr,ii);
        }
        /*
        * Llamamos al resolvedor
        */
        for (ii = 1; ii < mi-1; ii++)
        {
            // TODO: Crear funcion para tomar columnas de una matriz como vectores
            // * Llamar a la funcion tri
            for (jj = 0; jj < nj; jj++)
            {
                temper[ii][jj] = resulty[jj][ii];
            }
        }
        /*
        * Se actualiza la temperatura de la iteracion anterior
        */
        for (ii = 0; ii < mi; ii++)
        {
            for (jj = 0; jj < nj; jj++)
            {
                temp_ant[ii][jj] = temper[ii][jj];
            }
        }
        
    }
    /*
    * Escritura de resultados
    */
   // Abrir el archivo para escribir
    FILE *file = fopen("fort.101", "w");

    if (file != NULL) {
        // Utilizar un bucle para imprimir y guardar los datos
        for (ii = 0; ii < mi; ii++) {
            for (jj = 0; jj < nj; jj++)
            {
                fprintf(file, "%f %f %f\n", xx[ii], yy[ii], temper[ii][jj]);
            }
        }
        // Cerrar el archivo después de escribir
        fclose(file);
    } else {
        printf("Error al abrir el archivo.\n");
    }
    return 0;
}
