#ifndef CONSTANTES2_H
#define CONSTANTES2_H

// const int mi = 6, nj = 5;
const int ni = 14;
const double pi = 3.1415926535;

#include "herramientas2.h"

// #pragma acc routine vector
// void ensambla_tdmax(
//     double **AI,
//     double **AC,
//     double **AD,
//     const double deltax,
//     const double deltay,
//     const double cond_ter,
//     const double temp_ini,
//     const double temp_fin,
//     const int jj)
// {
//     int iin;

//     /*
//     * Se definen las condiciones de frontera
//     */
//     AI[0][jj] = -0.0* cond_ter / (deltax * deltax);
//     // AD[0][jj] = 0.0;
//     AD[mi - 1][jj] = -0.0* cond_ter / (deltax * deltax);
//     /*
//     * Ensamblado de la matriz tridiagonal y del vector de resultados
//     */
//     #pragma acc loop vector
//     for (iin = 1; iin < mi - 1; iin++)
//     {
//         AI[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
//         AC[iin][jj] = 2.0 * cond_ter * (1.0 / (deltax * deltax) + 1.0 / (deltay * deltay));
//         AD[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
//     }
// }

// #pragma acc routine vector
// void ensambla_tdmay(
//     double **BI,
//     double **BC,
//     double **BD,
//     const double deltax,
//     const double deltay,
//     const double cond_ter,
//     const double flux_aba,
//     const double flux_arr,
//     const int jj)
// {
//     int jjn;

//     /*
//      * Se definen las condiciones de frontera
//      */
//     BI[0][jj] = -0.0 * cond_ter/ (deltay * deltay);
//     BD[nj - 1][jj] = -3.0* cond_ter / (deltay * deltay);        // 1.0 / deltay;
//     /*
//     * Ensamblado de la matriz tridiagonal y del vector de resultados
//     */
//     #pragma acc loop vector
//     for (jjn = 1; jjn < nj - 1; jjn++)
//     {
//         BI[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
//         BD[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
//     }
// }

// #pragma acc routine
// void tri(
//     double **a,
//     double **b,
//     double **c,
//     double **r,
//     const int filas,
//     const int columna)
// {
//     int i;

//     // Eliminacion elementos bajo la matriz - Eliminacion Gaussiana
//     for (i = 1; i < filas; i++)
//     {
//         r[i][columna] = r[i][columna] - (a[i][columna] / b[i - 1][columna]) * r[i - 1][columna];
//         b[i][columna] = b[i][columna] - (a[i][columna] / b[i - 1][columna]) * c[i - 1][columna];
//     }
//     // Solucion para r - Substitucion hacia atras
//     r[filas - 1][columna] = r[filas - 1][columna] / b[filas - 1][columna];
//     for (i = filas - 2; i >= 0; i--)
//     {
//         r[i][columna] = (r[i][columna] - c[i][columna] * r[i + 1][columna]) / b[i][columna];
//     }
// }
#endif