#ifndef CONSTANTES2_H
#define CONSTANTES2_H

#include "herramientas2.h"

const int mi = 6, nj = 5, nn = 1024;
const double pi = 3.1415926535;

#pragma acc routine vector
void ensambla_tdmax(
    double **AI,
    double **AC,
    double **AD,
    double **resultx,
    const double deltax,
    const double deltay,
    double **temp_ant,
    const double cond_ter,
    const double temp_ini,
    const double temp_fin,
    const int jj)
{
    int iin;

    /*
    * Se definen las condiciones de frontera
    */
    AI[0][jj] = 2.0;
    // AC[0][jj] = 1.0;
    AD[0][jj] = 0.0;
    resultx[0][jj] = temp_ini;
    AI[mi - 1][jj] = 0.0;
    // AC[mi - 1][jj] = 1.0;
    resultx[mi - 1][jj] = temp_fin;
    /*
    * Ensamblado de la matriz tridiagonal y del vector de resultados
    */
    #pragma acc loop vector
    for (iin = 1; iin < mi - 1; iin++)
    {
        AI[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        AC[iin][jj] = 2.0 * cond_ter * (1.0 / (deltax * deltax) + 1.0 / (deltay * deltay));
        AD[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        resultx[iin][jj] = cond_ter / (deltay * deltay) * temp_ant[iin][jj + 1] + cond_ter / (deltay * deltay) * temp_ant[iin][jj - 1];
    }
}

#pragma acc routine vector
void ensambla_tdmay(
    double **BI,
    double **BC,
    double **BD,
    double **resulty,
    const double deltax,
    const double deltay,
    double **temper,
    const double cond_ter,
    const double flux_aba,
    const double flux_arr,
    const int jj)
{
    int jjn;

    /*
     * Se definen las condiciones de frontera
     */
    BI[0][jj] = 1.0;
    // BC[0][jj] = 1.0;             // -1.0 / deltay;
    BD[0][jj] = 0.0;             // 1.0 / deltay;
    resulty[0][jj] = 308.0;      // flux_aba;
    BI[nj - 1][jj] = 0.0;        // -1.0 / deltay;
    BD[nj - 1][jj] = 3.0;        // 1.0 / deltay;
    resulty[nj - 1][jj] = 308.0; // flux_arr
    /*
    * Ensamblado de la matriz tridiagonal y del vector de resultados
    */
    #pragma acc loop vector
    for (jjn = 1; jjn < nj - 1; jjn++)
    {
        BI[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
        BC[jjn][jj] = 2.0 * cond_ter * (1.0 / (deltay * deltay) + 1.0 / (deltax * deltax));
        BD[jjn][jj] = -1.0 * cond_ter / (deltay * deltay);
        resulty[jjn][jj] = cond_ter / (deltax * deltax) * temper[jj + 1][jjn] + cond_ter / (deltax * deltax) * temper[jj - 1][jjn];
    }
}

#pragma acc routine
void tri(
    double **a,
    double **b,
    double **c,
    double **r,
    const int filas,
    const int columna)
{
    int i;

    // Eliminacion elementos bajo la matriz - Eliminacion Gaussiana
    for (i = 1; i < filas; i++)
    {
        r[i][columna] = r[i][columna] - (a[i][columna] / b[i - 1][columna]) * r[i - 1][columna];
        b[i][columna] = b[i][columna] - (a[i][columna] / b[i - 1][columna]) * c[i - 1][columna];
    }
    // Solucion para r - Substitucion hacia atras
    r[filas - 1][columna] = r[filas - 1][columna] / b[filas - 1][columna];
    for (i = filas - 2; i >= 0; i--)
    {
        r[i][columna] = (r[i][columna] - c[i][columna] * r[i + 1][columna]) / b[i][columna];
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
                b[kk] += BI[jj-1][ii];
            }

            if (ii == 1)
            {
                b[kk] += AI[ii-1][jj];
            }

            if (ii == mi-2)
            {
                b[kk] += AD[ii+1][jj];
            }

            if (jj == nj-2)
            {
                b[kk] += BD[jj+1][ii];
            }
            
            kk++;
        }
        
    }
    
}
#endif