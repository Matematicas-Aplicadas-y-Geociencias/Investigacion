#ifndef CONSTANTES_H
#define CONSTANTES_H

#include <string.h>

const int mi = 1024, nj = 1024, nn = 1024;
const double pi = 3.1415926535;

void ensambla_tdmax(
    double **AI,
    double **AC,
    double **AD,
    double **resultx,
    const double deltax,
    const double deltay,
    const double **temp_ant,
    const double cond_ter,
    const double temp_ini,
    const double temp_fin,
    const int jj)
{
    int iin;

    /*
    * Se definen las condiciones de frontera
    */
    AC[0][jj] = 1.0;
    AD[0][jj] = 0.0;
    resultx[0][jj] = temp_ini;
    AI[mi-1][jj] = 0.0;
    AC[mi-1][jj] = 1.0;
    resultx[mi-1][jj] = temp_fin;
    /*
    * Ensamblado de la matriz tridiagonal y del vector de resultados
    */
    for (iin = 1; iin < mi-1; iin++)
    {
        AI[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        AC[iin][jj] = 2.0 * cond_ter * (1.0 / (deltax * deltax) + 1.0 / (deltay * deltay));
        AD[iin][jj] = -1.0 * cond_ter / (deltax * deltax);
        resultx[iin][jj] = cond_ter / (deltay * deltay) * temp_ant[iin][jj+1] + cond_ter / (deltay * deltay) * temp_ant[iin][jj-1]; 
    }
}

void ensambla_tdmay(
    double **BI,
    double **BC,
    double **BD,
    double **resulty,
    const double deltax,
    const double deltay,
    const double **temper,
    const double cond_ter,
    const double flux_aba,
    const double flux_arr,
    const int ii
)
{
    int jjn;

    /*
    * Se definen las condiciones de frontera
    */
    BC[0][ii]         = 1.0;    // -1.0 / deltay;
    BD[0][ii]         = 0.0;    // 1.0 / deltay;
    resulty[0][ii]    = 308.0;  // flux_aba;
    BI[nj-1][ii]      = 0.0;    // -1.0 / deltay;
    BC[nj-1][ii]      = 1.0;    // 1.0 / deltay;
    resulty[nj-1][ii] = 308.0;  // flux_arr
    /*
    * Ensamblado de la matriz tridiagonal y del vector de resultados
    */
    for (jjn = 1; jjn < nj-1; jjn++)
    {
        BI[jjn][ii]      = -1.0 * cond_ter / (deltay * deltay);
        BC[jjn][ii]      = 2.0 * cond_ter * (1.0 / (deltay * deltay) + 1.0 / (deltax * deltax));
        BD[jjn][ii]      = -1.0 * cond_ter / (deltay * deltay);
        resulty[jjn][ii] = cond_ter / (deltax * deltax) * temper[ii+1][jjn] + cond_ter / (deltax * deltax) * temper[ii-1][jjn]; 
    }
}

void tri(
    const double *a,
    double *b,
    const double *c,
    double *r,
    const int n
)
{
    int i;

    // Eliminacion elementos bajo la matriz - Eliminacion Gaussiana
    for (i = 1; i < n; i++)
    {
        r[i] = r[i] - (a[i] / b[i-1]) * r[i-1];
        b[i] = b[i] - (a[i] / b[i-1]) * c[i-1];
    }
    // Solucion para r - Substitucion hacia atras
    r[n-1] = r[n-1] / b[n-1];
    for (i = n-2; i >= 0; i--)
    {
        r[i] = (r[i] - c[i] * r[i+1]) / b[i];
    }
}

#endif