#ifndef HERRAMIENTAS_H
#define HERRAMIENTAS_H

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