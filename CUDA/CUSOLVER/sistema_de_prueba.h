#ifndef SISTEMA_DE_PRUEBA_H
#define SISTEMA_DE_PRUEBA_H

#include <math.h>

// Matriz de coeficientes
void crear_matriz_de_prueba(double *matrix, int rows, int columns)
{
	int nn_element = rows * columns - 1;
    for (int i = 0; i < rows; i++)
    {
        int diagonalIndex = rows * i + i;
        for (int j = 0; j < columns; j++)
        {
            int index = rows * i + j;
            if (index == diagonalIndex)
            {
                matrix[index] = 3.0;
                if (index == 0) 
                {
                    matrix[index+1] = -1.0;
                    j = index + 1;
                }
                else if (index == nn_element) 
                {
                    matrix[index-1] = -1.0;
                }
                else 
                {
                    matrix[index-1] = -1.0;
                    matrix[index+1] = -1.0;
                    j = index + 1;
                }
            }
            else
            {
                matrix[index] = 0.0;
            }
        }
    }
    for (int i = 0; i < rows; i++)
    {
        int index = (rows - 1) * (i + 1);
        if (matrix[index] == 0.0) 
        {
            matrix[index] = 0.5;
        }
    }
}

// Vector del lado derecho
void crear_vector_de_prueba(double *vector, int vector_size)
{
	int halfVector = (int)floor(vector_size/2);

    vector[0] = 2.5;
	vector[vector_size - 1] = 2.5;
	vector[halfVector - 1] = 1.0;
	vector[halfVector] = 1.0;
    
    for (int i = 1; i <= halfVector-2; i++)
    {
        vector[i] = 1.5;
	}
    
    for (int i = halfVector+1; i <= vector_size-2; i++)
    {
        vector[i] = 1.5;
	}
}   

#endif