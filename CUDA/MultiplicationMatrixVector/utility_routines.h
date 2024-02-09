#ifndef UTILITY_ROUTINES_H
#define UTILITY_ROUTINES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// C language standard error handler
void errorHandler(char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
	exit(1);
}

// Allocate an double vector 
double *allocateMemoryVector(int vector_size)
{
    double *vector; 
    size_t sizeInBytes = vector_size * sizeof(double);
    vector = (double *)malloc(sizeInBytes);
    if (!vector)
    {
        errorHandler("Allocation failure in allocateMemoryVector");
    }
    return vector;
}

// Allocate an double matrix 
double *allocateMemoryMatrix(int rows, int columns)
{
    double *matrix;
    size_t sizeInBytes = rows * columns * sizeof(double); 
    matrix = (double *)malloc(sizeInBytes);
    if (!matrix)
    {
        errorHandler("Allocation failure in allocateMemoryMatrix");
    }
    return matrix;
}

// Print a Vector
void printVector(double *vector, int vector_size)
{
    for (int i = 0; i < vector_size ; ++i)
    {
        printf("%f ", vector[i]);
        printf("\n");
    }
    printf("\n");
}

// Print a Matrix
void printMatrix(double *matrix, int rows, int columns)
{
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j) 
        {
            printf("%f ", matrix[rows * i + j]);

        }
        printf("\n");
    }
    printf("\n");
}

/* Build Linear Equation System of test */

void matrixTest(double *matrix, int dimension_matrix)
{
	int nn_element = dimension_matrix * dimension_matrix - 1;
    for (int i = 0; i < dimension_matrix; i++)
    {
        int diagonalIndex = dimension_matrix * i + i;
        for (int j = 0; j < dimension_matrix; j++)
        {
            int index = dimension_matrix * i + j;
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
    for (int i = 0; i < dimension_matrix; i++)
    {
        int index = (dimension_matrix - 1) * (i + 1);
        if (matrix[index] == 0.0) 
        {
            matrix[index] = 0.5;
        }
    }
}

void vectorTest(double *vector, int vector_size)
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