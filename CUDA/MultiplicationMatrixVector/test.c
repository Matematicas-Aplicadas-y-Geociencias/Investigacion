#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
#include <cublas.h>

#include "cuda_utility.h"
#include "utility_routines.h"

int main(int argc, char const *argv[])
{
    double *matrixA;
    matrixA = allocateMemoryMatrix(8, 8);

    double *vectorB;
    vectorB = allocateMemoryVector(8);

    vectorTest(vectorB, 8);
    matrixTest(matrixA, 8);

    printMatrix(matrixA, 8, 8);
    printVector(vectorB, 8);
    return 0;
}
