#ifndef HERRAMIENTAS_H
#define HERRAMIENTAS_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Error handler: device memory allocation check
void cudaMemoryCheckError(const cudaError_t error)
{
    if (error != cudaSuccess)
    {
        printf("Error: %s:%d, ", __FILE__, __LINE__);
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));
        exit(EXIT_FAILURE);
    }
}

// Allocate memory to an vector of double data type in the device
double *allocateMemoryVectorDevice(const int vector_size)
{
    size_t size_in_bytes; 
    cudaError_t error;
    double *device_pointer;

    size_in_bytes = vector_size * sizeof(double);
    error = cudaMalloc((void **)&device_pointer, size_in_bytes);
    cudaMemoryCheckError(error);
    return device_pointer;
}

// Allocate memory to an vector of integer data type in the device
int *allocateMemoryVectorDeviceInteger(const int vector_size)
{
    size_t size_in_bytes; 
    cudaError_t error;
    int *device_pointer;

    size_in_bytes = vector_size * sizeof(int);
    error = cudaMalloc((void **)&device_pointer, size_in_bytes);
    cudaMemoryCheckError(error);
    return device_pointer;
}

// Data transfer of a vector, from host to device
void transfer_vector_host_to_device(double *host_array, double *device_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(double);
    error = cudaMemcpy(device_array, host_array, size_in_bytes, cudaMemcpyHostToDevice);
    cudaMemoryCheckError(error);
}

// Data transfer of a vector, from host to device
void transfer_vector_host_to_device_integer(int *host_array, int *device_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(int);
    error = cudaMemcpy(device_array, host_array, size_in_bytes, cudaMemcpyHostToDevice);
    cudaMemoryCheckError(error);
}

// Data transfer of a vector, from device to host
void transfer_vector_device_to_host(double *device_array, double *host_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(double);
    error = cudaMemcpy(host_array, device_array, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaMemoryCheckError(error);
}

// Data transfer of a vector, from device to host
void transfer_vector_device_to_host_integer(int *device_array, int *host_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(int);
    error = cudaMemcpy(host_array, device_array, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaMemoryCheckError(error);
}

// Error handler
void error_handler(const char error_string[])
{
    fprintf(stderr, "C language runtime error...\n");
    fprintf(stderr, "%s\n", error_string);
    fprintf(stderr, "...now exiting to system...\n");
    exit(EXIT_FAILURE);
}

// Allocate memory to an vector of int data type
int *allocate_memory_vector_integer(int vector_size)
{
    int *vector;
    size_t size_in_bytes;

    size_in_bytes = vector_size * sizeof(int);
    vector = (int *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        error_handler("Allocation failure in allocate_memory_vector");
    }
    return vector;
}

// Allocate memory to an vector of double data type
double *allocate_memory_vector(int vector_size)
{
    double *vector;
    size_t size_in_bytes;

    size_in_bytes = vector_size * sizeof(double);
    vector = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!vector)
    {
        error_handler("Allocation failure in allocate_memory_vector");
    }
    return vector;
}

// Allocate memory to an matrix of double data type
double *allocate_memory_matrix(int rows, int columns)
{
    double *matrix;
    size_t size_in_bytes;

    size_in_bytes = rows * columns * sizeof(double);
    matrix = (double *)malloc(size_in_bytes);

    // Check that the memory allocation was successful.
    if (!matrix)
    {
        error_handler("Allocation failure in allocate_memory_matrix");
    }
    return matrix;
}

// Print to console a Vector
void print_vector_integer(int *vector, int vector_size)
{
    for (int i = 0; i < vector_size; ++i)
    {
        printf("%d ", vector[i]);
        printf("\n");
    }
    printf("\n");
}

// Print to console a Vector
void print_vector(double *vector, int vector_size)
{
    for (int i = 0; i < vector_size; ++i)
    {
        printf("%f ", vector[i]);
        printf("\n");
    }
    printf("\n");
}

// Print to console a Matrix
void print_matrix(double *matrix, int rows, int columns)
{
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            printf("%f ", matrix[rows * i + j]);
        }
        printf("\n");
    }
    printf("\n");
}

int obtener_nnzMatriz(double *matriz, int rows, int columns)
{
    int nnzMatriz = 0;
    for (int i = 0; i < rows; i++)
    {
        for (int j = 0; j < columns; j++)
        {
            int index = i * rows + j;
            if (matriz[index] != 0.0)
            {
                nnzMatriz++;
            }
        }
    }
    return nnzMatriz;
}

void obtener_formato_csr(double *matriz, int rows, int columns, double *csrValMatriz, int *csrRowPtrMatriz, int *csrColIndMatriz)
{
    int idx = 0;
    for (int i = 0; i < rows; i++)
    {
        csrRowPtrMatriz[i] = idx;
        for (int j = 0; j < columns; j++)
        {
            int index = i * rows + j;
            if (matriz[index] != 0.0)
            {
                csrValMatriz[idx] = matriz[index];
                csrColIndMatriz[idx] = j;
                idx++;
            }
        }
    }
    csrRowPtrMatriz[rows] = idx;
}

#endif