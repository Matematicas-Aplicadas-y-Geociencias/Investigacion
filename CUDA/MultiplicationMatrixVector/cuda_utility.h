#ifndef CUDA_UTILITY_H
#define CUDA_UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>
void cudaCheckError(const cudaError_t error);

void cudaCheckError(const cudaError_t error)
{
    if (error != cudaSuccess)
    {
        printf("Error: %s:%d, ", __FILE__, __LINE__);
        printf("code:%d, reason: %s\n", error, cudaGetErrorString(error));
        exit(1);
    }
}

double *allocateMemoryVectorDevice(const int vector_size)
{
    size_t size_in_bytes; 
    cudaError_t error;
    double *device_pointer;

    size_in_bytes = vector_size * sizeof(double);
    error = cudaMalloc((void **)&device_pointer, size_in_bytes);
    cudaCheckError(error);
    return device_pointer;
}

double *allocateMemoryMatrixDevice(const int rows, const int columns)
{
    size_t size_in_bytes; 
    cudaError_t error;
    double *device_pointer;

    size_in_bytes = rows * columns * sizeof(double);
    error = cudaMalloc((void **)&device_pointer, size_in_bytes);
    cudaCheckError(error);
    return device_pointer;
}

void transfer_vector_host_to_device(double *host_array, double *device_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(double);
    error = cudaMemcpy(device_array, host_array, size_in_bytes, cudaMemcpyHostToDevice);
    cudaCheckError(error);
}

void transfer_vector_device_to_host(double *device_array, double *host_array, const int array_size)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = array_size * sizeof(double);
    error = cudaMemcpy(host_array, device_array, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaCheckError(error);
}

void transfer_matrix_host_to_device(double *host_array, double *device_array, const int rows, const int columns)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = rows * columns * sizeof(double);
    error = cudaMemcpy(device_array, host_array, size_in_bytes, cudaMemcpyHostToDevice);
    cudaCheckError(error);
}

void transfer_matrix_device_to_host(double *device_array, double *host_array, const int rows, const int columns)
{
    size_t size_in_bytes; 
    cudaError_t error;

    size_in_bytes = rows * columns * sizeof(double);
    error = cudaMemcpy(host_array, device_array, size_in_bytes, cudaMemcpyDeviceToHost);
    cudaCheckError(error);
}

#endif