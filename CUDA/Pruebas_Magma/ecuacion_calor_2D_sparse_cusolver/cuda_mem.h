#ifndef CUDA_MEM_H
#define CUDA_MEM_H

#include <stdio.h>
#include <stdlib.h>
#include <cuda_runtime.h>

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
double *allocateMemoryVectorDeviceDouble(const int vector_size)
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
void transfer_vector_host_to_device_double(double *host_array, double *device_array, const int array_size)
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
void transfer_vector_device_to_host_double(double *device_array, double *host_array, const int array_size)
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

void print_memory_info(const char *event)
{
    size_t free_memory;
    size_t total_memory;
    cudaError_t cuda_error;

    cuda_error = cudaMemGetInfo(&free_memory, &total_memory);

    cudaMemoryCheckError(cuda_error);

    printf("\n%s: \n", event);
    printf("Memoria libre = %.2f GB\n", free_memory / (1024.0 * 1024 * 1024.0));
    printf("Memoria Total = %.2f GB\n", total_memory / (1024.0 * 1024.0 * 1024.0));
}

#endif