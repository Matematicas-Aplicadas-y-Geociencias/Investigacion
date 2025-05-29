# Version minima de CMake requerido
cmake_minimum_required(VERSION 3.22)

# Nombre del proyecto y lenguaje de programacion a utilizar
project(CALOR2DACC_ARR1D LANGUAGES Fortran)

# Configuracion para usar gpu
option(ENABLE_GPU "Enable GPU use" OFF)

if(ENABLE_GPU)
    message(STATUS "Using OpenAcc GPU")

    # Arquitectura de GPU a utilizar
    set(CUDA_ARCHITECTURE ${COMPUTE_CAPABILITY})
    message(STATUS "CUDA Architecture: ${CUDA_ARCHITECTURE}")

    # Definir el nombre del ejecutable
    set(EXEC_NAME CALOR2DACC_ARR1D-gpu)

    # Banderas para utilizar OpenAcc con GPU
    set(ACC_GPU "-fast -acc=gpu -gpu=${CUDA_ARCHITECTURE} -Minfo=accel")
else()
    message(STATUS "Using OpenAcc MULTICORE")

    # Definir el nombre del ejecutable
    set(EXEC_NAME CALOR2DACC_ARR1D-multicore)

    # Banderas para utilizar OpenAcc con Multicore
    set(ACC_MULTICORE "-fast -acc=multicore -Minfo=accel")
endif()

# Especificar el compilador de fortran a utilizar
set(CMAKE_Fortran_COMPILER nvfortran)

# Definir las banderas que va a utilizar el compilador
if(ENABLE_GPU)
    set(CMAKE_Fortran_FLAGS ${ACC_GPU})
else()
    set(CMAKE_Fortran_FLAGS ${ACC_MULTICORE})
endif()

# Variable para listar los archivos fuentes
set(SOURCE_FILES
    constantes.f90
    ecCalor2D_tdma_parallel_arreglos1D.f90)

# Creacion del ejecutable
add_executable(${EXEC_NAME} ${SOURCE_FILES})
