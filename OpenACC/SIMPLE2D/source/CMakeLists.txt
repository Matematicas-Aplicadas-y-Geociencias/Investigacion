# Version minima de CMake requerido
cmake_minimum_required(VERSION 3.22)

# Nombre del proyecto y lenguaje de programacion a utilizar
project(SIMPLE2D LANGUAGES Fortran)

# Configuracion para usar gpu
option(ENABLE_GPU "Enable GPU use" OFF)

if(ENABLE_GPU)
    message(STATUS "Using OpenAcc GPU")

    # Arquitectura de GPU a utilizar
    set(CUDA_ARCHITECTURE ${COMPUTE_CAPABILITY})
    message(STATUS "CUDA Architecture: ${CUDA_ARCHITECTURE}")

    # Definir el nombre del ejecutable
    set(EXEC_NAME SIMPLE2D-gpu)

    # Banderas para utilizar OpenAcc con GPU
    set(ACC_GPU "-fast -acc=gpu -gpu=${CUDA_ARCHITECTURE} -Minfo=accel")
else()
    message(STATUS "Using OpenAcc MULTICORE")

    # Definir el nombre del ejecutable
    set(EXEC_NAME SIMPLE2D-multicore)

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
    mod_cond_frontera.f90
    mod_ec_energia.f90
    mod_ec_momento.f90
    mod_front_inmersa.f90
    mod_malla.f90
    mod_postproceso.f90
    mod_residuos.f90
    mod_solucionador.f90
    mod_variables.f90
    main.f90)

# Creacion del ejecutable
add_executable(${EXEC_NAME} ${SOURCE_FILES})
