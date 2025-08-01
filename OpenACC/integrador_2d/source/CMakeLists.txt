cmake_minimum_required(VERSION 3.22)
project(Integrador3D LANGUAGES Fortran)
add_executable(integrador main.f90 integrador_2D.f90 integrate_2D.f90 mod_constantes.f90 functions.f90)
