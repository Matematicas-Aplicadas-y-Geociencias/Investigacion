cmake_minimum_required(VERSION 3.22)
project(Integrador3D LANGUAGES Fortran)
add_executable(integrador main.f90 sp_integrador_2d.f90 get_coef.f90 integrador_3d.f90 mod_constantes.f90 function.f90)
