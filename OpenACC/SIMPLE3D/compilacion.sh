nvfortran -acc=multicore  mod_malla.f90 mod_ec_continuidad.f90 mod_ec_momento.f90 mod_ec_energia.f90 mod_solucionador.f90 mod_postproceso.f90 main.f90 -o ixchel3D -Minfo=accel,acc
#gfortran -fopenmp mod_malla.f90 mod_ec_continuidad.f90 mod_ec_momento.f90 mod_ec_energia.f90 mod_solucionador.f90 mod_postproceso.f90 main.f90 -o test -ffree-line-length-512 -O2 -fbounds-check
