gfortran -fopenmp mod_malla.f90 mod_ec_continuidad.f90 mod_ec_momento.f90 mod_ec_energia.f90 mod_solucionador.f90 mod_postproceso.f90 IXCHEL.f90 -o ixchel3D -ffree-line-length-512 -O3 

#gfortran -fopenmp -foffload=-lgfortran -fopenmp  mod_malla.f90 mod_ec_continuidad.f90 mod_ec_momento.f90 mod_ec_energia.f90 mod_solucionador.f90 mod_postproceso.f90 IXCHEL.f90 -o ixchel3D -O3
#nvfortran -mp=gpu -gpu=cc89 mod_malla.f90 mod_ec_continuidad.f90 mod_ec_momento.f90 mod_ec_energia.f90 mod_solucionador.f90 mod_postproceso.f90 IXCHEL.f90 -o ixchel3D -O3 -Minfo=accel
