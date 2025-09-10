nvfortran -fast -acc=gpu -gpu=cc75,managed -Minfo=all *.f90 -o test
nvfortran -mp=gpu -gpu=cc120 constantes.f90 ecCalor2D_tdma_parallel.f90 -o ec2dacc -Minfo=accel
