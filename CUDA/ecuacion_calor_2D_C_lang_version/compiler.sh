nvc -fast -acc=gpu -gpu=cc75 -Minfo=all *.c -o test
nvc -acc=gpu -gpu=cc75 -Minfo=all ec_calor_2d_cusolver.c -o calor2dacc -lcudart -lcublas -lcusparse -lcusolver
