# MKLPATH=$MKLROOT/lib/em64t
# MKLPATH=$MKLROOT/lib/ia32
MKLPATH=$MKLROOT/lib/intel64

# MKLINCLUDE=$MKLROOT/include/em64t/lp64
# MKLINCLUDE=$MKLROOT/include/ia32
MKLINCLUDE=$MKLROOT/include/intel64/lp64

echo $OMP_NUM_THREADS
# ifort *.f90 -o SIMPLE_MKLS -O2 -L$MKLPATH -I$MKLINCLUDE $MKLPATH/libmkl_lapack95.a -Wl,--start-group $MKLPATH/libmkl_intel.a $MKLPATH/libmkl_sequential.a $MKLPATH/libmkl_core.a -Wl,--end-group -lpthread
ifort *.f90 -o SIMPLE_MKL -O3 -L$MKLPATH -I$MKLINCLUDE $MKLPATH/libmkl_lapack95_lp64.a -Wl,--start-group $MKLPATH/libmkl_intel_lp64.a $MKLPATH/libmkl_intel_thread.a $MKLPATH/libmkl_core.a -Wl,--end-group -openmp -lpthread
# ifort *.f90 -o SIMPLE_MKL -O2 -L$MKLPATH -I$MKLINCLUDE $MKLPATH/libmkl_lapack95.a -Wl,--start-group $MKLPATH/libmkl_intel.a $MKLPATH/libmkl_intel_thread.a $MKLPATH/libmkl_core.a -Wl,--end-group -openmp -lpthread

# ifort -openmp -lpthread *.f90 -Wl,-stack_size,0x10000000 -o SIMPLE
