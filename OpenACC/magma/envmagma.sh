export NVHPC_ROOT=/opt/nvidia/hpc_sdk/Linux_x86_64/25.3
export CUDA_ROOT=$NVHPC_ROOT/cuda
export MAGMA_ROOT=/usr/local

export PATH=$NVHPC_ROOT/compilers/bin:$CUDA_ROOT/bin:$PATH
export LD_LIBRARY_PATH=$MAGMA_ROOT/lib:$CUDA_ROOT/lib64:$NVHPC_ROOT/compilers/lib:$LD_LIBRARY_PATH
