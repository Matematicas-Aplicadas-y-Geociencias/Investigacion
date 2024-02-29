#!/bin/bash
#SBATCH --partition=64crs
#SBATCH --job-name=SIMPLE
#SBATCH --output=simpout.out
#SBATCH --error=simpout.err
#SBATCH --nodes=1
# #SBATCH --ntasks=64
# #SBATCH --mem=100G
#SBATCH --ntasks-per-node=24
#SBATCH --time=00:50:00
# #SBATCH --nodelist=nodo[02]
# srun --mpi=pmi2 ./Alya.x cilindro

export OMP_NUM_THREADS=32
ulimit -s unlimited
time ./SIMPLE2D 

