#!/bin/bash
#SBATCH --partition=64crs
#SBATCH --job-name=cil100
#SBATCH --output=simpout.out
#SBATCH --error=simpout.err
#SBATCH --nodes=1
# #SBATCH --ntasks=64
# #SBATCH --mem=100G
#SBATCH --ntasks-per-node=32
#SBATCH --time=02:30:00
#SBATCH --nodelist=nodo[02]
# srun --mpi=pmi2 ./Alya.x cilindro

ulimit -s unlimited
time ./SIMPLE3D 

