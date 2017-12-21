#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -N 1 -n 1                                                              
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                             
#SBATCH --gres=gpu:1                                                            
cd $SLURM_SUBMIT_DIR

module load cuda
./myODE_cuda 32 8
./myODE_cuda 32 16
./myODE_cuda 32 32
./myODE_cuda 32 64
./myODE_cuda 32 128
./myODE_cuda 64 8
./myODE_cuda 64 16
./myODE_cuda 64 32
./myODE_cuda 64 64
./myODE_cuda 64 128
./myODE_cuda 128 8
./myODE_cuda 128 16
./myODE_cuda 128 32
./myODE_cuda 128 64
./myODE_cuda 128 128
./myODE_cuda 256 8
./myODE_cuda 256 16
./myODE_cuda 256 32
./myODE_cuda 256 64
./myODE_cuda 256 128
./myODE_cuda 512 8
./myODE_cuda 512 16
./myODE_cuda 512 32
./myODE_cuda 512 64
./myODE_cuda 512 128
./myODE_cuda 1024 8
./myODE_cuda 1024 16
./myODE_cuda 1024 32
./myODE_cuda 1024 64
./myODE_cuda 1024 128


mv job_out cuda.txt
mv job_err cuda_err
