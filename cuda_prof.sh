#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -N 1 -n 1                                                              
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                             
#SBATCH --gres=gpu:1                                                            
cd $SLURM_SUBMIT_DIR

module load cuda
nvprof -o profile.out -f ./myODE_cuda_prof 128 128
nvprof -i profile.out --print-gpu-trace

mv job_out cuda_out.txt
mv job_err cuda_prof.txt
