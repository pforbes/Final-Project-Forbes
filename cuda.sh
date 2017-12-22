#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -N 1 -n 1                                                              
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                             
#SBATCH --gres=gpu:1                                                            
cd $SLURM_SUBMIT_DIR

module load cuda
./myODE_cuda 512 128


mv job_out cuda.txt
mv job_err cuda_err
