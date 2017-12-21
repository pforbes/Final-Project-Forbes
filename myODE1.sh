#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -c 1                                                            
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                                                                         
cd $SLURM_SUBMIT_DIR

./myODE1

mv job_out myODE1.txt
mv job_err myODE_err
