#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -c 1                                                            
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                                                                         
cd $SLURM_SUBMIT_DIR

./myODE_seq
./myODE_vec

mv job_out seq_vec.txt
mv job_err myODE_err
