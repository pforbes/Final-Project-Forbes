#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -c 20                                                            
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                                                                         
cd $SLURM_SUBMIT_DIR

./myODE_openMP 10
./myODE_openMP 20

mv job_out openMP.txt
mv job_err omp_err
