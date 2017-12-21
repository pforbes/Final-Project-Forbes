#!/bin/bash                                                                     
#SBATCH --partition=slurm_shortgpu                                             
#SBATCH -c 32                                                            
#SBATCH -o job_out                                                             
#SBATCH -e job_err                                                                                                         
cd $SLURM_SUBMIT_DIR

./myODE_openMP 1
./myODE_openMP 2
./myODE_openMP 3
./myODE_openMP 4
./myODE_openMP 5
./myODE_openMP 6
./myODE_openMP 7
./myODE_openMP 8
./myODE_openMP 9
./myODE_openMP 10
./myODE_openMP 11
./myODE_openMP 12
./myODE_openMP 13
./myODE_openMP 14
./myODE_openMP 15
./myODE_openMP 16
./myODE_openMP 17
./myODE_openMP 18
./myODE_openMP 19
./myODE_openMP 20
./myODE_openMP 21
./myODE_openMP 22
./myODE_openMP 23
./myODE_openMP 24
./myODE_openMP 25
./myODE_openMP 26
./myODE_openMP 27
./myODE_openMP 28
./myODE_openMP 29
./myODE_openMP 30
./myODE_openMP 31
./myODE_openMP 32
./myODE_openMP 33
./myODE_openMP 34
./myODE_openMP 35
./myODE_openMP 36
./myODE_openMP 37
./myODE_openMP 38
./myODE_openMP 39
./myODE_openMP 40

mv job_out openMP.txt
mv job_err omp_err
