# Final-Project-Forbes

Directions for running code

1. compile files using the Makfile
2. Run desired batch.sh files using sbatch as described below
   - Programs will be run on short GPU partition with required settings
   - output will be sent to txt file in the working directory
   - format and name of output file is described for each batch file below

***PROGRAMS***

1.  myODE1.sh
    -- instantaneous run time, requires 1 CPU-- 
    - runs myODE1 
      - Used to generate Figure 1 and Figure 2
      - generates Gain vs Length and all ODE data for one {b0,C} input pair
    --> output --> myODE1.txt
    	format -->
	       <Length(x)> <real(f(x),imag(f(x))> <real(f'(x)),imag(f'(x))> <real(f"(x)),imag(f"(x))> <gain> <gain[dB]>

2. sequential.sh //probably not worth your time for 2 lines of output//
   -- 7 minutes or more, requires 1 CPU --
   - used to generate data for Report section 3.1  
   - runs 2 files
     	  -myODE_seq  --> the sequential version hardcoded to solve 5500 iterations
	  -myODE_vec  --> the "vectorized" version solves the same 5500 iterations
   --> output --> seq_vec.txt
       format -- Two lines--

3. openMP.sh   
   -- 40 second run time, requires 20 CPU cores --
   - sample of what was used to generate Figure 5 and Figure 6
   - runs myODE_openMP to solve 5500 iterations using 10 cores and 20 cores
   --> output --> openMP.txt
       format -->
              <n_processors> <Time[ms]> <iterations> <time/iteration> <time[ms]>

4. cuda.sh  
   -- 40 second run time, requires 1 GPU --
   - sample of what was used to generate Figure 7 and Figure 8
   - stores only timing inormation
   - runs myODE_cuda 512 128
        - for 128000 iterations using 512 thread/block
   --> output --> cuda.txt
       format -->
              <Time[ms]> <iterations> <time/iteration[ms]> <threads/block>

5. cuda_prof.sh
   -- 40 second run time, requires 1 GPU --
   - used to generate profile information, Figure 3, and Figure 4
   - runs profiler nvprof and stores {b0,C} sweep output
   - uses myODE_cuda_prof 128 128
     - set up for 128000 iterations and 128 threads/block
   --> output (sweep output) --> cuda_out.txt
       format -->
     	     <C> <b0> <Max-Gain-Length> <Max-Gain-[dB]>
   --> output (profiler output) --> cuda_prof.txt
       	      < self explanatory > 