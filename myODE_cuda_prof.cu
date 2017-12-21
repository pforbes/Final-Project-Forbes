// INPUT is ./ODE numThreads N*501
#include <iostream>
#include <vector>
#include "cuda_complex.hpp"
#include <math.h>
#include <random>
#include <chrono>
#include "stopwatch.hpp"
#include <cuda.h>
	
const double b_sigma = 0.001; // normal distribution std deviation for b
/* set b0 values to increment: # points = (MAX-MIN)/INCR+1 */
const float MIN_b0 = 1.5125;
const float MAX_b0 = 1.5225;
const float INCR_b0 =0.00001;
const int SIZE_b0 = int((MAX_b0-MIN_b0)/INCR_b0)+1;

/* set C values to increment: # points = (MAX-MIN)/INCR+1 */
const float MIN_C = 0.0175;
const float MAX_C = 0.0225;
double INCR_C= 0.0001;
int SIZE_C = int((MAX_C-MIN_C)/INCR_C)+1;
  
/* Set resolution in length along device ('time' for ODE solvers) */
const float L_max = 150.0;
//const float L_inc = 0.001;

const int B_LENGTH = int(L_max) + 1;  // length of the b vector
//const double d = 0;
//const double QC = 1.85*0.25;


__global__ void ODE_Kernel(float* db0_in, float* dC_in,float *db_rand, float* dLength,double* dGain, int num_loops, int b_length)
{
  double dt = 0.001;
  int idx = blockIdx.x*blockDim.x + threadIdx.x;
  
  //
  // could implement shared memory method for interpolating b(x)
  //
  
  if (idx < num_loops)
    {
      double b0 = double(db0_in[idx]);
      double C  = double(dC_in[idx]);
      double QC = 0.25*1.85;  // change space charge constant consistent with changes in C but keeping Q fixed
      double d = 0.0;
      
      //state initialization
      complex<double> x0(0.0,0.0);
      complex<double> x1(0.0,0.0);
      complex<double> x2(1.0,0.0);

      complex<double> i(0,1);

      // attempted to minimize 150000*4*4 repeated multiplications but made no diffrence or took longer
      //complex<double> const1 = -i*C*C*C;
      
      // Solve ODE using
	for(int n = 0; n < 150000; n++) {

	//[ this should be replace with shared memory above
	// calculate current b 
	int T_low = int(floor(n*dt));
	int T_high = T_low + 1;
	if (T_low >= b_length-1) {
	  T_high = T_low;
	}
	double b =  b0 + db_rand[T_low] + double(n*dt-T_low)*(db_rand[T_high]-db_rand[T_low]);	
	//]

	
	// perhaps k,l,m can be simplified to reduce the number of registers
	complex<double> k0 = dt*x1;
	complex<double> l0 = dt*x2;
	complex<double> m0 = dt*(-i*C*C*C*(4.0*QC*(b+i*d)-1.0)*x0 - 4.0*QC*C*C*x1 - i*C*(b + i*d)*x2);
	complex<double> k1 = dt*(x1+0.5*l0);
	complex<double> l1 = dt*(x2+0.5*m0);
	complex<double> m1 = dt*(-i*C*C*C*(4.0*QC*(b+i*d)-1.0)*(x0+0.5*k0) - 4.0*QC*C*C*(x1+0.5*l0) - i*C*(b + i*d)*(x2+0.5*m0));;
	complex<double> k2 = dt*(x1+0.5*l1);
	complex<double> l2 = dt*(x2+0.5*m1);
	complex<double> m2 = dt*(-i*C*C*C*(4.0*QC*(b+i*d)-1.0)*(x0+0.5*k1) - 4.0*QC*C*C*(x1+0.5*l1) - i*C*(b + i*d)*(x2+0.5*m1));;
	complex<double> k3 = dt*(x1+0.5*l2);
	complex<double> l3 = dt*(x2+0.5*m2);
	complex<double> m3 = dt*(-i*C*C*C*(4.0*QC*(b+i*d)-1.0)*(x0+0.5*k2) - 4.0*QC*C*C*(x1+0.5*l2) - i*C*(b + i*d)*(x2+0.5*m2));;
	x0 = x0 + 1.0/6.0*(k0 + k1+k1 + k2+k2 + k3);
	x1 = x1 + 1.0/6.0*(l0 + l1+l1 + l2+l2 + l3);
	x2 = x2 + 1.0/6.0*(m0 + m1+m1 + m2+m2 + m3);
	//]
	
	// calculate gain(x)
	double GAIN = std::pow(abs( 1.0/(x2 + 4.0*QC*C*C*x0)),2);
	
	// test for best gain
	if ( dGain[idx] < GAIN)
	  {
	    dLength[idx] = float(n*dt); // store current length if best gain
	    dGain[idx] = GAIN; // store Gain if best gain
	  }
	}
	// change gain to dB units
      	dGain[idx] =  10*std::log10(dGain[idx]);
    }
}


int main(int  argc , char*  argv[]  )
{
	int num_threads = 128; //default num_threads
	if (argc >= 2) {
	   num_threads = atoi(argv[1]);
	}
	if (argc >= 3) {
	   SIZE_C = atoi(argv[2]);
	   INCR_C = (MAX_C-MIN_C)/double(SIZE_C-1);
	}
  
  int num_loops = SIZE_b0*SIZE_C;
  // vectors for input data
  float *b0_in = (float*)malloc(sizeof(float)*SIZE_b0);
  float *C_in = (float*)malloc(sizeof(float)*SIZE_C);

  // create b0_in and C_in vectors
  b0_in[0] = MIN_b0;
  C_in[0] = MIN_C;
  for (int b = 1; b < SIZE_b0; b++) {
    b0_in[b] = b0_in[b-1] + INCR_b0;
  }
  for (int c = 1; c < SIZE_C; c++) {
    C_in[c] = C_in[c-1] + INCR_C;
  }

  // vectors for output data
  float *C_out = (float*)malloc(sizeof(float)*num_loops);
  float *b0_out = (float*)malloc(sizeof(float)*num_loops);
  float *b_rand = (float*)malloc(sizeof(float)*B_LENGTH);
  double *gain_dB_out = (double*)malloc(sizeof(double)*num_loops);
  float *length_out = (float*)malloc(sizeof(float)*num_loops);

  // create C_out and b0_out vectors
  int idx = 0;
  for (int b = 0; b<SIZE_b0; b++) {
    for (int c = 0; c<SIZE_C; c++) {
      b0_out[idx] = b0_in[b];
      C_out[idx] = C_in[c];
      idx++;
    }
  }

  // set up random b(x) distribution
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<float> distribution(0,b_sigma);
  for (int n = 0; n < B_LENGTH; n++) {
    b_rand[n] = distribution(generator);
  }

//  std::cout << SIZE_b0 << ' '  << ' ' << SIZE_C << '\n';
//    std::cout << b0_in[0] << ' ' << b0_in[SIZE_b0-1] << '\t' << C_in[0] << '\t' << C_in[SIZE_C -1] <<  '\n';
  
  // Allocate device memory 
  float *db0_in, *dC_in, *db_rand, *dLength;
  double *dGain;

  // START TIMER
  cudaEvent_t startTime, stopTime;
  float time;
  cudaEventCreate(&startTime);
  cudaEventCreate(&stopTime);
  cudaEventRecord(startTime,0);

  cudaMalloc((void**)&db0_in, sizeof(float)*num_loops);
  cudaMalloc((void**)&dC_in, sizeof(float)*num_loops);
  cudaMalloc((void**)&db_rand, sizeof(float)*B_LENGTH);
  cudaMalloc((void**)&dLength, sizeof(float)*num_loops);
  cudaMalloc((void**)&dGain, sizeof(double)*num_loops);

  // send all data to GPU 
  cudaMemcpy(db0_in, b0_out, sizeof(float)*num_loops, cudaMemcpyHostToDevice);
  cudaMemcpy(dC_in, C_out, sizeof(float)*num_loops, cudaMemcpyHostToDevice);
  cudaMemcpy(db_rand, b_rand, sizeof(float)*B_LENGTH, cudaMemcpyHostToDevice);
  cudaMemset(dGain, 0, sizeof(double)*num_loops);
  cudaMemset(dLength, 0, sizeof(float)*num_loops);

  // Call Cuda Kernels to solve ODE
  // set by argv[] //int num_threads = 512;
  int num_blocks = int(num_loops/num_threads)+1;

  ODE_Kernel<<<num_blocks,num_threads>>>(db0_in,dC_in,db_rand,dLength,dGain, num_loops, B_LENGTH);

  cudaMemcpy(length_out, dLength, sizeof(float)*num_loops, cudaMemcpyDeviceToHost);
  cudaMemcpy(gain_dB_out, dGain, sizeof(double)*num_loops, cudaMemcpyDeviceToHost);  
  
  // stop inclusive timer
  cudaEventRecord(stopTime,0);
  cudaEventSynchronize(stopTime);
  cudaEventElapsedTime(&time, startTime, stopTime);
  cudaEventDestroy(startTime);
  cudaEventDestroy(stopTime);

  // print output data to terminal
//  std::cout  << 'C' << '\t' << "b0" << '\t' << "maxLength" << '\t' << "maxGain" << '\t' << "maxGain[dB]" << std::endl;  
//  std::cout << "---------------------------------------------------" << std::endl;
  for(int n = 0; n < num_loops; n++) {
    std::cout << C_out[n] << '\t' << b0_out[n] << '\t' << length_out[n] << '\t' << gain_dB_out[n] << '\t'<< n << '\n';
  }   
  //std::cout << "time[ms]" << '\t' << "num_loops" << '\t'<< "ms/loop" << '\t' << "num_threads" << '\n';
  std::cout << time << '\t' << num_loops << '\t'<< time / float(num_loops) << '\t' << num_threads << '\n';

  // free device memory
  cudaFree(db0_in);
  cudaFree(dC_in);
  cudaFree(db_rand);
  cudaFree(dLength);
  cudaFree(dGain);

  // free host memory
  free(length_out);
  free(gain_dB_out);
  free(C_in);
  free(C_out);
  free(b0_in);
  free(b0_out);
  free(b_rand);
}
