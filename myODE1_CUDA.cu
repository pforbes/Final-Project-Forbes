/*
 Copyright 2010-2012 Karsten Ahnert
 Copyright 2011-2013 Mario Mulansky
 Copyright 2013 Pascal Germroth
 Distributed under the Boost Software License, Version 1.0.
 (See accompanying file LICENSE_1_0.txt or
 copy at http://www.boost.org/LICENSE_1_0.txt)
 */


#include <iostream>
#include <vector>
#include <complex>
#include <boost/numeric/odeint.hpp>
#include <math.h>
#include <random>
#include <chrono>
#include "stopwatch.hpp"
#include <omp.h>
#include <cuda.h>

using namespace std;
using namespace boost::numeric::odeint;


//[ rhs_function
/* The type of container used to hold the state vector */
typedef vector< complex< double > > state_type;

	
const double b_sigma = 0.01; // normal distribution std deviation for b
/* set b0 values to increment: # points = (MAX-MIN)/INCR+1 */
const double MIN_b0 = 1.515;
const double MAX_b0 = 1.52;
const double INCR_b0 =0.00001;

/* set C values to increment: # points = (MAX-MIN)/INCR+1 */
const double MIN_C = 0.019;
const double MAX_C = 0.020;
const double INCR_C= 0.0001;

/* Set resolution in length along device ('time' for ODE solvers) */
const double L_max = 150.0;
const double L_inc = 0.001;

const int B_LENGTH = int(L_max) + 1;  // length of the b vector
const double d = 0;
const double QC = 1.85*0.25;

const complex<double> i(0,1);


//[ rhs_class
/* The rhs of x' = f(x) defined as a class */
class ODE_SOLVER{

  double& m_b;
  const int mb_size;
  const double mC;
  const double mb0;


public:
  ODE_SOLVER( double& b, int b_length, double C , double b0) : m_b(b), mb_size(b_length),  mC(C), mb0(b0) { }

  void operator()( const state_type &x , state_type &dxdt , double  t  )
  {
    int T_low = int(floor(t));
    int T_high = T_low + 1;
    if (T_low >= mb_size-1) {
       T_high = T_low;
    }
    double b_interp =  mb0 + m_b[T_low] + double(t-T_low)*(m_b[T_high]-m_b[T_low]);
    dxdt[0] = x[1];
    dxdt[1] = x[2];
    dxdt[2] = -i*mC*mC*mC*((4.0*QC)*(b_interp+i*d)-1.0)*x[0] - mC*mC*4.0*QC*x[1] -i*mC*(b_interp+i*d)*x[2];
  }
};
//]

struct max_observer
{
  double mC;
  double &mGain;
  double &mTime;
  
  max_observer(double C, double &maxGain, double &maxTime)
    :  mC( C ), mGain(maxGain), mTime( maxTime ) { }

  void operator()( const state_type &x , double t )
    {
      double GAIN = std::pow(abs( 1.0/(x[2] + 4.0*QC*mC*mC*x[0])),2);
      if ( mGain < GAIN)
	{
	  mTime = t;
	  mGain = GAIN;
	}
    }
};

__global__ void ODE_Kernel(double* db0_in, double* dC_in,double *db_rand, double* dLength,double* dGain, int num_loops, int b_length)
{

  int idx = blockIdx.x + threadIdx.x;
  if (idx < num_loops)
    {
      double b0 = db0_in[idx];
      double C  = dC_in[idx];
      
	//state initialization
      state_type x(3);
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 1.0;
      
      // //[ define_const_stepper
      runge_kutta4< state_type> stepper;
      
	// initiate ODE solver class
      ODE_SOLVER myODE(db_rand, b_length, C, b0);
      
      // Solve ODE using class
      size_t steps = integrate_const(stepper,  myODE ,
				     x , 0.0 , L_max , L_inc,
				     max_observer(C, dGain[idx], dLength[idx]) );
      
      dGain[idx] = 10*std::log10(dGain[idx]);
    }
}


// struct write_state
// {
//     void operator()( const state_type &x ) const
//     {
//         std::cout << x[0] << "\t" << x[1] << "\n";
//     }
// };


int main(int  argc , char*  argv[]  )
{
  int num_procs;
  if( argc == 2) {
    num_procs = atoi(argv[1]);
  } else {
    num_procs = omp_get_num_procs();
  }

  stopwatch<std::milli, float> sw_incl;
  sw_incl.start();
  
  vector <double> b0_in;
  for (double b0 = MIN_b0; b0 <= MAX_b0; b0 += INCR_b0) {
    b0_in.push_back(b0);
  }
  
  vector <double> C_in;
  for (double C = MIN_C; C <= MAX_C; C += INCR_C) {
    C_in.push_back(C);
  }
  
  vector <double> C_out;
  vector <double> b0_out;
  for (double b0 = MIN_b0; b0 <= MAX_b0; b0 += INCR_b0) {
    for (double C = MIN_C; C <= MAX_C; C += INCR_C) {
      C_out.push_back(C);
      b0_out.push_back(b0);
    }
  }
  
  int num_loops = C_in.size()*b0_in.size();

  // vectors to store output
  vector <double> gain_dB_out(num_loops);
  vector <double> length_out(num_loops);

  // set up random b(x) distribution
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(0,b_sigma);
  std::vector <double> b_rand;
  for (int n = 0; n < B_LENGTH; n++) {
    b_rand.push_back(distribution(generator));
  }

  
  /* Allocate device memory */
  double *db0_in, *dC_in, *db_rand, *dLength, *dGain;
  cudaMalloc((void**)&db0_in, sizeof(double)*num_loops);
  cudaMalloc((void**)&dC_in, sizeof(double)*num_loops);
  cudaMalloc((void**)&db_rand, sizeof(double)*b_rand.size());
  cudaMalloc((void**)&dLength, sizeof(double)*num_loops);
  cudaMalloc((void**)&dGain, sizeof(double)*num_loops);
  
  /* send all data to GPU */
  cudaMemcpy(&b0_in, db0_in, sizeof(double)*num_loops, cudaMemcpyDeviceToHost);
  cudaMemcpy(&C_in, dC_in, sizeof(double)*num_loops, cudaMemcpyDeviceToHost);
  cudaMemcpy(&b_rand, db_rand, sizeof(double)*b_rand.size(), cudaMemcpyDeviceToHost);
  cudaMemset(dGain, 0, sizeof(double)*num_loops);
  cudaMemset(dLength, -1, sizeof(double)*num_loops);

  /* Call Cuda Kernels to solve ODE */
  ODE_Kernel<<<16,16>>>(db0_in,dC_in,db_rand,dLength,dGain, num_loops, B_LENGTH);
  
  /* send results back to GPU */
  cudaMemcpy(&gain_dB_out, dGain, sizeof(double)*num_loops, cudaMemcpyDeviceToHost);
  cudaMemcpy(&length_out, dLength, sizeof(double)*num_loops, cudaMemcpyDeviceToHost);

  cudaFree(db0_in);
  cudaFree(dC_in);
  cudaFree(db_rand);
  cudaFree(dLength);
  cudaFree(dGain);
	   
  // stop inclusive timer
  sw_incl.stop();
  float time = sw_incl.count();

  // print output data to terminal
  cout  << 'C' << '\t' << "b0" << '\t' << "maxLength" << '\t' << "maxGain" << '\t' << "maxGain[dB]" << std::endl;  
  cout << "---------------------------------------------------" << std::endl;
  for(int n = 0; n < num_loops; n++) {
    cout << C_out[n] << '\t' << b0_out[n] << '\t' << length_out[n] << '\t' << gain_dB_out[n] << '\n';
  }   
  cout << "n_procs" << '\t' << "time[ms]" << '\t' << "num_loops" << '\t'<< "ms/loop" << '\n';
  cout << num_procs << '\t' << time << '\t' << num_loops << '\t'<< time/float(num_loops) << '\n';
    
}
