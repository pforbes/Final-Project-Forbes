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

using namespace std;
using namespace boost::numeric::odeint;


//[ rhs_function
/* The type of container used to hold the state vector */
typedef vector< complex< double > > state_type;

/* set b0 values to increment: # points = (MAX-MIN)/INCR+1 */
const double MIN_b0 = 1.515;
const double MAX_b0 = 1.52;
const double INCR_b0 =0.00001;
const int SIZE_b0 = int((MAX_b0-MIN_b0)/INCR_b0)+1;

/* set C values to increment: # points = (MAX-MIN)/INCR+1 */
const double MIN_C = 0.019;
const double MAX_C = 0.020;
const double INCR_C= 0.0001;
const int SIZE_C = int((MAX_C-MIN_C)/INCR_C)+1;

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

  const vector < double > &m_b;
  const double mC;

public:
  ODE_SOLVER( vector <double>& b, double C ) : m_b(b) , mC(C) { }

  void operator() ( const state_type &x , state_type &dxdt , double  t  )
  {
    int T_low = int(floor(t));
    int T_high = T_low + 1;
    if (T_low >= B_LENGTH-1) {
       T_high = T_low;
    }
    double b_interp =  m_b[T_low] + double(t-T_low)*(m_b[T_high]-m_b[T_low]);
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


// struct write_state
// {
//     void operator()( const state_type &x ) const
//     {
//         std::cout << x[0] << "\t" << x[1] << "\n";
//     }
// };


int main(int  argc , char*  argv[]  )
{
  int num_procs = 1;

  stopwatch<std::milli, float> sw_incl;
  sw_incl.start();
  
  vector <double> b0_in(SIZE_b0);
  vector <double> C_in(SIZE_C);

  // create b0_in and C_in vectors
  b0_in[0] = MIN_b0;
  C_in[0] = MIN_C;
  for (int b = 1; b < SIZE_b0; b++) {
    b0_in[b] = b0_in[b-1] + INCR_b0;
  }
  for (int c = 1; c < SIZE_C; c++) {
    C_in[c] = C_in[c-1] + INCR_C;
  }

  int num_loops = SIZE_b0*SIZE_C;

  // vectors to store output
  vector <double> C_out(num_loops);
  vector <double> b0_out(num_loops);
  vector <double> gain_dB_out(num_loops);
  vector <double> length_out(num_loops);
	
  double b_sigma = 0.01; // normal distribution std deviation for b
  // set up random b(x) distribution
  unsigned seed = std::chrono::system_clock::now().time_since_epoch().count();
  std::default_random_engine generator (seed);
  std::normal_distribution<double> distribution(0,b_sigma);
  vector <double> b_rand;
  for (int n = 0; n < B_LENGTH; n++) {
    b_rand.push_back(distribution(generator));
  }


  // start timer 
  stopwatch<std::milli, float> sw;
  sw.start();

    
  /* iterate over b0 */ 
  for (int b0idx = 0; b0idx < b0_in.size(); b0idx++) {

    double b0 = b0_in[b0idx];
    // define b vector for variation along length
    vector< double> b (B_LENGTH,b0);
    std::transform(b.begin( ), b.end( ), b_rand.begin( ), b.begin( ),std::plus<double>( ));

    /* Iterate over C */
    for(int Cidx = 0; Cidx < C_in.size(); Cidx++) {
      double C = C_in[Cidx];
      int IDX = b0idx*C_in.size() + Cidx; // used to store results in x_out arrays

      state_type x(3);
      //[ state_initialization
      x[0] = 0.0;
      x[1] = 0.0;
      x[2] = 1.0;
      //]

      // //[ define_const_stepper
      runge_kutta4< state_type> stepper;
      
      // define observer vectors
      double maxGain = 0;  // max gain for b0,C pair
      double maxTime = -1;  // Length at which max gain occurs
	
      // initiate ODE solver class
      ODE_SOLVER myODE(b, C);

      // Solve ODE using class
      size_t steps = integrate_const(stepper,  myODE ,
					 x , 0.0 , L_max , L_inc,
					 max_observer(C, maxGain, maxTime) );

      /* output */
      //cout << times[n] << '\t' << gain[n] << '\t' << '\t' << std::endl;
      // cout << C << '\t' << b0 << '\t' << maxTime << '\t' << maxGain << '\t' << 10*log10(maxGain) << '\n';
      C_out[IDX] = C;
      b0_out[IDX] = b0;
      gain_dB_out[IDX] = 10*std::log10(maxGain);
      length_out[IDX] = maxTime;
    }
  }

  // stop timer
  sw.stop();
  float time = sw.count();

      // stop inclusive timer
  sw_incl.stop();
  float time_incl = sw_incl.count();

  // print output data to terminal
  // cout  << 'C' << '\t' << "b0" << '\t' << "maxLength" << '\t' << "maxGain" << '\t' << "maxGain[dB]" << std::endl;  
    //cout << "---------------------------------------------------" << std::endl;
  //for(int n = 0; n < num_loops; n++) {
    //  cout << C_out[n] << '\t' << b0_out[n] << '\t' << length_out[n] << '\t' << gain_dB_out[n] << '\n';
    // }  
  //  cout << "n_procs" << '\t' << "time[ms]" << '\t' << "num_loops" << '\t'<< "ms/loop" << '\t' << "time_incl[ms]" <<'\n';
  cout << num_procs << '\t' << time << '\t' << num_loops << '\t'<< time/float(num_loops)<< '\t' << time_incl << '\n';
    
}
