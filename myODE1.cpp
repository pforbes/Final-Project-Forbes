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


using namespace std;
using namespace boost::numeric::odeint;

//[ rhs_function
/* The type of container used to hold the state vector */
typedef vector< complex< double > > state_type;

const double L_max = 150;
const double L_inc = 0.001;
const double C = 0.0211615;
const double d = 0;
const double QC = 1.85*0.25;
const double b = 1.51766;
const complex <double> i (0,1);


/* The rhs of x' = f(x) */
void harmonic_oscillator( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = x[2];
    dxdt[2] = -i*C*C*C*((4*QC)*(b+i*d)-1.0)*x[0] - C*C*4*QC*x[1] -i*C*(b+i*d)*x[2];
}
//]

/* The rhs of x' = f(x) */
void FW( const state_type &x , state_type &dxdt , const double /* t */ )
{
    dxdt[0] = x[1];
    dxdt[1] = x[2];
    dxdt[2] = (-i*C*C*C*(4*QC)*(b-i*d)+1.0)*x[0] - C*C*4*QC*x[1] -i*C*(b-i*d)*x[2];
}
//]





// //[ rhs_class
// /* The rhs of x' = f(x) defined as a class */
// class harm_osc {

//     double m_gam;

// public:
//     harm_osc( double gam ) : m_gam(gam) { }

//     void operator() ( const state_type &x , state_type &dxdt , const double /* t */ )
//     {
//         dxdt[0] = x[1];
//         dxdt[1] = -x[0] - m_gam*x[1];
//     }
// };
// //]





//[ integrate_observer
struct push_back_state_and_time
{

    std::vector< state_type >& m_states;
    std::vector< double >& m_times;
    std::vector< double>& m_gain;

  push_back_state_and_time( std::vector< state_type > &states , std::vector< double > &times , std::vector< double > & gain)
    : m_states( states ) , m_times( times ) , m_gain( gain ){ }

    void operator()( const state_type &x , double t )
    {
        m_states.push_back( x );
        m_times.push_back( t );
	//m_gainFW.push_back(std::pow(abs((x[2] + 4*QC*C*C*x[0])),2));
       	m_gain.push_back(std::pow(abs( 1.0/(x[2] + 4*QC*C*C*x[0])),2));
    }
};
//]

struct write_state
{
    void operator()( const state_type &x ) const
    {
        std::cout << x[0] << "\t" << x[1] << "\n";
    }
};


int main(int /* argc */ , char** /* argv */ )
{

    //[ state_initialization
    state_type x(3);
    x[0] = 0.0;
    x[1] = 0.0;
    x[2] = 1.0;
    //]



    //[ integration
    //   size_t steps = integrate( harmonic_oscillator ,
    //        x , 0.0 , 10.0 , 0.1 );
    //]



    //[ integration_class
 //   harm_osc ho(0.15);
 //   steps = integrate( ho ,
 //           x , 0.0 , 10.0 , 0.1 );
    //]


    //[ integrate_observ
    vector<state_type> x_vec;
    vector<double> times;
    vector<double> gain;
    
//    steps = integrate( harmonic_oscillator ,
//            x , 0.0 , 10.0 , 0.1 ,
//            push_back_state_and_time( x_vec , times ) );


    //[ define_const_stepper
    runge_kutta4< state_type> stepper;
    size_t steps = integrate_const( stepper , harmonic_oscillator , x , 0.0 , L_max , L_inc, push_back_state_and_time(x_vec, times, gain) );
    //]

    /* output */
    for( size_t n=0; n<=steps; n++ )
    {
      //cout << times[n] << '\t' << gain[n] << '\t' << '\t' << std::endl;
      cout << times[n] << '\t' << x_vec[n][0] << '\t' << x_vec[n][1] << '\t' << x_vec[n][2] << '\t' << gain[n] << '\t' << 10*log10(gain[n]) << '\n';
    }
    //]

}
