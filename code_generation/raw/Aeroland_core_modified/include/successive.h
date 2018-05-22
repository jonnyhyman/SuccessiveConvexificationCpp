#include <boost/numeric/odeint/external/eigen/eigen_algebra.hpp>
#include <boost/numeric/odeint.hpp>

using namespace boost::numeric::odeint;
runge_kutta4<DiscretizationODE::state_type, double, DiscretizationODE::state_type, double, vector_space_algebra> stepper;