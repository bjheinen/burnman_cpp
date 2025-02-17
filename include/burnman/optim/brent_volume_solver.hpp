/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_OPTIM_BRENT_VOLUME_SOLVER_HPP_INCLUDED
#define BURNMAN_OPTIM_BRENT_VOLUME_SOLVER_HPP_INCLUDED

#include "burnman/utils/eos.hpp"

namespace brent {
  /**
    * @brief Finds root of function in region [x_lo,x_hi] via Brent's method.
    *
    * Uses GSL implementation of Brent's method to a find a root of function
    * defined in a GSL function wrapper. Used here to find volume by
    * minimising P(V,T,...) - P.
    *
    * @param gsl_wrapper Pointer to a GSL function wrapper.
    * @param bm_params Parameter object for the GSL functor.
    * @param x_lo Low limit of starting interval.
    * @param x_hi High limit of starting interval.
    *
    * @return Thermal energy in [J/mol].
    */
  double find_root(
    double (*gsl_wrapper)(double, void*),
    const ParamsGSL::SolverParams_P bm_params,
    double x_lo,
    double x_hi);
}

#endif // BURNMAN_OPTIM_BRENT_VOLUME_SOLVER_HPP_INCLUDED
