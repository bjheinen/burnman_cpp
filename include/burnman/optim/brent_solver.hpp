/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_BRENT_SOLVER_HPP_INCLUDED
#define BURNMAN_OPTIM_BRENT_SOLVER_HPP_INCLUDED

#include <cassert>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_roots.h>
#include "burnman/utils/eos.hpp"
#include "burnman/utils/constants.hpp"

#include <iostream>

namespace brent {
  /**
    * @brief Finds root of function in region [x_lo,x_hi] via Brent's method.
    *
    * Uses GSL implementation of Brent's method to a find a root of function
    * defined in a GSL function wrapper. 
    * The main use here to find volume by minimising P(V,T,...) - P
    *
    * @param gsl_wrapper Pointer to a GSL function wrapper.
    * @param bm_params Parameter object for the GSL functor.
    * @param x_lo Low limit of starting interval.
    * @param x_hi High limit of starting interval.
    *
    * @return root x at f(x) = 0
    */
  template <typename SolverParamsType>
  double find_root(
    double (*gsl_wrapper)(double, void*),
    SolverParamsType gsl_wrapper_params,
    double x_lo,
    double x_hi) 
  {
    // Find root in [P(V) - P] to find V that fits P
    // Set up GSL function object and params
    gsl_function obj_func;
    obj_func.function = gsl_wrapper;
    obj_func.params = &gsl_wrapper_params;

    // Define and allocate the solver
    const gsl_root_fsolver_type* T = gsl_root_fsolver_brent;
    gsl_root_fsolver* solver = gsl_root_fsolver_alloc(T);
    // Set the GSL solver
    gsl_root_fsolver_set(solver, &obj_func, x_lo, x_hi);
    // Iterate the solver
    int iter = 0;
    int maxiter = 100;
    int status;
    double root;
    double atol = constants::precision::abs_tolerance;
    double rtol = constants::precision::rel_tolerance_eps;
    do
      {
        iter++;
        status = gsl_root_fsolver_iterate(solver);
        root = gsl_root_fsolver_root(solver);
        x_lo = gsl_root_fsolver_x_lower(solver);
        x_hi = gsl_root_fsolver_x_upper(solver);
        status = gsl_root_test_interval(x_lo, x_hi,
                                        atol, rtol);
        // Leaving in print from GSL example to debug
        // Debug print statements
        //if (status == GSL_SUCCESS)
        //  printf ("Converged:\n");
        //printf ("%5d [%.7f, %.7f] %.7f %.7f\n",
        //        iter, x_lo, x_hi,
        //        root,
        //        x_hi - x_lo);
      } while (status == GSL_CONTINUE && iter < maxiter);

    // Get final root (can delet root getting above when tested)
    root = gsl_root_fsolver_root(solver);
    // Clean up
    gsl_root_fsolver_free(solver);

    if (status == GSL_SUCCESS) {
        return root;
    } else {
        return -1; // throw an error? "Error: Root-finding did not converge!"
    }
  }

  /**
    * @brief Expands a bracketing interval [x_lo,x_hi] until valid.
    *
    * Used to check if a bracketing interval used as an input to Brent's
    * method is valid (i.e. a sign change exists between f(x_lo), f_x_hi).
    * For invalid intervals the bracket is expanded in the 0 direction by
    * a small amount (dx). The step increases with each iteration by a
    * contant factor (expansion_factor).
    * @see brent::find_root
    *
    * @param gsl_wrapper Pointer to a GSL style function wrapper.
    * @param bm_params Parameter object for the GSL functor.
    * @param[in,out] x_lo Low limit of starting interval.
    * @param[in,out] x_hi High limit of starting interval.
    * @param dx Initial step size for bracket expansion.
    * @param expansion_factor Factor to increase step size each iteration.
    * @param max_iter Maximum number of iterations
    *
    * @return success Flag.
    */
  template <typename SolverParamsType>
  bool bracket_root(
    double (*gsl_wrapper)(double, void*),
    SolverParamsType gsl_wrapper_params,
    double& x_lo,
    double& x_hi,
    double dx,
    double expansion_factor = 1.2,
    int max_iter = 100)
  {
    // TODO: Use higher max_iter?
    assert(dx > 0 && "dx must be +ve!");
    // Evaluate function at interval bounds
    double f_lo = gsl_wrapper(x_lo, &gsl_wrapper_params);
    double f_hi = gsl_wrapper(x_hi, &gsl_wrapper_params);
    // Expand towards f(x)=0 if no sign change
    int n_iter = 0;
    while (f_lo * f_hi >= 0 && n_iter < max_iter) {
      if (std::abs(f_lo) < std::abs(f_hi)) {
        x_lo -= dx;
        f_lo = gsl_wrapper(x_lo, &gsl_wrapper_params);
      } else {
        x_hi += dx;
        f_hi = gsl_wrapper(x_hi, &gsl_wrapper_params);
      }
      // TODO Too risky to use factor? (overshoot leads to error at invalid V)
      //dx *= expansion_factor;
      n_iter++;
    }
    // Return true if valid bracket
    return (f_lo * f_hi < 0);
  }
}
#endif // BURNMAN_OPTIM_BRENT_SOLVER_HPP_INCLUDED
