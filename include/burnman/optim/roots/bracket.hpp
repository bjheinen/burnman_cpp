/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_OPTIM_ROOTS_BRACKET_HPP_INCLUDED
#define BURNMAN_OPTIM_ROOTS_BRACKET_HPP_INCLUDED

#include <cassert>
#include <cmath>

namespace optim {
namespace roots {
  /**
    * @brief Expands a bracketing interval [x_lo,x_hi] until valid.
    *
    * Used to check if a bracketing interval used as an input to Brent's
    * method is valid (i.e. a sign change exists between f(x_lo), f_x_hi).
    * For invalid intervals the bracket is expanded downhill by
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
    * @param max_iter Maximum number of iterations.
    *
    * @return true if valid root found in interval, false if not.
    */
  template <typename SolverParamsType>
  bool bracket_root(
    double (*gsl_wrapper)(double, void*),
    SolverParamsType gsl_wrapper_params,
    double& x_lo,
    double& x_hi,
    double dx,
    double expansion_factor = 1.0,
    int max_iter = 100)
  {
    assert(dx > 0 && "dx must be +ve!");
    // Evaluate function at interval bounds
    double f_lo = gsl_wrapper(x_lo, &gsl_wrapper_params);
    double f_hi = gsl_wrapper(x_hi, &gsl_wrapper_params);
    // Expand towards f(x)=0 if no sign change
    int n_iter = 0;
    bool bracket_found = (f_lo * f_hi < 0);
    while (!bracket_found && n_iter < max_iter) {
      if (std::abs(f_lo) < std::abs(f_hi)) {
        x_lo -= dx;
        f_lo = gsl_wrapper(x_lo, &gsl_wrapper_params);
      } else {
        x_hi += dx;
        f_hi = gsl_wrapper(x_hi, &gsl_wrapper_params);
      }
      bracket_found = (f_lo * f_hi < 0);
      // TODO Too risky to use factor? (overshoot leads to error at invalid V)
      dx *= expansion_factor;
      n_iter++;
    }
    // Return true if valid bracket
    return bracket_found;
  }
} // namespace roots
} // namespace optim

#endif // BURNMAN_OPTIM_ROOTS_BRACKET_HPP_INCLUDED
