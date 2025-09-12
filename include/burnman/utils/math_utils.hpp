/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_UTILS_MATH_UTILS_INCLUDED
#define BURNMAN_UTILS_MATH_UTILS_INCLUDED

#include <cmath>
#include <Eigen/Dense>
#include "burnman/utils/constants.hpp"

namespace burnman {
namespace utils {

  /**
   * @brief Safe inverse to prevent infinities at x=0
   *
   * Returns 1/x, with a 1st order expansion of 1/x about eps:
   *   2/eps - x/eps/eps
   * Work for 1D & 2D arrays.
   * @see `constants::precision' for eps value.
   *
   * @param x Eigen array to calculate inverse
   * @return Inverse(-ish) of x
   */
  template <typename Derived>
  typename Derived::PlainObject inverseish(const Eigen::ArrayBase<Derived>& x) {
    using ArrayType = typename Derived::PlainObject;
    // Grab eps from constants
    double eps = constants::precision::inverseish_eps;
    ArrayType inverse_x = (2.0 / eps - x / (eps * eps)); // To work with float and double could add .template cast<typename Derived::Scalar>();
    inverse_x = (x > eps).select(x.inverse(), inverse_x);
    return inverse_x;
  }

  /**
   * @brief Safe log to prevent infinities at x=0
   * 
   * Return log(x), with 2nd order expansion of log(x) about eps:
   *  log(eps) - sum_k=1^infty (f_eps)^k / k
   * Works for 1D & 2D arrays.
   * @see `constants::precision' for eps value.
   * 
   * @param x Eigen array to calculate log.
   * @return log-ish(x).
   */
  template <typename Derived>
  typename Derived::PlainObject logish(const Eigen::ArrayBase<Derived>& x) {
    using ArrayType = typename Derived::PlainObject;
    double eps = constants::precision::logish_eps;
    ArrayType f_eps = 1.0 - x/eps;  
    ArrayType log_x = (std::log(eps) - f_eps - f_eps.square() / 2.0);
    log_x = (x > eps).select(x.log(), log_x);
    return log_x;
  }

} // namespace utils
} // namespace burnman

#endif // BURNMAN_UTILS_MATH_UTILS_INCLUDED
