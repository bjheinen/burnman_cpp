/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#include "burnman/tools/averaging/averaging_utils.hpp"
#include "burnman/utils/constants.hpp"

namespace averaging{

namespace utils{

  // blank namespace for internals
  namespace {

    double hs_bulk_fn(
        const Eigen::ArrayXd& volumes,
        const Eigen::ArrayXd& bulk_moduli,
        const Eigen::ArrayXd& shear_moduli,
        bool n
    ) {
      Eigen::ArrayXd vol_frac = volumes / volumes.sum();
      double K_n = n ? bulk_moduli.maxCoeff() : bulk_moduli.minCoeff();
      double G_n = n ? shear_moduli.maxCoeff() : shear_moduli.minCoeff();
      double alpha_n = -3.0 / (3.0 * K_n + 4.0 * G_n);
      Eigen::Array<bool, Eigen::Dynamic, 1> mask = (bulk_moduli != K_n);
      Eigen::ArrayXd denominator = (1.0 / (bulk_moduli - K_n)) - alpha_n;
      Eigen::ArrayXd A_terms = vol_frac / denominator;
      A_terms = mask.select(A_terms, 0.0);
      double A_n = A_terms.sum();
      return K_n + A_n / (1.0 + alpha_n * A_n);
    }

    double hs_shear_fn(
        const Eigen::ArrayXd& volumes,
        const Eigen::ArrayXd& bulk_moduli,
        const Eigen::ArrayXd& shear_moduli,
        bool n
    ) {
      Eigen::ArrayXd vol_frac = volumes / volumes.sum();
      double K_n = n ? bulk_moduli.maxCoeff() : bulk_moduli.minCoeff();
      double G_n = n ? shear_moduli.maxCoeff() : shear_moduli.minCoeff();
      double beta_n = -3.0 * (K_n + 2.0 * G_n) / (5.0 * G_n * (3.0 * K_n + 4.0 * G_n));
      Eigen::Array<bool, Eigen::Dynamic, 1> mask = (shear_moduli != G_n);
      Eigen::ArrayXd denominator = (1.0 / (2.0 * (shear_moduli - G_n)) - beta_n);
      Eigen::ArrayXd B_terms = vol_frac / denominator;
      B_terms = mask.select(B_terms, 0);
      double B_n = B_terms.sum();
      return G_n + 0.5 * B_n / (1.0 + beta_n * B_n);
    }

  }

  double voigt_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X
  ) {
    Eigen::ArrayXd vol_frac = phase_volumes / phase_volumes.sum();
    return (vol_frac * X).sum();
  }

  double reuss_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X
  ) {
    Eigen::ArrayXd vol_frac = phase_volumes / phase_volumes.sum();
    double eps = constants::precision::abs_tolerance;
    // Warn when X <= 0 and |vol_frac| >= eps 
    Eigen::Array<bool, Eigen::Dynamic, 1> problematic =
        (X <= 0.0).select(vol_frac.abs() > eps, false);
    if (problematic.any()) {
        // Warn!
        // Could use cerr?
        // std::cerr << "Oops, called reuss_average with Xi <= 0!" << std::endl;
        return 0.0;
    } else {
      return 1.0 / (vol_frac / X).sum();
    }
  }

  double voigt_reuss_hill_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X
  ) {
    return (
        voigt_fn(phase_volumes, X)
        + reuss_fn(phase_volumes, X)
      ) / 2.0;
  }

  double lower_hs_bulk_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli
  ) {
    return hs_bulk_fn(volumes, bulk_moduli, shear_moduli, 0);
  }

  double lower_hs_shear_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli
  ) {
    return hs_shear_fn(volumes, bulk_moduli, shear_moduli, 0);
  }

  double upper_hs_bulk_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli
  ) {
    return hs_bulk_fn(volumes, bulk_moduli, shear_moduli, 1);
  }

  double upper_hs_shear_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli
  ) {
    return hs_shear_fn(volumes, bulk_moduli, shear_moduli, 1);
  }

} // namespace utils

} // namespace averaging
