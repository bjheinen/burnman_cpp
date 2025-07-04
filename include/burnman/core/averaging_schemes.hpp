/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED
#define BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED

#include <Eigen/Dense>

/**
 * @class Averaging
 * @brief Base interface for averaging shemes.
 *
 * Functions return scalar values for a list of properties and
 * volume fractions.

 * Derived classes should implement averaging of elastic moduli.
 */
class Averaging {

 public:

  virtual ~Averaging() = default;

  // Moduli averages need to be implemented in subclasses.
  /**
   * @brief Average bulk modulus, `K', for a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param bulk_moduli Bulk moduli of each phase in [Pa].
   * @param shear_moduli Shear moduli of each phase in [Pa].
   * @return Average bulk modulus, K, in [Pa].
   */
  virtual double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const = 0;

  /**
   * @brief Average shear modulus, `G', for a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param bulk_moduli Bulk moduli of each phase in [Pa].
   * @param shear_moduli Shear moduli of each phase in [Pa].
   * @return Average shear modulus, G, in [Pa].
   */
  virtual double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const = 0;

  // Common default implementations
  /**
   * @brief Density average of a composite.
   *
   * Average assuming: \f$ \rho = \frac{\sum_i \rho_i V_i}{\sum_i V_i} \f$.
   *
   * @param volumes Volumes of each phase in [m^3].
   * @param densities Densities of each phase in [kg/m^3].
   * @return \f$ \rho \f$ in [kg/m^3].
    */
  virtual double average_density(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& densities) const;

  /**
   * @brief Average thermal expansion coefficient of a composite.
   * @param volumes Volumes of each phase in [m^3].
   * @param alphas Thermal expansivity of each phase in [1/K].
   * @return \f$ \alpha \f$ in [1/K].
    */
  virtual double average_thermal_expansivity(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& alphas) const;

  /**
   * @brief Average C_v of a composite.
   *
   * Averages C_v by molar fractions as in Ita, et al., 1992 [eq. 16].

   * @param fractions Molar fractions of each phase.
   * @param c_v Heat capacities at constant volume of each phase in [J/K/mol].
   * @return Average C_v in [J/K/mol].
    */
  virtual double average_heat_capacity_v(
    const Eigen::ArrayXd& fractions,
    const Eigen::ArrayXd& c_v) const;

  /**
   * @brief Average C_p of a composite.
   * @param fractions Molar fractions of each phase.
   * @param c_p Heat capacities at constant volume of each phase in [J/K/mol].
   * @return Average C_p in [J/K/mol].
    */
  virtual double average_heat_capacity_p(
    const Eigen::ArrayXd& fractions,
    const Eigen::ArrayXd& c_p) const;

  // Static utility functions for V & R averages
  /**
   * @brief Calculated Voigt (iso-strain) average.
   *
   * @param phase_volumes Partial volumes or molar fractions
   * @param X Phase moduli to average
   *
   * @return Voigt average
   */
  static double voigt_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X);

  /**
   * @brief Calculated Reuss (iso-stress) average.
   *
   * @param phase_volumes Partial volumes or molar fractions
   * @param X Phase moduli to average
   *
   * @return Reuss average
   */
  static double reuss_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X);

  /**
   * @brief Calculated Voigt-Reuss-Hill average.
   *
   * VRH average is the arithmetic mean of the V & R bounds.
   *
   * @param phase_volumes Partial volumes or molar fractions
   * @param X Phase moduli to average
   *
   * @return VRH average
   */
  static double voigt_reuss_hill_fn(
    const Eigen::ArrayXd& phase_volumes,
    const Eigen::ArrayXd& X);

 protected:

  /**
   * @see `HashinShtrikmanLower::average_bulk_moduli'.
   */
  static double lower_hs_bulk_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli);

  /**
   * @see `HashinShtrikmanLower::average_shear_moduli'.
   */
  static double lower_hs_shear_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli);

  /**
   * @see `HashinShtrikmanUpper::average_bulk_moduli'.
   */
  static double upper_hs_bulk_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli);

  /**
   * @see `HashinShtrikmanUpper::average_shear_moduli'.
   */
  static double upper_hs_shear_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli);

 private:

  /**
   * @brief Compute HS bulk modulus bound.
   */
  static double hs_bulk_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli,
    bool n);

  /**
   * @brief Compute HS shear modulus bound.
   */
  static double hs_shear_fn(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli,
    bool n);
};


/**
 * @class Voigt
 * @brief Voigt averaging scheme.
 *
 * Class for computing the Voigt (iso-strain) bound for elastic
 * properties.
 *
 * Overrides bulk and shear moduli averaging in `Averaging'.
 */
class Voigt : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};


/**
 * @class Reuss
 * @brief Reuss averaging scheme.
 *
 * Class for computing the Reuss (iso-stress) bound for elastic
 * properties.
 *
 * Overrides bulk and shear moduli averaging in `Averaging'.
 */
class Reuss : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};


/**
 * @class VoigtReussHills
 * @brief VRH averaging scheme.
 *
 * Class for computing the Voigt-Reuss-Hill average for elastic properties.
 * Defined as arithmetic mean of Voigt and Reuss bounds.
 * Overrides bulk and shear moduli averaging in `Averaging'.
 */
class VoigtReussHill : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};


/**
 * @class HashinShtrikmanLower
 * @brief Lower Hashin-Shtrikman bound.
 *
 * Class for computing the lower Hashin-Shtrikman bound for
 * elastic properties.

 * Overrides bulk and shear moduli averaging in `Averaging'.
 * Uses formulae from Watt, 1976.
 *
 * @note The Hashin-Shtrikman bounds are tighter than the Voigt and Reuss
 * bounds because they make the additional assumption that the orientation
 * of the phases are statistically isotropic. In some cases this may be a
 * good assumption, and in others it may not be.
 */
class HashinShtrikmanLower : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};


/**
 * @class HashinShtrikmanUpper
 * @brief Upper Hashin-Shtrikman bound.
 *
 * Class for computing the upper Hashin-Shtrikman bound for
 * elastic properties.

 * Overrides bulk and shear moduli averaging in `Averaging'.
 * Uses formulae from Watt, 1976.
 *
 * @note The Hashin-Shtrikman bounds are tighter than the Voigt and Reuss
 * bounds because they make the additional assumption that the orientation
 * of the phases are statistically isotropic. In some cases this may be a
 * good assumption, and in others it may not be.
 */
class HashinShtrikmanUpper : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};


/**
 * @class HashinShtrikman
 * @brief Hashin-Shtrikman averaging scheme.
 *
 * Class for computing the Hashin-Shtrikman average for elastic properties.
 * Defined as arithmetic mean of upper and lower bounds.
 * Overrides bulk and shear moduli averaging in `Averaging'.
 */
class HashinShtrikman : public Averaging {
 public:
  double average_bulk_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;

  double average_shear_moduli(
    const Eigen::ArrayXd& volumes,
    const Eigen::ArrayXd& bulk_moduli,
    const Eigen::ArrayXd& shear_moduli) const override;
};

#endif // BURNMAN_CORE_AVERAGING_SCHEMES_HPP_INCLUDED
