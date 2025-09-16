/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_EQUATION_OF_STATE_HPP_INCLUDED
#define BURNMAN_CORE_EQUATION_OF_STATE_HPP_INCLUDED

#include <string>
#include "burnman/utils/types/mineral_params.hpp"

namespace burnman {

/**
 * @class EquationOfState
 * @brief Base class for Equations of State.
 *
 * Provides the interface for an equation of state that a Mineral object
 * uses to determine its properties at a given \f$ P, T \f$.
 *
 * New equations of state should override the virtual functions
 * defined here. If called, the default implementations will
 * throw `exceptions::NotImplementedError`.
 *
 * All functions should accept and return values in SI units.
 *
 * In general, EquationOfState methods are functions of pressure,
 * temperature, and/or volume, as well as a `types::MineralParams' object.
 *
 */
class EquationOfState {

 public:
  // Virtual destructor for cleanup in derived classes
  virtual ~EquationOfState() = default;

  // Implemented EOS functions
  /**
   * @brief Computes the density of the mineral.
   *
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Density in [kg/m^3].
   */
  double compute_density(double volume, const types::MineralParams& params) const;

  // Function to override in derived classes with EOS specific implementations

  // Helper functions
  /**
   * @brief Validates the MineralParams object for EOS
   * TODO: documentation after implementation done
   *       function is not const
   *       MineralParams reference is not const (sometime need to set defaults)
   *       return type bool/void? depends on throwing error inside/outside func
   * @note Default implementation passes quietly
   *       Derived classes should override method.
   */
  virtual bool validate_parameters(types::MineralParams& params);

  // EOS relations
  /**
   * @brief Computes the molar volume of the mineral.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Molar volume in [m^3].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_volume(
    double pressure,
    double temperature,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the pressure (incl. cold and thermal).
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Pressure in [Pa].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_pressure(
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the Grueneisen parameter.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Grueneisen parameter [unitless].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_grueneisen_parameter(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes K_T.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Isothermal bulk modulus in [Pa].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_isothermal_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes K_S.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Isentropic bulk modulus in [Pa].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_isentropic_bulk_modulus_reuss(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the shear modulus, G.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Shear modulus in [Pa].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_shear_modulus(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes molar heat capacity at constant volume, C_v.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Heat capacity at constant volume in [J/K/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_molar_heat_capacity_v(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes molar heat capacity at constant pressure, C_p.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Heat capacity at constant pressure in [J/K/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_molar_heat_capacity_p(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the thermal expansivity.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Thermal expansivity in [1/K].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_thermal_expansivity(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the gibbs free energy.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Gibbs free energy in [J/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_gibbs_free_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;


  /**
   * @brief Computes the Helmholtz free energy.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Helmholtz free energy in [J/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_helmholtz_free_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the entropy.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Entropy in [J/K/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_entropy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the enthalpy.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Enthalpy in [J/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_enthalpy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

  /**
   * @brief Computes the molar internal energy of the mineral.
   *
   * @note Default implementation throws `exceptions::NotImplementedError'.
   *       Derived classes should override method.
   *
   * @param pressure The pressure to evaluate [Pa].
   * @param temperature The temperature to evaluate [K].
   * @param volume Molar volume of the mineral [m^3].
   * @param params Mineral parameters object of type types::MineralParams
   *
   * @return Internal energy in [J/mol].
   * @throws `exceptions::NotImplementedError' if default implementation called.
   */
  virtual double compute_molar_internal_energy(
    double pressure,
    double temperature,
    double volume,
    const types::MineralParams& params) const;

 private:

  // Private functions for throwing exceptions
  /**
   * @brief Helper function to throw `exceptions::NotImplementedError'
   *
   * Use __func__ for method name
   */
  [[noreturn]] void throw_not_implemented_error(const std::string& method) const;

  /**
   * @brief Helper function to get name of class (incl. derived)
   */
  std::string get_class_name() const;

};

} // namespace burnman

#endif // BURNMAN_CORE_EQUATION_OF_STATE_HPP_INCLUDED
