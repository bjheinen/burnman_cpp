/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_MINERAL_HPP_INCLUDED
#define BURNMAN_CORE_MINERAL_HPP_INCLUDED

#include <memory>
#include <string>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/utils/types/mineral_params.hpp"
#include "burnman/eos/components/excess_params.hpp"
#include "burnman/core/equation_of_state.hpp"
#include "burnman/core/material.hpp"

namespace burnman {

/**
 * @class Mineral
 * @brief Base class for minerals (endmembers).
 *
 * This class provides an interface for computing mineral properties.
 * Mineral objects define the endmembers in Solutions, or act as stand-alone
 * phases. Solid solutions should use `Solution' not `Mineral'.
 *
 * States and properties of the mineral can only be queried after
 * setting the P, T or V, T condition with set_state.
 * Use `set_method` to set the EOS for computing properties.
 * @see types::EOSType for predefined equations of state, or pass a pointer to a
 * class derived from EquationOfState.
 *
 * All material parameters expected in SI units.
 * Elastic moduli in Pa NOT GPa, Debye temperature etc. in K not C.
 * Additionally, the reference volume should be in m^3/(mol molecule) and not
 * the unit cell volume in Angstrom^3.
 * To convert: V_uc * 1e-30 * N_a / Z, where N_a is Avogadro's number and Z is
 * the number of formula units per unit cell.
 *
 */
class Mineral : public Material{
 public:

  /**
   * @brief Parameters for the mineral.
   */
  types::MineralParams params;

  /**
   * @brief Shared pointer to the EquationOfState used to compute properties.
   */
  std::shared_ptr<EquationOfState> eos_method;

  /**
   * @brief Sets the parameters for thermodynamic property modifiers.
   *
   * The property modifiers are used to calculate excess thermodynamic properties based
   * on a model, and summed to give total modification excesses.
   *
   * Possible models are:
   *   - second order transitions (landau, landau_slb_2022, landau_hp)
   *   - order-disorder (bragg_williams)
   *   - magnetism (magnetic_chs)
   *   - linear
   * @see burnman::eos::excesses
   *
   * Usage:
   * @code{.cpp}
   * using namespace burnman::eos::excesses;
   * ExcessParamVector p = {
   *   LinearParams{delta_V, delta_S, delta_E},
   *   LandauParams{Tc_0, V_D, S_D},
   *   DebyeParams{Cv_inf, Theta_0}
   * };
   * mineral.set_property_modifier_params(p);
   * @endcode
   *
   * @param excess_params List of excess parameter structs.
   * @note Calling this function will reset any previously set property modifiers and recompute where possible.
   */
  void set_property_modifier_params(eos::excesses::ExcessParamVector excess_params);

  /**
   * @brief Retrieves the current property modifiers.
   *
   * Returns the summed excesses from all property modifiers currently set.
   * @see `compute_property_modifiers()` for details on how these are calculated.
   *
   * @return Excesses struct containing the summed excess properties.
   */
  eos::excesses::Excesses get_property_modifiers() const;

  // Override public methods
  std::string get_name() const override;
  void set_state(double new_pressure, double new_temperature) override;
  void set_method(types::EOSType new_method) override;
  void set_method(std::shared_ptr<EquationOfState> new_method) override;

 protected:

  /**
   * @brief Computes and sums excesses from all property modifiers.
   *
   * Uses vector or property modifiers defined in property_modifier_params.
   * Stores values in the Excesses struct property_modifier_excesses.
   *
   * To calculate thermodynamic properties use the following functions,
   * where _o suffix implies the original value,
   * and e is property_modifier_excesses
   *
   * Gibbs = Gibbs_o + e.G
   * S = S_o - e.dGdT
   * V = V_o + e.dGdP
   * K_T = V / ((V_o / K_T_o) - e.d2GdP2)
   * C_p = C_p_o - T * e.d2GdT2
   * alpha = ((alpha_o * V_o) + e.d2GdPdT) / V
   *
   * @note Modifies property_modifier_excesses.
   */
  void compute_property_modifiers();

  // Overridden compute functions
  double compute_molar_internal_energy() const override;
  double compute_molar_gibbs() const override;
  double compute_molar_helmholtz() const override;
  double compute_molar_mass() const override;
  double compute_molar_volume() const override;
  double compute_molar_volume_unmodified() const override;
  double compute_density() const override;
  double compute_molar_entropy() const override;
  double compute_molar_enthalpy() const override;
  double compute_isothermal_bulk_modulus_reuss() const override;
  double compute_isentropic_bulk_modulus_reuss() const override;
  double compute_isothermal_compressibility_reuss() const override;
  double compute_isentropic_compressibility_reuss() const override;
  double compute_shear_modulus() const override;
  double compute_p_wave_velocity() const override;
  double compute_bulk_sound_velocity() const override;
  double compute_shear_wave_velocity() const override;
  double compute_grueneisen_parameter() const override;
  double compute_thermal_expansivity() const override;
  double compute_molar_heat_capacity_v() const override;
  double compute_molar_heat_capacity_p() const override;
  double compute_isentropic_thermal_gradient() const override;
  types::FormulaMap compute_formula() const override;

 private:
  // Excesses struct for property modifiers
  // Values initialise to 0
  eos::excesses::Excesses property_modifier_excesses;
  eos::excesses::ExcessParamVector property_modifier_params;

};

} // namespace burnman

#endif // BURNMAN_CORE_MINERAL_HPP_INCLUDED
