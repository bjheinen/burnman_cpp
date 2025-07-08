/*
 * Copyright (c) 2025 Benedict Heinen
 *
 * This file is part of burnman_cpp and is licensed under the
 * GNU General Public License v3.0 or later. See the LICENSE file
 * or <https://www.gnu.org/licenses/> for details.
 *
 * burnman_cpp is based on BurnMan: <https://geodynamics.github.io/burnman/>
 */
#ifndef BURNMAN_CORE_COMPOSITE_HPP_INCLUDED
#define BURNMAN_CORE_COMPOSITE_HPP_INCLUDED

#include <vector>
#include <memory>
#include <optional>
#include <Eigen/Dense>
#include "burnman/core/material.hpp"
#include "burnman/core/mineral.hpp"
#include "burnman/core/solution.hpp"
#include "burnman/core/averaging_schemes.hpp"

/**
 * @class Composite
 * @brief Base class for a composite material.
 *
 * Base class for all composite materials.
 *
 * The static phases can be of any class derived from `Material', e.g.
 * Mineral, Solution, etc. As `Composite' is derived from `Material',
 * composite can be nested arbitrarily.
 *
 * The fractions of phases can be input as 'molar' or 'mass' and modified
 * using `set_fractions'.
 *
 */
class Composite : public Material {

 public:

  virtual ~Composite() = default;

  // Override of reset to include additional solution properties
  void reset() override;

  // Utility functions
  /**
   * @brief Utility function to map phase properties to Eigen::Array.
   *
   * Virtual dispatch should allow this to work for `Material' and any
   * derived classes.
   * 
   * Usage:
   *  map_phases_to_array(&Material::get_property);
   *    where, get_property() is the function needed
   *  alternatively, with a lambda function, e.g.
   * map_phases_to_array([](const Material& m) {
   *  return m.get_property() + 10.0;
   * });
   *
   * @note Properties expected as type double.
   *
   * @param func Pointer to public getter function in `Material'
   * @return properties Array of endmember properties.
   */
  template <typename Func>
  Eigen::ArrayXd map_phases_to_array(Func&& func) const {
    Eigen::ArrayXd mapped_properties(phases.size());
    std::transform(
      phases.begin(), phases.end(), mapped_properties.data(),
      [&func](const auto& phase) {
        return ((*phase).*func)();
    });
    return mapped_properties;
  }

  // Public getters for extra Composite functions
  /**
   * @brief Retrieves n_i * V_i.
   *
   * Uses a cached value if available, or calls
   * `Composite::compute_volume_fractions()` and caches the result.
   *
   * @note Use `Composite::reset()` to clear cached values.
   *
   * @return Molar fractions * volumes of phases in composite.
   */
  Eigen::ArrayXd get_volume_fractions() const;

  /**
   * @brief Retrieves the number of different elements in the composite.
   */
  int get_n_elements() const;

  /**
   * @brief Retrieves the total number of endmembers in the composite.
   */
  int get_n_endmembers() const;

  /**
   * @brief Retrieves a list of the number of endmembers in each phase.
   */
  std::vector<int> get_endmembers_per_phase() const;

  /**
   * @brief Retrieves a list of different elements in the composite.
   *
   * Elements are sorted into IUPAC order.
   */
  std::vector<std::string> get_elements() const;

  /**
   * @brief Returns a list of readable endmember names
   */
  std::vector<std::string> get_endmember_names() const;

  /**
   * Retrieves chemical formulae of all endmembers in the composite.
   */
  std::vector<FormulaMap> get_endmember_formulae() const;


 protected:

  // Overrides of defaults from Material
  double compute_molar_internal_energy() const override;
  double compute_molar_gibbs() const override;
  double compute_molar_helmholtz() const override;
  double compute_molar_mass() const override;
  double compute_molar_volume() const override;
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

  // Additional Composite compute functions
  Eigen::ArrayXd compute_volume_fractions() const;

 private:

  // Shared pointer to solution model class
  std::vector<std::shared_ptr<Material>> phases;

  // Unique pointer to averaging scheme
  std::unique_ptr<Averaging> averaging_scheme;

  // Molar fractions - define public getter if Composite subclassed.
  Eigen::ArrayXd molar_fractions;

  // Cached properties (can be reset)
  // 1D Eigen arrays (could be made Vectors)
  mutable std::optional<Eigen::ArrayXd> volume_fractions;

  // Stored (cached) properties not reset
  mutable std::optional<int> n_elements;
  mutable std::optional<int> n_endmembers;
  mutable std::optional<std::vector<int>> endmembers_per_phase;
  mutable std::optional<std::vector<std::string>> elements;
  mutable std::optional<std::vector<std::string>> endmember_names;
  mutable std::optional<std::vector<FormulaMap>> endmember_formulae;

  // Compute / setup functions for Composite properties
  void setup_endmember_properties() const;


};

#endif // BURNMAN_CORE_COMPOSITE_HPP_INCLUDED
