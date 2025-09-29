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

#include <cstddef>
#include <algorithm>
#include <initializer_list>
#include <iterator>
#include <memory>
#include <optional>
#include <type_traits>
#include <utility>
#include <vector>
#include <Eigen/Dense>
#include "burnman/utils/types/simple_types.hpp"
#include "burnman/tools/averaging/averaging_schemes.hpp"
#include "burnman/core/material.hpp"
#include "burnman/core/composite_material.hpp"

namespace burnman {

/**
 * @class Assemblage
 * @brief Base class for a composite material (assemblage).
 *
 * Base class for all assemblages.
 *
 * The static phases can be of any class derived from `Material', e.g.
 * Mineral, Solution, etc. As `Assemblage' is derived from `Material',
 * Assemblage can be nested arbitrarily.
 *
 * The fractions of phases can be input as 'molar' or 'mass' and modified
 * using `set_fractions'.
 *
 * See `CompositeMaterial' for stoichiometry specific properties.
 *
 */
class Assemblage : public CompositeMaterial {

 public:

  virtual ~Assemblage() = default;

  // Override of reset_cache to include additional solution properties
  void reset_cache() override;

  // Override of clear_computed_properties() to include additional assemblage properties
  void clear_computed_properties() override;

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

  /**
   * @brief Utility function to map endmember properties to std::vector.
   *
   * Usage: see `Assemblage::map_endmembers_to_array'.
   *
   * @param func Pointer to public getter function in `Mineral'
   * @return vector of endmember properties (of type T)
   */
  template <typename T, typename Func>
  std::vector<T> map_phases_to_vector(Func&& func) const {
    std::vector<T> mapped_properties;
    mapped_properties.reserve(phases.size());
    std::transform(
      phases.begin(), phases.end(), std::back_inserter(mapped_properties),
      [&func](const auto& phase) {
        return ((*phase).*func)();
      });
    return mapped_properties;
  }

  // Public setters for assemblage properties

  // Set phases

  /**
   * @brief Appends a single phase to the phase list
   *
   * Overload to add a phase by value (no external access to pointer).
   * Consumes phase object
   */
  template <typename PhaseType>
  void add_phase(PhaseType&& phase) {
    static_assert(
      std::is_base_of_v<Material, std::remove_reference_t<PhaseType>>,
      "PhaseType must derive from Material!"
    );
    this->phases.push_back(
      std::make_shared<std::remove_reference_t<PhaseType>>(std::move(phase))
    );
  }

  /**
   * @brief Appends a single phase to the phase list
   *
   * Overload to add a phase pointer (for external access to pointer).
   * Do:
   *   auto ol_ptr = std::make_shared<Solution>(olivine);
   *   assemblage.add_phase(ol_ptr);
   */
  void add_phase(const std::shared_ptr<Material>& phase_ptr);

  // Multiple phases

  /**
   * @brief Sets the phase list
   *
   * Overload to add phases by value (no external access to pointers).
   * This overload works with variadic template types so you can do:
   *  add_phases(bdg, fper, capv, ...);
   * Storing objects in a container first isn't supported as we want
   * to allow mixed types including Assemblage (so run into recursion issues).
   * Consumes phase objects.
   */
  template<typename... Ts>
  void add_phases(Ts&&... args) {
      (add_phase(std::forward<Ts>(args)), ...);
  }

  /**
   * @brief Sets the phase list
   *
   * Overload to add phase pointers (for external access to pointers).
   * Construct a std::vector<std::shared_ptr<Material>> first.
   */
  void add_phases(
    const std::vector<std::shared_ptr<Material>>& phase_ptr_list);

  /**
   * @brief Sets the phase list
   *
   * Overload to add phase pointers (for external access to pointers).
   * Overload to pass an initialiser list.
   */
  void add_phases(
    std::initializer_list<std::shared_ptr<Material>> phase_ptr_list);

  /**
   * @brief Sets the fraction of each phase in the assemblage.
   *
   * Fractions can be types::FractionType::Molar or types::FractionType::Mass.
   * Mass fractions will be converted to molar fractions automatically.
   */
  void set_fractions(
    const Eigen::ArrayXd& fractions,
    const types::FractionType fraction_type = types::FractionType::Molar);


  /**
   * @brief Sets the fraction of each phase in the assemblage.
   *
   * Fractions can be types::FractionType::Molar or types::FractionType::Mass.
   * Mass fractions will be converted to molar fractions automatically.
   *
   * Overload to pass an initialiser list: set_fractions({0.5, 0.5});
   */
  void set_fractions(
    std::initializer_list<double> fractions,
    const types::FractionType fraction_type = types::FractionType::Molar);

  /**
   * @brief Sets the averaging scheme to use for computing properties.
   *
   * Sets the averaging scheme as an types::AveragingType, can be one of:
   *  Voigt
   *  Reuss
   *  VRH
   *  HashinShtrikmanLower
   *  HashinShtrikmanUpper
   *  HashinShtrikman
   *
   */
  void set_averaging_scheme(types::AveragingType scheme_type);

  /**
   * @brief Sets the averaging scheme to use for computing properties.
   *
   * Sets the averaging scheme to a custom class. Must be derived from `AveragingScheme'.
   *
   */
  void set_averaging_scheme(std::shared_ptr<averaging::AveragingScheme> custom_scheme);

  // Override public methods
  void set_state(double new_pressure, double new_temperature) override;
  void set_method(types::EOSType new_method) override;
  void set_method(std::shared_ptr<EquationOfState> new_method) override;

  /**
   * @brief Set n_moles (Used to convert mole fractions / absolute phase amounts).
   */
  void set_n_moles(double new_n_moles);

  /**
   * @brief Set the assemblage equilibrium tolerance.
   */
  void set_equilibrium_tolerance(double new_equilibrium_tolerance);

  // Public getters for extra Assemblage functions

  /**
   * @brief Retrieves molar fraction array.
   */
  Eigen::ArrayXd get_molar_fractions() const;

  /**
   * @brief Retrieves n_i * V_i.
   *
   * Uses a cached value if available, or calls
   * `Assemblage::compute_volume_fractions()` and caches the result.
   *
   * @note Use `Assemblage::reset_cache()` to clear cached values.
   *
   * @return Molar fractions * volumes of phases in assemblage.
   */
  Eigen::ArrayXd get_volume_fractions() const;

  /**
   * @brief Retrieves a list of the number of endmembers in each phase.
   */
  std::vector<int> get_endmembers_per_phase() const;

  /**
   * @brief Retrieves the phase at specified index.
   *
   * Returns a shared_ptr<Material>, use a dynamic cast to access behaviour
   * specific to derived classes (Mineral, Solution, etc.).
   */
  std::shared_ptr<Material> get_phase(std::size_t index) const;

  /**
   * @brief Retrieves the phase at specified index.
   *
   * Returns a shared_ptr<T> where T is the phase type.
   * Usage example: assemblage.get_phase<Solution>(0);
   */
  template <typename T>
  std::shared_ptr<T> get_phase(std::size_t index) const {
    return std::dynamic_pointer_cast<T>(this->phases.at(index));
  }

  /**
   * @brief Get number of phases in the assemblage.
   */
  Eigen::Index get_n_phases() const;

  /**
   * @brief Get n_moles (Used to convert mole fractions / absolute phase amounts).
   */
  double get_n_moles() const;

  /**
   * @brief Returns the equilibrium tolerance (J/reaction).
   *
   * @note Use `Assemblage::set_equilibtrium_tolerance' to change.
   * The default value of 1.0e-3 is rest on calls to `Assemblage::reset_cache()'.
   */
  double get_equilibrium_tolerance() const;

  /**
   * @brief Returns the reaction affinities vector.
   */
  Eigen::VectorXd get_reaction_affinities() const;

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
  types::FormulaMap compute_formula() const override;

  // Overrides of pure virtual functions from CompositeMaterial
  Eigen::Index compute_n_endmembers() const override;
  void setup_endmember_names() const override;
  void setup_endmember_formulae() const override;
  Eigen::ArrayXd compute_partial_gibbs() const override;

  // Additional Assemblage compute functions
  Eigen::ArrayXd compute_volume_fractions() const;
  Eigen::VectorXd compute_reaction_affinities() const;

 private:

  // Vector of pointers to phases in the assemblage
  std::vector<std::shared_ptr<Material>> phases;

  // Unique pointer to averaging scheme
  std::shared_ptr<averaging::AveragingScheme> averaging_scheme;

  Eigen::ArrayXd molar_fractions;

  // Cached properties (can be reset)
  // 1D Eigen arrays (could be made Vectors)
  mutable std::optional<Eigen::ArrayXd> volume_fractions;
  mutable std::optional<Eigen::VectorXd> reaction_affinities;

  // Stored (cached) properties not reset
  mutable std::optional<std::vector<int>> endmembers_per_phase;

  mutable std::optional<double> n_moles;
  double equilibrium_tolerance = 1.0e-3; // J/reaction

  // Compute / setup functions for Composite properties
  void setup_endmember_properties() const;

  // Setters for cached properties
  void set_endmembers_per_phase(std::vector<int> v) const;

  /**
   * @brief Convert from mass fractions to molar fractions.
   */
  Eigen::ArrayXd convert_mass_to_molar_fractions(const Eigen::ArrayXd& mass_fractions) const;

};

} // namespace burnman

#endif // BURNMAN_CORE_COMPOSITE_HPP_INCLUDED
