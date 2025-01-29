/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED
#define BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED

#include <string>
#include <optional>

/**
 * Base class for materials.

  TODO
  Python doc:

    Base class for all materials. The main functionality is unroll() which
    returns a list of objects of type :class:`~burnman.Mineral` and their molar
    fractions. This class is available as ``burnman.Material``.

    The user needs to call set_method() (once in the beginning) and set_state()
    before querying the material with unroll() or
 */ 
class Material {

 public:
  // Virtual destructor may be needed for cleanup in derived classes
  virtual ~Material() = default;

  /**
   * @brief Clears the property cache
   *
   * Resets cached material properties to empty `std::optional` so that
   * they are freshly computed on next query.
   */
  void reset();

  // Public getter functions
  /**
   * @brief Returns current pressure
   *
   * @note: Use Material::set_state() to set pressure
   *
   * @return Pressure in [Pa]
   */
  double get_pressure() const;

  /**
   * @brief Returns current temperature
   * 
   * @note: Use Material::set_state() to set temperature
   *
   * @return Temperature in [K]
   */
  double get_temperature() const;

  /**
   * @brief Retrieves the molar internal energy of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_internal_energy()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Internal energy in [J/mol].
   */
  double get_molar_internal_energy() const;

  /**
   * @brief Retrieves the molar Gibbs free energy of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_gibbs()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Gibbs free energy in [J/mol].
   */
  double get_molar_gibbs() const;

  /**
   * @brief Retrieves the molar Helmholtz free energy of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_helmholtz()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Helmholtz free energy in [J/mol].
   */
  double get_molar_helmholtz() const;

  /**
   * @brief Retrieves the molar mass of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_mass()` and caches the result.
   *
   * @note `Use Material::reset()` to clear cached values.
   *
   * @return Molar mass in [kg/mol].
   */
  double get_molar_mass() const;

    /**
   * @brief Retrieves the molar volume of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_volume()` and caches the result.
   *
   * @note `Use Material::reset()` to clear cached values.
   *
   * @return Molar volume in [m^3/mol].
   */
  double get_molar_volume() const;

  /**
   * @brief Retrieves the density of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_density()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Density in [kg/m^3].
   */
  double get_density() const;

  /**
   * @brief Retrieves the molar entropy of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_entropy()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Entropy in [J/K/mol].
   */
  double get_molar_entropy() const;

  /**
   * @brief Retrieves the molar enthalpy of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_enthalpy()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Enthalpy in [J/mol].
   */
  double get_molar_enthalpy() const;

  /**
   * @brief Retrieves the isothermal bulk modulus of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_isothermal_bulk_modulus_reuss()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Isothermal bulk modulus in [Pa].
   */
  double get_isothermal_bulk_modulus_reuss() const;

  /**
   * @brief Retrieves the isentropic bulk modulus of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_isentropic_bulk_modulus_reuss()`
   * and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Isentropic bulk modulus in [Pa].
   */
  double get_isentropic_bulk_modulus_reuss() const;
  
  /**
   * @brief Retrieves the isothermal compressibility of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_isothermal_compressibility_reuss()`
   * and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Isothermal compressibility in [1/Pa].
   */
  double get_isothermal_compressibility_reuss() const;

  /**
   * @brief Retrieves the isentropic compressibility of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_isentropic_compressibility_reuss()`
   * and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Adiabatic compressibility in [1/Pa].
   */
  double get_isentropic_compressibility_reuss() const;

  /**
   * @brief Retrieves the shear modulus of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_shear_modulus()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Shear modulus in [Pa].
   */
  double get_shear_modulus() const;

  /**
   * @brief Retrieves the P wave velocity of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_p_wave_velocity()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return P wave velocity in [m/s].
   */
  double get_p_wave_velocity() const;

  /**
   * @brief Retrieves the bulk sound velocity of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_bulk_sound_velocity()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Bulk sound velocity in [m/s].
   */
  double get_bulk_sound_velocity() const;

  /**
   * @brief Retrieves shear wave velocity the  of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_shear_wave_velocity()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Shear wave velocity in [m/s].
   */
  double get_shear_wave_velocity() const;

  /**
   * @brief Retrieves the grueneisen parameter of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_grueneisen_parameter()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Grueneisen parameter [unitless].
   */
  double get_grueneisen_parameter() const;

  /**
   * @brief Retrieves the thermal expansion coefficient of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_thermal_expansivity()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Thermal expansivity in [1/K].
   */
  double get_thermal_expansivity() const;

  /**
   * @brief Retrieves the molar heat capacity at constant volume of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_heat_capacity_v()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Isochoric heat capacity in [J/K/mol].
   */
  double get_molar_heat_capacity_v() const;

  /**
   * @brief Retrieves molar heat capacity at constant pressure of the mineral.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_molar_heat_capacity_p()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return Isobaric heat capacity in [J/K/mol].
   */
  double get_molar_heat_capacity_p() const;

  /**
   * @brief Retrieves change in temperature with pressure at constant entropy.
   *
   * Uses a cached value if available, or calls
   * `Material::compute_isentropic_thermal_gradient()` and caches the result.
   *
   * @note Use `Material::reset()` to clear cached values.
   *
   * @return dTdP at constant entropy in [Pa/K].
   */
  double get_isentropic_thermal_gradient() const;


 protected:

  // protected compute functions to override in derived classes
  /**
   * @brief Computes the molar internal energy of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Internal energy in [J/mol].
   * @throws NotImplemntedError if default implementation called.
   */
  virtual double compute_molar_internal_energy() const;

  /**
   * @brief Computes the molar Gibbs free energy of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Gibbs free energy in [J/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_gibbs() const;

  /**
   * @brief Computes the molar Helmholtz free energy of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Helmholtz free energy in [J/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_helmholtz() const;

  /**
   * @brief Computes the molar mass of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Molar mass in [kg/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_mass() const;

    /**
   * @brief Computes the molar volume of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Molar volume in [m^3/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_volume() const;

  /**
   * @brief Computes the density of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Density in [kg/m^3].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_density() const;

  /**
   * @brief Computes the molar entropy of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Entropy in [J/K/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_entropy() const;

  /**
   * @brief Computes the molar enthalpy of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Enthalpy in [J/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_enthalpy() const;

  /**
   * @brief Computes the isothermal bulk modulus of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Isothermal bulk modulus in [Pa].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_isothermal_bulk_modulus_reuss() const;

  /**
   * @brief Computes the isentropic bulk modulus of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Isentropic bulk modulus in [Pa].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_isentropic_bulk_modulus_reuss() const;
  
  /**
   * @brief Computes the isothermal compressibility of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Isothermal compressibility in [1/Pa].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_isothermal_compressibility_reuss() const;

  /**
   * @brief Computes the isentropic compressibility of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Adiabatic compressibility in [1/Pa].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_isentropic_compressibility_reuss() const;

  /**
   * @brief Computes the shear modulus of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Shear modulus in [Pa].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_shear_modulus() const;

  /**
   * @brief Computes the P wave velocity of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return P wave velocity in [m/s].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_p_wave_velocity() const;

  /**
   * @brief Computes the bulk sound velocity of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Bulk sound velocity in [m/s].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_bulk_sound_velocity() const;

  /**
   * @brief Computes shear wave velocity the  of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Shear wave velocity in [m/s].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_shear_wave_velocity() const;

  /**
   * @brief Computes the grueneisen parameter of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Grueneisen parameter [unitless].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_grueneisen_parameter() const;

  /**
   * @brief Computes the thermal expansion coefficient of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Thermal expansivity in [1/K].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_thermal_expansivity() const;

  /**
   * @brief Computes the molar heat capacity at constant volume of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Isochoric heat capacity in [J/K/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_heat_capacity_v() const;

  /**
   * @brief Computes molar heat capacity at constant pressure of the mineral.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return Isobaric heat capacity in [J/K/mol].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_molar_heat_capacity_p() const;

  /**
   * @brief Computes change in temperature with pressure at constant entropy.
   *
   * @note Default implementation throws NotImplementedError.
   *       Derived classes should override method.
   *
   * @return dTdP at constant entropy in [Pa/K].
   * @throws NotImplementedError if default implementation called.
   */
  virtual double compute_isentropic_thermal_gradient() const;


 private:
  // std::optional used for caching
  mutable std::optional<double> pressure;
  mutable std::optional<double> temperature;
  mutable std::optional<double> molar_internal_energy;
  mutable std::optional<double> molar_gibbs;
  mutable std::optional<double> molar_helmholtz;
  mutable std::optional<double> molar_mass;
  mutable std::optional<double> molar_volume;
  mutable std::optional<double> density;
  mutable std::optional<double> molar_entropy;
  mutable std::optional<double> molar_enthalpy;
  mutable std::optional<double> isothermal_bulk_modulus_reuss;
  mutable std::optional<double> isentropic_bulk_modulus_reuss;
  mutable std::optional<double> isothermal_compressibility_reuss;
  mutable std::optional<double> isentropic_compressibility_reuss;
  mutable std::optional<double> shear_modulus;
  mutable std::optional<double> p_wave_velocity;
  mutable std::optional<double> bulk_sound_velocity;
  mutable std::optional<double> shear_wave_velocity;
  mutable std::optional<double> grueneisen_parameter;
  mutable std::optional<double> thermal_expansivity;
  mutable std::optional<double> molar_heat_capacity_v;
  mutable std::optional<double> molar_heat_capacity_p;
  mutable std::optional<double> isentropic_thermal_gradient;

  // Private functions for throwing exceptions
  /**
   * @brief Helper function to throw NotImplementedError
   *
   * Use __func__ for method name
   */
  [[noreturn]] void throw_not_implemented_error(const std::string& method) const;
  
  /**
   * @brief Helper function to get name of class (incl. derived)
   */
  std::string get_class_name() const;

};

#endif // BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED