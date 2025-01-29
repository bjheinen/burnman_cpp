/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED
#define BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED

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

};

#endif // BURNMAN_CLASSES_MATERIAL_HPP_INCLUDED