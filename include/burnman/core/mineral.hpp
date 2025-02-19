/*
  TODO: Copyright Notice!
*/
#ifndef BURNMAN_CORE_MINERAL_HPP_INCLUDED
#define BURNMAN_CORE_MINERAL_HPP_INCLUDED

#include <string>
#include <optional>
#include "burnman/core/material.hpp"
#include "burnman/core/equation_of_state.hpp"
#include "burnman/utils/eos.hpp"

/**
 * Base class for materials.

  TODO
  Python doc:

    This is the base class for all minerals. States of the mineral
    can only be queried after setting the pressure and temperature
    using set_state(). The method for computing properties of
    the material is set using set_method(). This is done during
    initialisation if the param 'equation_of_state' has been defined.
    The method can be overridden later by the user.

    This class is available as ``burnman.Mineral``.

    If deriving from this class, set the properties in self.params
    to the desired values. For more complicated materials you
    can overwrite set_state(), change the params and then call
    set_state() from this class.

    All the material parameters are expected to be in plain SI units.  This
    means that the elastic moduli should be in Pascals and NOT Gigapascals,
    and the Debye temperature should be in K not C.  Additionally, the
    reference volume should be in m^3/(mol molecule) and not in unit cell
    volume and 'n' should be the number of atoms per molecule.  Frequently in
    the literature the reference volume is given in Angstrom^3 per unit cell.
    To convert this to m^3/(mol of molecule) you should multiply by 10^(-30) *
    N_a / Z, where N_a is Avogadro's number and Z is the number of formula units per
    unit cell. You can look up Z in many places, including www.mindat.org

  TODO
  Funcs:
    __init__(self, params=None, property_modifiers=None):
    set_method
    to_string
    debug_print
    unroll
    set_state
    _molar_volume_unmodified
    formula
 */ 
class Mineral : public Material{
 public:

  // Parameter object
  MineralParams params;

  // EOS Class for equation of state to use
  EquationOfState eos_method;

  // Override set_method etc.
  void set_state(double new_pressure, double new_temperature) override;

 protected:

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

 private:
  // Excesses struct for property modifiers
  // Values initialise to 0
  excesses::Excesses property_modifier_excesses;
  excesses::ExcessParamVector property_modifier_params;

};


#endif // BURNMAN_CORE_MINERAL_HPP_INCLUDED