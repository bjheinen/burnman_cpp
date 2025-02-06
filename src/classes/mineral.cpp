/*
  TODO: Copyright Notice!
*/
#include "../../include/burnman/classes/mineral.hpp"

// store _property_modifiers['etc'] under property_modifier_excesses.etc

// EOS properties - in P,T form

double Mineral::compute_molar_gibbs() const {
  // calls method.gibbs_free_energy()

  get_pressure();
  get_temperature();
  

  // adds ._property_modifiers['G']
  property_modifier_excesses.G;
}

// should probably cache this as well...
double Mineral::compute_molar_volume_unmodified() const{
  return method.volume(get_pressure(), get_temperature(), params);
}

double Mineral::compute_molar_volume() const {
  return compute_molar_volume_unmodified() + property_modifier_excesses.dGdP;
}

double Mineral::compute_molar_enthalpy() const {
  return method.entropy(
    get_pressure(),
    get_temperature(),
    compute_molar_volume_unmodified(),
    params
  ) - property_modifier_excesses.dGdT;
}




// TODO: 
//       implement set_method so know how to call EOS funcs
//         should EOS be functional to reduce overheads??
//       implement params as a struct
//       implement property_modifiers as a struct (or with params?)




// Mineral parameter properties